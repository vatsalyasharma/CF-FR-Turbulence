#include "NavierStokes/TransitionTransport.hh"
#include "NavierStokes/TransitionArithmetic.hh"
#include "NavierStokes3DGReKLogOSSTPuvt.hh"
#include "Common/BadValueException.hh"
#include <typeinfo>
#include "GReKO.hh"
#include "NavierStokes3DGReKLogOPuvt.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/ShouldNotBeHereException.hh"
#include "Common/NotImplementedException.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::KOmega;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GReKO {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NavierStokes3DGReKLogOPuvt, DiffusiveVarSet,
			    GReKOModule, 2>
ns3DGReKLogOPuvtProvider("NavierStokes3DGReKLogOPvt");

//////////////////////////////////////////////////////////////////////////////

NavierStokes3DGReKLogOPuvt::NavierStokes3DGReKLogOPuvt
(const std::string& name, SafePtr<PhysicalModelImpl> model) :
  NavierStokesKLogOmegaPvt<NavierStokesKLogOmegaVarSet<NavierStokes3DVarSet, 0> >(name, model),
  _unperturbedFluxGa(),
  _unperturbedFluxRe()
{
  const CFuint nbTurbEquations = _eulerModel->getNbScalarVars(0);
  
  vector<std::string> names(5 + nbTurbEquations);
  names[0] = "p";
  names[1] = "u";
  names[2] = "v";
  names[3] = "w";
  names[4] = "T";
  
  // Names for turbulent variables
  names[5] = "K";
  names[6] = "Omega";
  
  // Names for transition onset variables
  names[7] = "Ga";
  names[8] = "Re";
  
  setVarNames(names);
  setModelCoefficients();
}

CFreal NavierStokes3DGReKLogOPuvt::getLaminarDynViscosityFromGradientVars(const RealVector& state)
{
  return transition_detail::laminarViscosity(getModel(), state[0], state[4],
      _eulerModel->getPressRef(), _eulerModel->getTempRef());
}

lm2009::Inputs NavierStokes3DGReKLogOPuvt::getLM2009Inputs(
    const RealVector& state, const std::vector<RealVector*>& gradients)
{
  lm2009::Inputs p;p.dim=3;
  if(typeid(*this)==typeid(NavierStokes3DGReKLogOSSTPuvt))p.model=lm2009::Inputs::SST;
  else if(typeid(*this)==typeid(NavierStokes3DGReKLogOBSLPuvt))p.model=lm2009::Inputs::BSL;
  else if(typeid(*this)==typeid(NavierStokes3DGReKLogOPuvt))p.model=lm2009::Inputs::Base;
  else throw Common::BadValueException(FromHere(), "Unsupported LM2009 constitutive model");
  for(unsigned e=0;e<9;++e){p.state[e]=state[e];for(unsigned d=0;d<3;++d)p.grad[e][d]=(*gradients[e])[d];}
  p.distance=_wallDistance;p.R=_eulerModel->getR();p.isentropicGamma=_eulerModel->getGamma();
  p.mu=getLaminarDynViscosityFromGradientVars(state);
  const DynamicViscosityLaw& law=getModel().getDynamicViscosityLaw();
  const CFreal T=state[4],scale=_eulerModel->getTempRef();
  if(typeid(law)==typeid(SutherlandDynViscosity)){
    const CFreal shift=static_cast<const SutherlandDynViscosity&>(law).getSutherlandConstant();
    if(scale*T>1.e-10)p.muT=p.mu*(1.5/T-scale/(scale*T+shift));
  } else if(typeid(law)==typeid(SimplerSutherlandDynViscosity)){
    const CFreal shift=static_cast<const SimplerSutherlandDynViscosity&>(law).getTemperatureShift();
    if(scale*T>1.e-10)p.muT=p.mu*(1.5/T-scale/(scale*T+shift));
  } else if(typeid(law)==typeid(FixedKinViscosity)){p.muP=p.mu/state[0];p.muT=-p.mu/T;}
  // Fixed/constant laws have zero pressure/temperature partials; unsupported
  // laws have already been rejected by checked laminar transport above.
  p.sigmaK1=_sigmaK1;p.sigmaK2=_sigmaK2;p.sigmaW1=_sigmaOmega1;p.sigmaW2=_sigmaOmega2;
  p.beta1=_beta1;p.beta2=_beta2;p.betaStar=_betaStar;p.gamma1=_gamma1;p.gamma2=_gamma2;p.coeffTau=getModel().getCoeffTau();
  return p;
}

CFreal NavierStokes3DGReKLogOPuvt::getDensity(const RealVector& state)
{
  transition_detail::validateThermodynamics(state[0], state[4], _eulerModel->getR(), _eulerModel->getGamma());
  return state[0]/(_eulerModel->getR()*state[4]);
}

CFreal NavierStokes3DGReKLogOPuvt::getRhoKOverMuT(const RealVector& state, const std::vector<RealVector*>& gradients)
{
  // Dispatch only the implemented constitutive laws; never infer a zero-mut
  // limit from a floor or a synthetic positive-k state.
  if (typeid(*this) == typeid(NavierStokes3DGReKLogOSSTPuvt))
    return static_cast<NavierStokes3DGReKLogOSSTPuvt*>(this)->getRhoKOverMuT(state, gradients);
  if (typeid(*this) == typeid(NavierStokes3DGReKLogOPuvt) || typeid(*this) == typeid(NavierStokes3DGReKLogOBSLPuvt))
    return std::exp(state[6]);
  throw Common::BadValueException(FromHere(), "Transition source viscosity ratio is not implemented for this diffusive variable set");
}

//////////////////////////////////////////////////////////////////////////////

NavierStokes3DGReKLogOPuvt::~NavierStokes3DGReKLogOPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

RealVector& NavierStokes3DGReKLogOPuvt::getFlux(const RealVector& values,
                                              const vector<RealVector*>& gradients,
                                              const RealVector& normal,
                                              const CFreal& radius)
{
  // here nsData is filled in
  getLaminarDynViscosityFromGradientVars(values);
  computeTransportProperties(values, gradients, normal);
  computeStressTensor(values, gradients, radius);
  setGradientState(values);
  
  const RealVector& gradT     = *gradients[4];
  const RealVector& gradK     = *gradients[5];
  const RealVector& gradOmega = *gradients[6];
  const RealVector& gradGa    = *gradients[7];
  const RealVector& gradRe    = *gradients[8];

  const CFreal avU = values[1];
  const CFreal avV = values[2];
  const CFreal avW = values[3];

  RealVector& nsData = getModel().getPhysicalData();
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal nz = normal[ZZ];
  
  const CFreal qFlux = -getModel().getCoeffQ()*nsData[NSTerm::LAMBDA]*
    MathFunctions::innerProd(*gradients[_TID], normal);
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  _flux.slice(_uID, dim) = _tau*normal;
  _flux[_TID] = MathFunctions::innerProd(_flux.slice(_uID,dim), _gradState.slice(_uID,dim)) - qFlux;
 
//  const CFreal tauXX = _tau(XX, XX);
//  const CFreal tauXY = _tau(XX, YY);
//  const CFreal tauYY = _tau(YY, YY);
//  _flux[1] = tauXX*nx + tauXY*ny;
//  _flux[2] = tauXY*nx + tauYY*ny;
//  _flux[3] = (tauXX*avU + tauXY*avV)*nx + (tauXY*avU + tauYY*avV)*ny - qFlux;
  
  //diff flux corresponding to K
  computeBlendingCoefFromGradientVars(values,gradK, gradOmega);
  CFreal coeffTauMu1 = getModel().getCoeffTau()*nsData[NSTurbTerm::MU];
  CFreal coeffTauMu2 = getModel().getCoeffTau()*nsData[NSTurbTerm::MUT]*getSigmaK();

  _flux[5] = (coeffTauMu1 + coeffTauMu2)*(gradK[XX]*nx + gradK[YY]*ny + gradK[ZZ]*nz);
  //diff flux corresponding to Omega
  coeffTauMu1 = getModel().getCoeffTau()*nsData[NSTurbTerm::MU];
  coeffTauMu2 = getModel().getCoeffTau()*nsData[NSTurbTerm::MUT]*getSigmaOmega();

  _flux[6] = (coeffTauMu1 + coeffTauMu2)*(gradOmega[XX]*nx + gradOmega[YY]*ny + gradOmega[ZZ]*nz);
  
  //diff flux corresponding to Gamma
  const CFreal sigmaf= 1.0;
  coeffTauMu1 = getModel().getCoeffTau()*nsData[NSTurbTerm::MU];
  coeffTauMu2 = getModel().getCoeffTau()*nsData[NSTurbTerm::MUT]/sigmaf;

  _flux[7] = (coeffTauMu1 + coeffTauMu2)*(gradGa[XX]*nx + gradGa[YY]*ny + gradGa[ZZ]*nz);
  
  //diff flux corresponding to Re
  const CFreal sigmatheta= 2.0;
  coeffTauMu1 = getModel().getCoeffTau()*nsData[NSTurbTerm::MU];
  coeffTauMu2 = getModel().getCoeffTau()*nsData[NSTurbTerm::MUT];

  _flux[8] = sigmatheta*(coeffTauMu1 + coeffTauMu2)*(gradRe[XX]*nx + gradRe[YY]*ny + gradRe[ZZ]*nz);
  //jacobian contribution
  if(_isPerturb)
  {
    if(_iPerturbVar != 5){
      _flux[5] = _unperturbedFluxK;
    }
    if(_iPerturbVar != 6){
      _flux[6] = _unperturbedFluxOmega;
    }
    if(_iPerturbVar != 7){
      _flux[7] = _unperturbedFluxGa;
    }
    if(_iPerturbVar != 8){
      _flux[8] = _unperturbedFluxRe;
    }
  }
  else
  {
    _unperturbedFluxK = _flux[5];
    _unperturbedFluxOmega = _flux[6];
    _unperturbedFluxGa    = _flux[7];
    _unperturbedFluxRe    = _flux[8];
  }
 //cout << _flux << endl;
  return _flux;
}

//////////////////////////////////////////////////////////////////////////////


    } // namespace GReKO

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
