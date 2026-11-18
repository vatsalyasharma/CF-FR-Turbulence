#include "NavierStokes/TransitionTransport.hh"
#include "NavierStokes/TransitionArithmetic.hh"
#include "NavierStokes2DGReKLogOSSTPuvt.hh"
#include "Common/BadValueException.hh"
#include <typeinfo>
#include "GReKO.hh"
#include "NavierStokes2DGReKLogOPuvt.hh"
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

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GReKO {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NavierStokes2DGReKLogOPuvt, DiffusiveVarSet,
			    GReKOModule, 2>
ns2DGReKLogOPuvtProvider("NavierStokes2DGReKLogOPuvt");

//////////////////////////////////////////////////////////////////////////////

NavierStokes2DGReKLogOPuvt::NavierStokes2DGReKLogOPuvt
(const std::string& name, SafePtr<PhysicalModelImpl> model) :
  NavierStokesKLogOmegaPvt<NavierStokesKLogOmegaVarSet<NavierStokes2DVarSet, 0> >(name, model),
  _unperturbedFluxGa(),
  _unperturbedFluxRe()
{
  const CFuint nbTurbEquations = _eulerModel->getNbScalarVars(0);
  
  vector<std::string> names(4 + nbTurbEquations);
  names[0] = "p";
  names[1] = "u";
  names[2] = "v";
  names[3] = "T";
  
  // Names for turbulent variables
  names[4] = "K";
  names[5] = "Omega";
  
  // Names for transition onset variables
  names[6] = "Ga";
  names[7] = "Re";
  
  setVarNames(names);
  setModelCoefficients();
}

CFreal NavierStokes2DGReKLogOPuvt::getLaminarDynViscosityFromGradientVars(const RealVector& state)
{
  return transition_detail::laminarViscosity(getModel(), state[0], state[3],
      _eulerModel->getPressRef(), _eulerModel->getTempRef());
}

lm2009::Inputs NavierStokes2DGReKLogOPuvt::getLM2009Inputs(
    const RealVector& state, const std::vector<RealVector*>& gradients)
{
  lm2009::Inputs p;p.dim=2;
  if(typeid(*this)==typeid(NavierStokes2DGReKLogOSSTPuvt))p.model=lm2009::Inputs::SST;
  else if(typeid(*this)==typeid(NavierStokes2DGReKLogOBSLPuvt))p.model=lm2009::Inputs::BSL;
  else if(typeid(*this)==typeid(NavierStokes2DGReKLogOPuvt))p.model=lm2009::Inputs::Base;
  else throw Common::BadValueException(FromHere(), "Unsupported LM2009 constitutive model");
  for(unsigned e=0;e<8;++e){p.state[e]=state[e];for(unsigned d=0;d<2;++d)p.grad[e][d]=(*gradients[e])[d];}
  p.distance=_wallDistance;p.R=_eulerModel->getR();p.isentropicGamma=_eulerModel->getGamma();
  p.mu=getLaminarDynViscosityFromGradientVars(state);
  const DynamicViscosityLaw& law=getModel().getDynamicViscosityLaw();
  const CFreal T=state[3],scale=_eulerModel->getTempRef();
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

CFreal NavierStokes2DGReKLogOPuvt::getDensity(const RealVector& state)
{
  transition_detail::validateThermodynamics(state[0], state[3], _eulerModel->getR(), _eulerModel->getGamma());
  return state[0]/(_eulerModel->getR()*state[3]);
}

CFreal NavierStokes2DGReKLogOPuvt::getRhoKOverMuT(const RealVector& state, const std::vector<RealVector*>& gradients)
{
  // Dispatch only the implemented constitutive laws; never infer a zero-mut
  // limit from a floor or a synthetic positive-k state.
  if (typeid(*this) == typeid(NavierStokes2DGReKLogOSSTPuvt))
    return static_cast<NavierStokes2DGReKLogOSSTPuvt*>(this)->getRhoKOverMuT(state, gradients);
  if (typeid(*this) == typeid(NavierStokes2DGReKLogOPuvt) || typeid(*this) == typeid(NavierStokes2DGReKLogOBSLPuvt))
    return std::exp(state[5]);
  throw Common::BadValueException(FromHere(), "Transition source viscosity ratio is not implemented for this diffusive variable set");
}

omega_detail::Partial NavierStokes2DGReKLogOPuvt::getLogOmegaSourcePartial(
    const RealVector& state, const std::vector<RealVector*>& gradients,
    const CFreal shearProduction, const CFreal dilatation,
    const omega_detail::Partial destructionFactor, const bool limitProduction)
{
  omega_detail::Inputs p;
  if(typeid(*this)==typeid(NavierStokes2DGReKLogOSSTPuvt)) p.model=omega_detail::Inputs::SST;
  else if(typeid(*this)==typeid(NavierStokes2DGReKLogOBSLPuvt)) p.model=omega_detail::Inputs::BSL;
  else if(typeid(*this)==typeid(NavierStokes2DGReKLogOPuvt)) p.model=omega_detail::Inputs::Base;
  else throw Common::BadValueException(FromHere(), "GR log-omega partials require a supported base/BSL/SST variable set");
  p.rho=getDensity(state);p.mu=getLaminarDynViscosityFromGradientVars(state);
  p.k=std::max(state[4],0.);p.omega=std::exp(state[5]);p.distance=_wallDistance;
  p.soundSpeedSquared=getSqSoundSpeed(state);
  CFreal velocityGradients[3][3]={};
  for(CFuint e=0;e<2;++e)for(CFuint d=0;d<2;++d)velocityGradients[e][d]=(*gradients[e+1])[d];
  p.strain=transition_detail::velocityGradientNorm(velocityGradients,2,false);
  p.gradKDotGradLogOmega=(*gradients[4])[XX]*(*gradients[5])[XX]+(*gradients[4])[YY]*(*gradients[5])[YY];
  p.gradLogOmegaSquared=(*gradients[5])[XX]*(*gradients[5])[XX]+(*gradients[5])[YY]*(*gradients[5])[YY];
  p.sigma1=_sigmaOmega1;p.sigma2=_sigmaOmega2;p.beta1=_beta1;p.beta2=_beta2;
  p.betaStar=_betaStar;p.gamma1=_gamma1;p.gamma2=_gamma2;p.coeffTau=getModel().getCoeffTau();
  p.shearProduction=shearProduction;p.dilatation=dilatation;
  p.destructionFactor=destructionFactor;p.limitProduction=limitProduction;
  return omega_detail::source(p);
}

//////////////////////////////////////////////////////////////////////////////

NavierStokes2DGReKLogOPuvt::~NavierStokes2DGReKLogOPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

RealVector& NavierStokes2DGReKLogOPuvt::getFlux(const RealVector& values,
                                              const vector<RealVector*>& gradients,
                                              const RealVector& normal,
                                              const CFreal& radius)
{
  // here nsData is filled in
  getLaminarDynViscosityFromGradientVars(values);
  computeTransportProperties(values, gradients, normal);
  computeStressTensor(values, gradients, radius);
  
  const RealVector& gradT     = *gradients[3];
  const RealVector& gradK     = *gradients[4];
  const RealVector& gradOmega = *gradients[5];
  const RealVector& gradGa    = *gradients[6];
  const RealVector& gradRe    = *gradients[7];

  const CFreal avU = values[1];
  const CFreal avV = values[2];

  RealVector& nsData = getModel().getPhysicalData();
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];

  const CFreal qFlux = -getModel().getCoeffQ()*nsData[NSTurbTerm::LAMBDA]*
    (gradT[XX]*nx + gradT[YY]*ny);
  
  const CFreal tauXX = _tau(XX, XX);
  const CFreal tauXY = _tau(XX, YY);
  const CFreal tauYY = _tau(YY, YY);
  _flux[1] = tauXX*nx + tauXY*ny;
  _flux[2] = tauXY*nx + tauYY*ny;
  _flux[3] = (tauXX*avU + tauXY*avV)*nx + (tauXY*avU + tauYY*avV)*ny - qFlux;
  
  //diff flux corresponding to K
  computeBlendingCoefFromGradientVars(values,gradK, gradOmega);
  CFreal coeffTauMu1 = getModel().getCoeffTau()*nsData[NSTurbTerm::MU];
  CFreal coeffTauMu2 = getModel().getCoeffTau()*nsData[NSTurbTerm::MUT]*getSigmaK();

  _flux[4] = (coeffTauMu1 + coeffTauMu2)*(gradK[XX]*nx + gradK[YY]*ny);
  //diff flux corresponding to Omega
  coeffTauMu1 = getModel().getCoeffTau()*nsData[NSTurbTerm::MU];
  coeffTauMu2 = getModel().getCoeffTau()*nsData[NSTurbTerm::MUT]*getSigmaOmega();

  _flux[5] = (coeffTauMu1 + coeffTauMu2)*(gradOmega[XX]*nx + gradOmega[YY]*ny);
  
  //diff flux corresponding to Gamma
  const CFreal sigmaf= 1.0;
  coeffTauMu1 = getModel().getCoeffTau()*nsData[NSTurbTerm::MU];
  coeffTauMu2 = getModel().getCoeffTau()*nsData[NSTurbTerm::MUT]/sigmaf;

  _flux[6] = (coeffTauMu1 + coeffTauMu2)*(gradGa[XX]*nx + gradGa[YY]*ny);
  
  //diff flux corresponding to Re
  const CFreal sigmatheta= 2.0;
  coeffTauMu1 = getModel().getCoeffTau()*nsData[NSTurbTerm::MU];
  coeffTauMu2 = getModel().getCoeffTau()*nsData[NSTurbTerm::MUT];

  _flux[7] = sigmatheta*(coeffTauMu1 + coeffTauMu2)*(gradRe[XX]*nx + gradRe[YY]*ny);
  //jacobian contribution
  if(_isPerturb)
  {
    if(_iPerturbVar != 4){
      _flux[4] = _unperturbedFluxK;
    }
    if(_iPerturbVar != 5){
      _flux[5] = _unperturbedFluxOmega;
    }
    if(_iPerturbVar != 6){
      _flux[6] = _unperturbedFluxGa;
    }
    if(_iPerturbVar != 7){
      _flux[7] = _unperturbedFluxRe;
    }
  }
  else
  {
    _unperturbedFluxK = _flux[4];
    _unperturbedFluxOmega = _flux[5];
    _unperturbedFluxGa    = _flux[6];
    _unperturbedFluxRe    = _flux[7];
  }
 //cout << _flux << endl;
  return _flux;
}

//////////////////////////////////////////////////////////////////////////////


    } // namespace GReKO

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
