#include "NavierStokes/TransitionTransport.hh"
#include "NavierStokes/TransitionArithmetic.hh"
#include "NavierStokes3DGammaAlphaSSTPuvt.hh"
#include "Common/BadValueException.hh"
#include <typeinfo>
#include "GammaAlpha.hh"
#include "NavierStokes3DGammaAlphaPuvt.hh"
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

    namespace GammaAlpha {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NavierStokes3DGammaAlphaPuvt, DiffusiveVarSet,
			    GammaAlphaModule, 2>
ns3DGammaAlphaPuvtProvider("NavierStokes3DGammaAlphaPvt");

//////////////////////////////////////////////////////////////////////////////

NavierStokes3DGammaAlphaPuvt::NavierStokes3DGammaAlphaPuvt
(const std::string& name, SafePtr<PhysicalModelImpl> model) :
  NavierStokesKLogOmegaPvt<NavierStokesKLogOmegaVarSet<NavierStokes3DVarSet, 0> >(name, model),
  _unperturbedFluxGa(),
  _unperturbedFluxAlpha()
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
  names[7] = "gamma";
  names[8] = "alpha";
  
  setVarNames(names);
  setModelCoefficients();
}

CFreal NavierStokes3DGammaAlphaPuvt::getLaminarDynViscosityFromGradientVars(const RealVector& state)
{
  return transition_detail::laminarViscosity(getModel(), state[0], state[4],
      _eulerModel->getPressRef(), _eulerModel->getTempRef());
}

CFreal NavierStokes3DGammaAlphaPuvt::getDensity(const RealVector& state)
{
  transition_detail::validateThermodynamics(state[0], state[4], _eulerModel->getR(), _eulerModel->getGamma());
  return state[0]/(_eulerModel->getR()*state[4]);
}

CFreal NavierStokes3DGammaAlphaPuvt::getRhoKOverMuT(const RealVector& state, const std::vector<RealVector*>& gradients)
{
  // Dispatch only the implemented constitutive laws; never infer a zero-mut
  // limit from a floor or a synthetic positive-k state.
  if (typeid(*this) == typeid(NavierStokes3DGammaAlphaSSTPuvt))
    return static_cast<NavierStokes3DGammaAlphaSSTPuvt*>(this)->getRhoKOverMuT(state, gradients);
  if (typeid(*this) == typeid(NavierStokes3DGammaAlphaBSLPuvt))
    return static_cast<NavierStokes3DGammaAlphaBSLPuvt*>(this)->getRhoKOverMuT(state, gradients);
  if (typeid(*this) == typeid(NavierStokes3DGammaAlphaPuvt))
    return std::exp(state[6]);
  throw Common::BadValueException(FromHere(), "Transition source viscosity ratio is not implemented for this diffusive variable set");
}

//////////////////////////////////////////////////////////////////////////////

NavierStokes3DGammaAlphaPuvt::~NavierStokes3DGammaAlphaPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

RealVector& NavierStokes3DGammaAlphaPuvt::getFlux(const RealVector& values,
                                              const vector<RealVector*>& gradients,
                                              const RealVector& normal,
                                              const CFreal& radius)
{
  // Check before transport/stress evaluation or writing the native flux buffer.
  const auto freestream = transition_detail::gammaAlphaFreestream(std::max(values[0],1.0e-5),
      _eulerModel->getR(), _eulerModel->getGamma(), _eulerModel->getMachInf(),
      _eulerModel->getVelInf(), _eulerModel->getStaticPressInf());
  // compute the flux in the base class
  //RealVector& flux = BASE::getFlux(values, gradients, normal, radius);
  //CFLog(INFO, "uID: " << _uID << ", vID: " << _vID << ", wID: " << _wID << ", TID: " << _TID << ", kID: " << _kID << "\n");
  setGradientState(values);
  getLaminarDynViscosityFromGradientVars(values);
  transition_detail::laminarViscosity(getModel(), values[0], freestream.temperature,
      _eulerModel->getPressRef(), _eulerModel->getTempRef());
  computeTransportProperties(values, gradients, normal);
  computeStressTensor(values, gradients, radius);
  
  const RealVector& nsData = getModel().getPhysicalData();
  const CFreal qFlux = -getModel().getCoeffQ()*nsData[NSTerm::LAMBDA]*
    MathFunctions::innerProd(*gradients[_TID], normal);
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  _flux.slice(_uID, dim) = _tau*normal;
  _flux[_TID] = MathFunctions::innerProd(_flux.slice(_uID,dim), _gradState.slice(_uID,dim)) - qFlux;
  
  //CFLog(INFO, "u v w: " << _gradState.slice(_uID,dim) << "\n");
  
  const RealVector& gradK     = *gradients[5];
  const RealVector& gradOmega = *gradients[6];
  const RealVector& gradGa    = *gradients[7];
  const RealVector& gradAlpha = *gradients[8];
  
  const CFreal p = max(values[0],1.0e-5);
  const CFreal gamma = min(max(values[7],0.01),0.99);
  
  //diff flux corresponding to K
  computeBlendingCoefFromGradientVars(values,gradK, gradOmega);
  
  const CFreal coeffTauMu1 = this->getModel().getCoeffTau()*nsData[NSTurbTerm::MU];
  const CFreal coeffTauMu2 = this->getModel().getCoeffTau()*nsData[NSTurbTerm::MUT]*getSigmaK();
  const CFreal coeffTauMu3 = this->getModel().getCoeffTau()*nsData[NSTurbTerm::MUT]*getSigmaOmega();
    
  _flux[5] = (coeffTauMu1 + coeffTauMu2)*MathFunctions::innerProd(gradK, normal);
  _flux[6] = (coeffTauMu1 + coeffTauMu3)*MathFunctions::innerProd(gradOmega, normal);
  
  //diff flux corresponding to Gamma
  //const CFreal sigmaf= 1.0;
  //const CFreal coeffTauMu4 = getModel().getCoeffTau()*nsData[NSTurbTerm::MUT]/sigmaf;
  
  const CFreal gammaIsentropic = _eulerModel->getGamma();
  const CFreal rhoInfLocal = freestream.density;
  const CFreal MInfLocal = freestream.mach;
  const CFreal uInfLocal = freestream.speed;
  const CFreal TInfLocal = freestream.temperature;
  
  // RealVector copy construction can retain a non-owning state view.
  // Allocate owned storage before changing the freestream temperature.
  RealVector copyState(values.size());
  copyState = values;
  copyState[4] = TInfLocal;
  
  const CFreal muInfLocal = getLaminarDynViscosityFromGradientVars(copyState);
  
  //const CFreal distance = std::max(_wallDistance, 1.0e-12);
  
  const CFreal fMuGamma = 1.0-exp(-256.0*(_wallDistance*uInfLocal*rhoInfLocal/muInfLocal)*(_wallDistance*uInfLocal*rhoInfLocal/muInfLocal));
  const CFreal fMMuGamma = (1.0+0.26*(gammaIsentropic-1.0)/2.0*MInfLocal*MInfLocal)*sqrt(1+0.38*pow(MInfLocal,0.6));
 //if (_wallDistance < 1.0e-4) CFLog(INFO, "wallDist: " << _wallDistance << "\n"); 
  //const CFreal gammaLim = min(max(0.01,gamma),0.99);
  const CFreal muGamma = 0.57*pow(-log(1-gamma),-5.0/6.0*(1.0-gamma))*fMuGamma*fMMuGamma*getModel().getCoeffTau()*nsData[NSTurbTerm::MU];

  _flux[7] = (muGamma)*MathFunctions::innerProd(gradGa, normal);
  
  //diff flux corresponding to alpha
  const CFreal sigmatheta= 2.0;
  const CFreal coeffTauMu4 = getModel().getCoeffTau()*nsData[NSTurbTerm::MUT];

  _flux[8] = sigmatheta*(coeffTauMu1 + coeffTauMu4)*MathFunctions::innerProd(gradAlpha, normal);
  
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
      _flux[8] = _unperturbedFluxAlpha;
    }
  }
  else
  {
    _unperturbedFluxK = _flux[5];
    _unperturbedFluxOmega = _flux[6];
    _unperturbedFluxGa    = _flux[7];
    _unperturbedFluxAlpha = _flux[8];
  }
  
  return _flux;
}

//////////////////////////////////////////////////////////////////////////////
  
void NavierStokes3DGammaAlphaPuvt::computeStressTensor(const RealVector& state,
							       const std::vector<RealVector*>& gradients,
							       const CFreal& radius)
{  
  using namespace std;
  using namespace COOLFluiD::Framework;
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  CFreal divTerm = 0.0;
  if (dim == DIM_2D && radius > MathTools::MathConsts::CFrealEps())  {
    // if the face is a boundary face, the radius could be 0
    // check against eps instead of 0. for safety
    divTerm = state[this->_vID]/radius;
  }
  else if (dim == DIM_3D) {
    const RealVector& gradW = *gradients[this->_wID];
    divTerm = gradW[ZZ];
  }
  
  CFreal twoThirdRhoK = 0.;
  if (getNbTurbVars() > 1) {
    // Spalart-Allmaras case
    const CFreal avK = state[_kID];
    
    const CFreal gamma = std::min(std::max(state[_kID+2],0.01),0.99);
    const CFreal gammaTerm = gamma + gamma*(1.0-gamma);
  
    twoThirdRhoK = this->_twoThird*this->getDensity(state)*avK*gammaTerm;
  }
  
  const RealVector& nsData = getModel().getPhysicalData();
  const CFreal coeffTauMu = getModel().getCoeffTau()*(nsData[NSTurbTerm::MU] + nsData[NSTurbTerm::MUT]);
  const RealVector& gradU = *gradients[this->_uID];
  const RealVector& gradV = *gradients[this->_vID]; 
 
  const CFreal twoThirdDivV = this->_twoThird*(gradU[XX] + gradV[YY] + divTerm);
  
  this->_tau(XX,XX) = coeffTauMu*(2.*gradU[XX] - twoThirdDivV)-twoThirdRhoK; 
  this->_tau(XX,YY) = this->_tau(YY,XX) = coeffTauMu*(gradU[YY] + gradV[XX]);
  this->_tau(YY,YY) = coeffTauMu*(2.*gradV[YY] - twoThirdDivV)-twoThirdRhoK;
  
  // KB & AL: This term should be wrong but KOmega in old COOLFLUID gives different results if this term is commented out!! 
  // this->_tau(XX,XX) += (1. - getModel().getCoeffTau()*nsData[NSTurbTerm::MUT]) * twoThirdRhoK;
  // this->_tau(YY,YY) += (1. - getModel().getCoeffTau()*nsData[NSTurbTerm::MUT]) * twoThirdRhoK;
  
  if (dim == DIM_3D) {
    const RealVector& gradW = *gradients[this->_wID];
    this->_tau(XX,ZZ) = this->_tau(ZZ,XX) = coeffTauMu*(gradU[ZZ] + gradW[XX]);
    this->_tau(YY,ZZ) = this->_tau(ZZ,YY) = coeffTauMu*(gradV[ZZ] + gradW[YY]);
    this->_tau(ZZ,ZZ) = coeffTauMu*(2.*gradW[ZZ] - twoThirdDivV)-twoThirdRhoK; 
    
    // KB & AL: This term should be wrong but KOmega in old COOLFLUID gives different results if this term is commented out!! 
    // this->_tau(ZZ,ZZ) += (1. - getModel().getCoeffTau()*nsData[NSTurbTerm::MUT]) * twoThirdRhoK;
  }
}

//////////////////////////////////////////////////////////////////////////////


    } // namespace GammaAlpha

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
