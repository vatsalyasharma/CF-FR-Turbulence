#include "NavierStokes/TransitionTransport.hh"
#include "NavierStokes/TransitionArithmetic.hh"
#include "NavierStokes2DGammaAlphaSSTPuvt.hh"
#include "Common/BadValueException.hh"
#include <typeinfo>
#include "GammaAlpha.hh"
#include "NavierStokes2DGammaAlphaPuvt.hh"
#include "GammaAlpha/GammaAlphaMuTFactor.hh"
#include "GammaAlpha/GammaAlphaArtificialDiffusion.hh"
#include "GammaAlpha/Euler2DGammaAlphaPuvt.hh"
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

    namespace GammaAlpha {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NavierStokes2DGammaAlphaPuvt, DiffusiveVarSet,
			    GammaAlphaModule, 2>
ns2DGammaAlphaPuvtProvider("NavierStokes2DGammaAlphaPuvt");

//////////////////////////////////////////////////////////////////////////////

NavierStokes2DGammaAlphaPuvt::NavierStokes2DGammaAlphaPuvt
(const std::string& name, SafePtr<PhysicalModelImpl> model) :
  NavierStokesKLogOmegaPvt<NavierStokesKLogOmegaVarSet<NavierStokes2DVarSet, 0> >(name, model),
  _unperturbedFluxGa(),
  _unperturbedFluxAlpha()
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
  names[6] = "gamma";
  names[7] = "alpha";
  
  setVarNames(names);
  setModelCoefficients();
}

CFreal NavierStokes2DGammaAlphaPuvt::getLaminarDynViscosityFromGradientVars(const RealVector& state)
{
  return transition_detail::laminarViscosity(getModel(), state[0], state[3],
      _eulerModel->getPressRef(), _eulerModel->getTempRef());
}

CFreal NavierStokes2DGammaAlphaPuvt::getDensity(const RealVector& state)
{
  transition_detail::validateThermodynamics(state[0], state[3], _eulerModel->getR(), _eulerModel->getGamma());
  return state[0]/(_eulerModel->getR()*state[3]);
}

CFreal NavierStokes2DGammaAlphaPuvt::getRhoKOverMuT(const RealVector& state, const std::vector<RealVector*>& gradients)
{
  // Dispatch only the implemented constitutive laws; never infer a zero-mut
  // limit from a floor or a synthetic positive-k state.
  if (typeid(*this) == typeid(NavierStokes2DGammaAlphaSSTPuvt))
    return static_cast<NavierStokes2DGammaAlphaSSTPuvt*>(this)->getRhoKOverMuT(state, gradients);
  if (typeid(*this) == typeid(NavierStokes2DGammaAlphaBSLPuvt))
    return static_cast<NavierStokes2DGammaAlphaBSLPuvt*>(this)->getRhoKOverMuT(state, gradients);
  if (typeid(*this) == typeid(NavierStokes2DGammaAlphaPuvt))
    return std::exp(state[5]);
  throw Common::BadValueException(FromHere(), "Transition source viscosity ratio is not implemented for this diffusive variable set");
}

//////////////////////////////////////////////////////////////////////////////

NavierStokes2DGammaAlphaPuvt::~NavierStokes2DGammaAlphaPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

RealVector& NavierStokes2DGammaAlphaPuvt::getFlux(const RealVector& values,
                                              const vector<RealVector*>& gradients,
                                              const RealVector& normal,
                                              const CFreal& radius)
{
  // Check before transport/stress evaluation or writing the native flux buffer.
  const bool v2025 = gammaAlphaModelVersion() == 2025;
  const auto freestream = v2025 ?
    transition_detail::gammaAlphaFreestreamInlet(std::max(values[0],1.0e-5),
      _eulerModel->getR(), _eulerModel->getGamma(), _eulerModel->getMachInf(),
      _eulerModel->getVelInf(), _eulerModel->getStaticPressInf()) :
    transition_detail::gammaAlphaFreestream(std::max(values[0],1.0e-5),
      _eulerModel->getR(), _eulerModel->getGamma(), _eulerModel->getMachInf(),
      _eulerModel->getVelInf(), _eulerModel->getStaticPressInf());
  // here nsData is filled in
  getLaminarDynViscosityFromGradientVars(values);
  transition_detail::laminarViscosity(getModel(), values[0], freestream.temperature,
      _eulerModel->getPressRef(), _eulerModel->getTempRef());
  computeTransportProperties(values, gradients, normal);
  computeStressTensor(values, gradients, radius);
  
  const RealVector& gradT     = *gradients[3];
  const RealVector& gradK     = *gradients[4];
  const RealVector& gradOmega = *gradients[5];
  const RealVector& gradGa    = *gradients[6];
  const RealVector& gradAlpha    = *gradients[7];

  const CFreal p = max(values[0],1.0e-5);
  const CFreal avU = values[1];
  const CFreal avV = values[2];
  const CFreal gamma = min(max(values[6],0.01),0.99);
  //const CFreal gamma = min(max(values[6],0.0),1.0);

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
  // KOmegaDiffusionSST (GammaAlpha2DSourceTerm): k and omega diffuse with the SST mu_t = mu_t*/f(gamma) (A4)
  const CFreal kOmegaMutFactor = gammaAlphaKOmegaDiffusionSST() ? 1.0/gammaAlphaMeanFlowMuTFactor(gamma) : 1.0;
  CFreal coeffTauMu1 = getModel().getCoeffTau()*nsData[NSTurbTerm::MU];
  CFreal coeffTauMu2 = getModel().getCoeffTau()*nsData[NSTurbTerm::MUT]*kOmegaMutFactor*getSigmaK();

  _flux[4] = (coeffTauMu1 + coeffTauMu2)*(gradK[XX]*nx + gradK[YY]*ny);
  //diff flux corresponding to Omega
  coeffTauMu1 = getModel().getCoeffTau()*nsData[NSTurbTerm::MU];
  coeffTauMu2 = getModel().getCoeffTau()*nsData[NSTurbTerm::MUT]*kOmegaMutFactor*getSigmaOmega();

  _flux[5] = (coeffTauMu1 + coeffTauMu2)*(gradOmega[XX]*nx + gradOmega[YY]*ny);
  
  //diff flux corresponding to Gamma
  const CFreal sigmaf= 1.0;
  coeffTauMu1 = getModel().getCoeffTau()*nsData[NSTurbTerm::MU];
  coeffTauMu2 = getModel().getCoeffTau()*nsData[NSTurbTerm::MUT]/sigmaf;
  
  const CFreal gammaIsentropic = _eulerModel->getGamma();
  const CFreal rhoInfLocal = freestream.density;
  const CFreal MInfLocal = freestream.mach;
  const CFreal uInfLocal = freestream.speed;
  const CFreal TInfLocal = freestream.temperature;
  
  // RealVector copy construction can retain a non-owning state view.
  // Allocate owned storage before changing the freestream temperature.
  RealVector copyState(values.size());
  copyState = values;
  copyState[3] = TInfLocal;
  
  const CFreal muInfLocal = getLaminarDynViscosityFromGradientVars(copyState);
  
  //const CFreal distance = std::max(_wallDistance, 1.0e-12);
  
  const CFreal fMuGamma = 1.0-exp(-256.0*(_wallDistance*uInfLocal*rhoInfLocal/muInfLocal)*(_wallDistance*uInfLocal*rhoInfLocal/muInfLocal));
  const CFreal fMMuGamma = (1.0+0.26*(gammaIsentropic-1.0)/2.0*MInfLocal*MInfLocal)*sqrt(1+0.38*pow(MInfLocal,0.6));
 //if (_wallDistance < 1.0e-4) CFLog(INFO, "wallDist: " << _wallDistance << "\n"); 
  //const CFreal gammaLim = min(max(0.01,gamma),0.99);
  // 2018: thesis eq. 7.65; 2025: Jacobs et al. eq. 31 (exact Klebanoff profile, C1 = 0.076, n = 5)
  const CFreal muGammaShape = v2025 ? 0.076*pow(gamma/pow(1.0-gamma,7.0),1.0/6.0)*(1.0-pow(gamma,5.0))
                                    : 0.57*pow(-log(1-gamma),-5.0/6.0*(1.0-gamma));
  const CFreal muGamma = muGammaShape*fMuGamma*fMMuGamma*getModel().getCoeffTau()*nsData[NSTurbTerm::MU];

  _flux[6] = (muGamma)*(gradGa[XX]*nx + gradGa[YY]*ny);
  // EGammaArtificialDiffusion (GammaAlpha2DSourceTerm): artificial diffusivity of the gamma equation at this point
  const CFreal gammaArtDiff = gammaArtificialDiffusion().current;
  if (gammaArtDiff != 0.0) _flux[6] += getModel().getCoeffTau()*gammaArtDiff*(gradGa[XX]*nx + gradGa[YY]*ny);
  
  //diff flux corresponding to alpha
  const CFreal sigmatheta= 2.0;
  coeffTauMu1 = getModel().getCoeffTau()*nsData[NSTurbTerm::MU];
  coeffTauMu2 = getModel().getCoeffTau()*nsData[NSTurbTerm::MUT];
  if (v2025) {
    // eq. A15 uses the turbulence-model mu_t; nsData holds the effective mu_t* (eq. A4) of the SST variable set
    const CFreal effective = gammaAlphaMeanFlowMuTFactor(gamma);
    coeffTauMu2 /= effective;
  }

  _flux[7] = sigmatheta*(coeffTauMu1 + coeffTauMu2)*(gradAlpha[XX]*nx + gradAlpha[YY]*ny);
  
  //if (values[6]<0) CFLog(INFO, "Flux update: " << _flux << "\n");
  
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
      _flux[7] = _unperturbedFluxAlpha;
    }
  }
  else
  {
    _unperturbedFluxK = _flux[4];
    _unperturbedFluxOmega = _flux[5];
    _unperturbedFluxGa    = _flux[6];
    _unperturbedFluxAlpha    = _flux[7];
  }
 //cout << _flux << endl;
  return _flux;
}

//////////////////////////////////////////////////////////////////////////////
  
void NavierStokes2DGammaAlphaPuvt::computeStressTensor(const RealVector& state,
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
