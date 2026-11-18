#include "NavierStokes/TransitionArithmetic.hh"
#include "Common/BadValueException.hh"
#include "GammaAlpha.hh"
#include "NavierStokes2DGammaAlphaSSTPuvt.hh"
#include "GammaAlpha/GammaAlphaMuTFactor.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GammaAlpha {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NavierStokes2DGammaAlphaSSTPuvt, DiffusiveVarSet,
	       GammaAlphaModule, 2>
ns2DGammaAlphaSSTPuvtProvider("NavierStokes2DGammaAlphaSSTPuvt");

//////////////////////////////////////////////////////////////////////////////

NavierStokes2DGammaAlphaSSTPuvt::NavierStokes2DGammaAlphaSSTPuvt
(const std::string& name, SafePtr<PhysicalModelImpl> model) :
  NavierStokes2DGammaAlphaBSLPuvt(name, model)
{
}

//////////////////////////////////////////////////////////////////////////////

NavierStokes2DGammaAlphaSSTPuvt::~NavierStokes2DGammaAlphaSSTPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DGammaAlphaSSTPuvt::setModelCoefficients()
{

  NavierStokes2DGammaAlphaPuvt::setModelCoefficients();

  //Modified k-omega coeficients for SST
  _sigmaK1 = 0.85;
 // _sigmaK2 = 1.0;
 // _sigmaOmega2 = 0.856;
 // _beta2 = 0.0828;
  
 // Menter SST model from 2003  
 // _gamma1= 5/9;
 // _gamma2= 0.44;
}

//////////////////////////////////////////////////////////////////////////////

CFreal NavierStokes2DGammaAlphaSSTPuvt::getTurbDynViscosityFromGradientVars(const RealVector& state, const vector<RealVector*>& gradients)
{
  return computeViscosityOrRatio(state, gradients, false);
}

CFreal NavierStokes2DGammaAlphaSSTPuvt::getRhoKOverMuT(const RealVector& state, const vector<RealVector*>& gradients)
{
  return computeViscosityOrRatio(state, gradients, true);
}

CFreal NavierStokes2DGammaAlphaSSTPuvt::computeViscosityOrRatio(const RealVector& state, const vector<RealVector*>& gradients, const bool ratio)
{

  if (ratio && !(_wallDistance > 0. && gradients.size() > 0))
    throw Common::BadValueException(FromHere(), "Transition viscosity ratio requires positive wall distance and supplied gradients");
  cf_assert(_wallDistance >= 0.);
  //cf_assert(gradients.size() > 0);

  CFreal mut = 0.;
  if((_wallDistance > 0.) && (gradients.size() > 0))
  {
    const CFreal a1 = 0.31;

    const CFreal rho = std::max(0., getDensity(state));
    const CFreal K = std::max(state[4],0.0);
    const CFreal Omega = std::exp(state[5]);
    
    //const CFreal vorticity = fabs((*(gradients[2]))[XX] - (*(gradients[1]))[YY]);
    CFreal velocityGradients[3][3] = {};
    for (CFuint e=0; e<2; ++e) for (CFuint d=0; d<2; ++d)
      velocityGradients[e][d] = (*gradients[e+1])[d];
    const CFreal strain = transition_detail::velocityGradientNorm(velocityGradients,2,false);

    const CFreal mu = getLaminarDynViscosityFromGradientVars(state);
    const CFreal avGa = std::min(std::max(state[6],0.01),0.99);
    //const CFreal avGa = std::min(std::max(state[6],0.01),1.0);

    ///@todo here there should be an adimensionalization coef
    // Cap each nonnegative ratio before arithmetic; min/max commute with
    // this cap, and the downstream tanh is already saturated here.
    const CFreal arg2_1 = transition_detail::cappedPositiveProductRatio({2., sqrt(K)}, {0.09, Omega, _wallDistance}, 8.);
    const CFreal arg2_2 = transition_detail::cappedPositiveProductRatio({500., mu}, {rho, _wallDistance, _wallDistance, Omega}, 8.);
    const CFreal arg2 = std::max(arg2_1,arg2_2);

    const CFreal F2 = transition_detail::tanhSquare(arg2);


    
    // compute effective eddy viscosity
    const CFreal cmut1 = 0.5;//0.7;//0.68;//
    const CFreal cmut2 = 0.25;
    const CFreal cmut3 = 0.1;
    // A4 (cmut1 = 0.5, cmut2 = 0.25, cmut3 = 0.1) times the laminar flow share (GammaAlphaMuTFactor.hh; 1 = A4)
    if (ratio) return transition_detail::finiteRatio(std::max(a1*Omega, strain*F2),(a1*gammaAlphaMeanFlowMuTFactor(avGa)));
    mut = (a1 * rho * K )/ std::max(a1*Omega , strain*F2);
    mut *= gammaAlphaMeanFlowMuTFactor(avGa);
  }

  return mut;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace GammaAlpha

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

