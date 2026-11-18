#include "NavierStokes/TransitionArithmetic.hh"
#include "Common/BadValueException.hh"
#include "GReKO.hh"
#include "NavierStokes3DGReKLogOSSTPuvt.hh"
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

    namespace GReKO {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NavierStokes3DGReKLogOSSTPuvt, DiffusiveVarSet,
	       GReKOModule, 2>
ns3DGReKLogOSSTPuvtProvider("NavierStokes3DGReKLogOSSTPvt");

//////////////////////////////////////////////////////////////////////////////

NavierStokes3DGReKLogOSSTPuvt::NavierStokes3DGReKLogOSSTPuvt
(const std::string& name, SafePtr<PhysicalModelImpl> model) :
  NavierStokes3DGReKLogOBSLPuvt(name, model)
{
}

//////////////////////////////////////////////////////////////////////////////

NavierStokes3DGReKLogOSSTPuvt::~NavierStokes3DGReKLogOSSTPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DGReKLogOSSTPuvt::setModelCoefficients()
{

  NavierStokes3DGReKLogOPuvt::setModelCoefficients();

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

CFreal NavierStokes3DGReKLogOSSTPuvt::getTurbDynViscosityFromGradientVars(const RealVector& state, const vector<RealVector*>& gradients)
{
  return computeViscosityOrRatio(state, gradients, false);
}

CFreal NavierStokes3DGReKLogOSSTPuvt::getRhoKOverMuT(const RealVector& state, const vector<RealVector*>& gradients)
{
  return computeViscosityOrRatio(state, gradients, true);
}

CFreal NavierStokes3DGReKLogOSSTPuvt::computeViscosityOrRatio(const RealVector& state, const vector<RealVector*>& gradients, const bool ratio)
{

  if (ratio && !(_wallDistance > 0. && gradients.size() > 0))
    throw Common::BadValueException(FromHere(), "Transition viscosity ratio requires positive wall distance and supplied gradients");
  cf_assert(_wallDistance >= 0.);
///@todo put this here??
  cf_assert(gradients.size() > 0);
//if (_wallDistance == 0.) CFLog(INFO, "wallDist: " << _wallDistance << "\n");
  CFreal mut = 0.;
  if((_wallDistance > 0.) && (gradients.size() > 0))
  {
    const CFreal a1 = 0.31;

    const CFreal rho = getDensity(state);
    const CFreal K = std::max(state[5],0.0);
    const CFreal Omega = std::exp(state[6]);
    
    //const CFreal vorticity = fabs((*(gradients[2]))[XX] - (*(gradients[1]))[YY]);
    CFreal velocityGradients[3][3] = {};
    for (CFuint e=0; e<3; ++e) for (CFuint d=0; d<3; ++d)
      velocityGradients[e][d] = (*gradients[e+1])[d];
    const CFreal strain = transition_detail::velocityGradientNorm(velocityGradients,3,false);
    
    const CFreal mu = getLaminarDynViscosityFromGradientVars(state);

    ///@todo here there should be an adimensionalization coef
    // Cap each nonnegative ratio before arithmetic; min/max commute with
    // this cap, and the downstream tanh is already saturated here.
    const CFreal arg2_1 = transition_detail::cappedPositiveProductRatio({2., sqrt(K)}, {0.09, Omega, _wallDistance}, 8.);
    const CFreal arg2_2 = transition_detail::cappedPositiveProductRatio({500., mu}, {rho, _wallDistance, _wallDistance, Omega}, 8.);
    const CFreal arg2 = std::max(arg2_1,arg2_2);

    const CFreal F2 = transition_detail::tanhSquare(arg2);

    if (ratio) return transition_detail::finiteRatio(std::max(a1*Omega, strain*F2),a1);
    mut = (a1 * rho * K )/ std::max(a1*Omega , strain*F2);
  }

  return mut;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace GReKO

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

