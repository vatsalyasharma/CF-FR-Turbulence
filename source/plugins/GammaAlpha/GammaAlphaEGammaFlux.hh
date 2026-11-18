#ifndef COOLFluiD_Physics_GammaAlpha_GammaAlphaEGammaFlux_hh
#define COOLFluiD_Physics_GammaAlpha_GammaAlphaEGammaFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <cmath>

#include "NavierStokes/TransitionArithmetic.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GammaAlpha {

//////////////////////////////////////////////////////////////////////////////

/// Mass-flux density rho*w of the intermittency dissipation of Jacobs, Van den Eynde & Steelant (2025, eq. 37,
/// boundary-layer branch) written as a transport of gamma away from the wall, as Jatayu v1.5 evaluates it:
///   E_gamma = r dgamma/dn,  r = C_Eg mu_gamma 2 c sqrt(max(0, 1 - c/c_e)) rho/mu,
///   mu_gamma = mu C_mug1 (g/(1-g)^7)^(1/6) (1 - g^5) f_M   (eq. 31, g = gamma kept in [0.01, 0.99]),
/// so that -E_gamma = -rho w.grad(gamma) = -div(rho w gamma) + gamma div(rho w) with rho w = r n. The first part
/// is transported by the convective operator (upwinded at the faces), the second is a local source. n is the flow
/// direction rotated by +90 degrees, (-v, u)/c: away from a wall on the right of the flow (the T3A orientation, as
/// the N_alpha flux form). The viscosity cancels, and the speed with the direction:
///   rho w = C_Eg C_mug1 shape(g) f_M 2 sqrt(max(0, 1 - c/c_e)) rho (-v, u).
/// Not included: the wall damping f_mugamma of eq. 29 (taken as 1; it departs from 1 only within
/// y < ~0.1 mu_inf/(rho_inf c_e) of the wall) and the stagnation branch of eq. 37 (zero with this n).
/// c_e and M_e: the edge state of Appendix B from the local pressure and the inlet totals.
inline void eGammaMassFlux(const CFreal cEg, const CFreal rho, const CFreal u, const CFreal v, const CFreal p,
                           const CFreal intermittency, const CFreal R, const CFreal kappa, const CFreal machInf,
                           const CFreal uInf, const CFreal pInf, CFreal& fx, CFreal& fy)
{
  const NavierStokes::transition_detail::GammaAlphaFreestream edge =
    NavierStokes::transition_detail::gammaAlphaFreestreamInlet(std::max(p, 1.0e-3), R, kappa, machInf, uInf, pInf);
  const CFreal c = std::sqrt(u*u + v*v);
  const CFreal layer = 2.0*std::sqrt(std::max(0.0, 1.0 - c/edge.speed));
  const CFreal g = std::min(std::max(intermittency, 0.01), 0.99);
  const CFreal shape = std::pow(g/std::pow(1.0 - g, 7.0), 1.0/6.0)*(1.0 - std::pow(g, 5.0));
  const CFreal fM = (1.0 + 0.26*0.5*(kappa - 1.0)*edge.mach*edge.mach)*std::sqrt(1.0 + 0.38*std::pow(edge.mach, 0.6));
  const CFreal coefficient = cEg*0.076*shape*fM*layer*rho;
  fx = -coefficient*v;
  fy = coefficient*u;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace GammaAlpha

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GammaAlpha_GammaAlphaEGammaFlux_hh
