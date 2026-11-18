#ifndef COOLFluiD_FluxReconstructionMethod_AUSMPlusUpFluxGammaAlpha2D_hh
#define COOLFluiD_FluxReconstructionMethod_AUSMPlusUpFluxGammaAlpha2D_hh

#include "Framework/CFSide.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler2DPuvt.hh"
#include "FluxReconstructionNavierStokes/AUSMPlusUpFlux.hh"

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

/**
 * AUSM+up flux for the 2D gamma-alpha model with the N_alpha term of the alpha equation in flux form.
 *
 * N_alpha = C_Nalpha rho |U| dalpha/dn (Van den Eynde & Steelant; Jacobs et al. 2025, eq. A16), with n the flow
 * direction rotated by +90 degrees, equals -rho w.grad(alpha) with w = C_Nalpha (v, -u): a wall-ward transport
 * of alpha. Evaluated pointwise as a source it is neither upwinded nor implicit and was unstable in FR (V&V runs
 * vv029-vv031). Here the conservative part rho*alpha*w is added to the alpha interface flux, upwinded on the
 * averaged w.n; Euler2DGammaAlphaPuvt adds it to the physical flux and GammaAlpha2DSourceTerm adds the remaining
 * alpha*div(rho w), so the alpha equation is unchanged. Active only when nAlphaFluxCoefficient() > 0.
 *
 * The same for E_gamma of the 2025 model (eq. 37) when eGammaFluxCoefficient() > 0: E_gamma = rho w_E.grad(gamma)
 * is a transport of gamma away from the wall (GammaAlphaEGammaFlux.hh); gamma*rho w_E is added to the gamma
 * interface flux, upwinded on the mean of the two sides' rho w_E.n, and GammaAlpha2DSourceTerm adds
 * gamma*div(rho w_E). As in Jatayu v1.5, where E_gamma is an upwind transport along the wall normal.
 */
class AUSMPlusUpFluxGammaAlpha2D :
  public AUSMPlusUpFlux<Physics::NavierStokes::MultiScalarVarSet<Physics::NavierStokes::Euler2DPuvt> > {
public:
  typedef AUSMPlusUpFlux<Physics::NavierStokes::MultiScalarVarSet<Physics::NavierStokes::Euler2DPuvt> > Base;

  AUSMPlusUpFluxGammaAlpha2D(const std::string& name);
  virtual ~AUSMPlusUpFluxGammaAlpha2D();

  using Base::computeFlux;
  virtual RealVector& computeFlux(Framework::State& lState, Framework::State& rState, const RealVector& normal);

  /// C_Egamma of the upwind non-conservative form of E_gamma (GammaAlpha2DSourceTerm, option EGammaUpwind): the
  /// source has -rho w_E.grad(gamma) with the central (BR2) gradient, and the gamma interface flux gets the upwind
  /// dissipation -|rho w_E.n|/2 (gamma_R - gamma_L), the difference between an upwind and a central flux. Both
  /// vanish for uniform gamma; at P0 the pair is first-order upwind. 0 = off (set by GammaAlpha2DSourceTerm).
  /// With option EGammaFaceDissipation the dissipation is added (times that factor) to the pointwise E_gamma.
  static CFreal& eGammaUpwindCoefficient();
};

  } // namespace FluxReconstructionMethod
} // namespace COOLFluiD

#endif
