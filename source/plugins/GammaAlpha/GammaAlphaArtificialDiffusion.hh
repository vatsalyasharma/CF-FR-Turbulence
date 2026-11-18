#ifndef COOLFluiD_Physics_GammaAlpha_GammaAlphaArtificialDiffusion_hh
#define COOLFluiD_Physics_GammaAlpha_GammaAlphaArtificialDiffusion_hh

//////////////////////////////////////////////////////////////////////////////

#include <vector>

#include "Common/COOLFluiD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GammaAlpha {

//////////////////////////////////////////////////////////////////////////////

/// Artificial diffusivity [kg/(m s)] of the gamma equation (GammaAlpha2DSourceTerm, option EGammaArtificialDiffusion):
/// the source term fills one value per state at the start of each residual evaluation; the turbulence diffusion
/// commands (FluxReconstructionTurb) select the value of the point where the diffusive flux is evaluated next, as they
/// do for the wall distance (solution point: its state; interior flux point: the mean of the two cells' solution
/// points closest to it; boundary flux point: zero), and NavierStokes2DGammaAlphaPuvt::getFlux adds current*grad(gamma).
/// Header-only: the inline accessor's static is shared by all libraries that include it (GNU unique symbol).
struct GammaArtificialDiffusion
{
  std::vector<CFreal> perState;
  CFreal current;
  bool active;
  GammaArtificialDiffusion() : perState(), current(0.), active(false) {}
};

inline GammaArtificialDiffusion& gammaArtificialDiffusion()
{
  static GammaArtificialDiffusion data;
  return data;
}

inline CFreal gammaArtificialDiffusionAt(const CFuint stateID)
{
  const GammaArtificialDiffusion& d = gammaArtificialDiffusion();
  return (d.active && stateID < d.perState.size()) ? d.perState[stateID] : 0.;
}

/// at a solution point
inline void setGammaArtificialDiffusionAtState(const CFuint stateID)
{
  gammaArtificialDiffusion().current = gammaArtificialDiffusionAt(stateID);
}

/// at an interior flux point, from the solution points of the two cells closest to it
inline void setGammaArtificialDiffusionAtFace(const CFuint stateIDL, const CFuint stateIDR)
{
  gammaArtificialDiffusion().current = 0.5*(gammaArtificialDiffusionAt(stateIDL) + gammaArtificialDiffusionAt(stateIDR));
}

/// at a boundary flux point: no artificial flux through boundaries
inline void setGammaArtificialDiffusionZero()
{
  gammaArtificialDiffusion().current = 0.;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace GammaAlpha

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GammaAlpha_GammaAlphaArtificialDiffusion_hh
