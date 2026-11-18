#ifndef COOLFluiD_Physics_GammaAlpha_GammaAlphaMuTFactor_hh
#define COOLFluiD_Physics_GammaAlpha_GammaAlphaMuTFactor_hh

//////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <cmath>

#include "GammaAlpha/Euler2DGammaAlphaPuvt.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GammaAlpha {

//////////////////////////////////////////////////////////////////////////////

/// mu_t*/mu_t of the mean flow: A4 of Jacobs, Van den Eynde & Steelant (2025) times the laminar flow share
/// (Jatayu v1.5 gamma_alpha_laminar_flow_share, gammaAlphaLaminarFlowShare(); 1 = A4 as printed): the share at
/// gamma = 0.01, all of mu_t* from gamma = 0.05 on, linear between. g is gamma kept in [0.01, 0.99]. The A4 part
/// is the expression of NavierStokes2DGammaAlphaSSTPuvt (bit for bit), so a share of 1 changes nothing.
inline CFreal gammaAlphaMeanFlowMuTFactor(const CFreal g)
{
  const CFreal a4 = g + g*(1.0-g)*0.5*(1.0+std::tanh((g-0.25)/0.1));
  const CFreal share = gammaAlphaLaminarFlowShare();
  if (!(share < 1.0)) return a4;
  return a4*std::min(std::max(share + (g - 0.01)/(0.05 - 0.01)*(1.0 - share), share), 1.0);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace GammaAlpha

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GammaAlpha_GammaAlphaMuTFactor_hh
