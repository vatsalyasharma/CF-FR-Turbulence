#ifndef COOLFluiD_FluxReconstructionMethod_GammaAlphaOnsetTable_hh
#define COOLFluiD_FluxReconstructionMethod_GammaAlphaOnsetTable_hh

#include <cmath>
#include <limits>
#include <map>
#include <utility>

#include "Common/COOLFluiD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Gamma-alpha onset margins at the wall flux points of this rank, (alpha - tau_w*)/alpha keyed by the flux point
 * coordinates (x, y). GammaAlphaTransitionCriterion fills it in the unperturbed boundary pass of every residual
 * evaluation; GammaAlpha2DSourceTerm (OnsetFunction = true) gathers it over the ranks and turns it into
 * F_onset (Jacobs, Van den Eynde & Steelant 2025, eq. A7: F_onset = 1 where tau_w* <= alpha_cr).
 */
inline std::map< std::pair<CFreal,CFreal>, CFreal >& gammaAlphaOnsetTable()
{
  static std::map< std::pair<CFreal,CFreal>, CFreal > table;
  return table;
}

/**
 * alpha and the pressure at the same wall flux points (GammaAlphaTransitionCriterion), keyed as
 * gammaAlphaOnsetTable(): the inputs of beta0 of the wall law (GammaAlpha2DSourceTerm, option WallLaw).
 */
inline std::map< std::pair<CFreal,CFreal>, std::pair<CFreal,CFreal> >& gammaAlphaWallPointTable()
{
  static std::map< std::pair<CFreal,CFreal>, std::pair<CFreal,CFreal> > table;
  return table;
}

/**
 * Wall value of gamma by the wall law at every wall flux point of all ranks, keyed as gammaAlphaOnsetTable().
 * Written by GammaAlpha2DSourceTerm (option WallLaw: the march along the wall in every residual evaluation), read by
 * BCNoSlipWallGammaAlpha2D (option WallLaw) in the next evaluation.
 */
inline std::map< std::pair<CFreal,CFreal>, CFreal >& gammaAlphaWallLawTable()
{
  static std::map< std::pair<CFreal,CFreal>, CFreal > table;
  return table;
}

/**
 * Entry of a table above at the point (x, y), matched within a relative tolerance: the flux point coordinates of the
 * convective and the diffusive boundary passes may differ in the last bits (seen after a P0 -> P1 upgrade).
 * Returns end() when no entry is that close.
 */
template <class Value>
inline typename std::map< std::pair<CFreal,CFreal>, Value >::const_iterator
findWallPoint(const std::map< std::pair<CFreal,CFreal>, Value >& table, const CFreal x, const CFreal y)
{
  const CFreal tolX = 1.0e-8*(1.0 + std::fabs(x));
  const CFreal tolY = 1.0e-8*(1.0 + std::fabs(y));
  typename std::map< std::pair<CFreal,CFreal>, Value >::const_iterator it =
    table.lower_bound(std::make_pair(x - tolX, -std::numeric_limits<CFreal>::max()));
  for (; it != table.end() && it->first.first <= x + tolX; ++it)
  {
    if (std::fabs(it->first.second - y) <= tolY) return it;
  }
  return table.end();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_GammaAlphaOnsetTable_hh
