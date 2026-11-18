// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
// Distributed under the GNU Lesser General Public License version 3 (LGPLv3).
#ifndef COOLFluiD_NewtonMethod_TurbulenceUpdatePolicy_hh
#define COOLFluiD_NewtonMethod_TurbulenceUpdatePolicy_hh

#include <cmath>
#include <string>
#include "Common/BadValueException.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/SpaceMethodData.hh"
#include "NewtonMethod/NewtonIteratorData.hh"

namespace COOLFluiD { namespace Numerics { namespace NewtonMethod {
namespace detail {

// Both raw and logarithmic variable sets label the component "Omega".
// Identify only known model/layout combinations; never infer from that label.
inline bool turbulenceStoresLogOmega(NewtonIteratorData& data,
                                    const bool transitionUpdater,
                                    const CFreal omegaInlet,
                                    const std::string& command)
{
  const CFuint dim = Framework::PhysicalModelStack::getActive()->getDim();
  const CFuint eqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  const std::string model = Framework::PhysicalModelStack::getActive()->getConvectiveName();
  const std::string vars = data.getCollaborator<Framework::SpaceMethod>()
    ->getSpaceMethodData()->getUpdateVarStr();
  const std::string prefix = dim == 2 ? "Euler2D" : "Euler3D";
  const bool primitive = (dim == 2 && vars == "Puvt") || (dim == 3 && vars == "Pvt");
  const bool rawKO = model == prefix + "KOmega";
  const bool logKO = model == prefix + "KLogOmega";
  const bool rawGR = model == prefix + "GReKO";
  const bool logGR = model == prefix + "GReKLogO";
  const bool ga = model == prefix + "GammaAlpha";
  const bool supported = transitionUpdater ? (rawGR || logGR) : (rawKO || logKO || rawGR || logGR || ga);
  const CFuint scalars = (rawKO || logKO) ? 2 : 4;
  if (!primitive || !supported || eqs != dim + 2 + scalars)
  {
    throw Common::BadValueException(FromHere(), command +
      ": unsupported turbulence update layout " + model + "/" + vars +
      ". Requires a known primitive p,velocity,T,k,omega/log(omega) layout; "
      "GReKOUpdateSol additionally requires gamma,ReTheta (not gamma,alpha).");
  }
  if (!std::isfinite(omegaInlet) || omegaInlet <= 0.)
  {
    throw Common::BadValueException(FromHere(), command +
      ": OmegaInlet must be finite and positive, specified as unlogged omega in solver units.");
  }
  return logKO || logGR || ga;
}

} // namespace detail
} } } // namespace COOLFluiD::Numerics::NewtonMethod
#endif
