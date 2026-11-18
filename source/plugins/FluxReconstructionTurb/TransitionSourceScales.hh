#ifndef COOLFluiD_FluxReconstructionMethod_TransitionSourceScales_hh
#define COOLFluiD_FluxReconstructionMethod_TransitionSourceScales_hh

#include <cmath>
#include <algorithm>
#include <limits>
#include <sstream>
#include "NavierStokes/TransitionArithmetic.hh"
#include "Common/BadValueException.hh"
#include "Common/COOLFluiD.hh"

namespace COOLFluiD {
namespace FluxReconstructionMethod {
namespace detail {

// Input admissibility, not a bound on every downstream correlation/product.
// Call for every cell state before writing residuals, and before clearing a
// point Jacobian. Reject rather than clip: min/max can silently hide NaNs.
// Omega must be normal and finite; its reciprocal is then representable.
// The upper log endpoint is excluded because log(max) can round upward.
template<class Vector>
inline void validateTransitionSourceState(const Vector& state, const CFuint dim,
                                          const CFuint stateID,
                                          const bool requireNonnegativeK = true)
{
  const CFuint temperature = dim + 1;
  const CFuint k = dim + 2;
  const CFuint logOmega = dim + 3;
  const char* reason = 0;
  CFuint component = 0;
  for (CFuint i = 0; i < state.size(); ++i) {
    if (!std::isfinite(state[i])) {
      reason = "nonfinite state component"; component = i; break;
    }
  }
  if (!reason && !(state[0] > 0.0))
    { reason = "pressure must be positive"; component = 0; }
  if (!reason && !(state[temperature] > 0.0))
    { reason = "temperature must be positive"; component = temperature; }
  // A signed-k model extension evaluates all k-dependent physics with max(k,0);
  // high-order polynomials may then undershoot zero (cf. negative SA).
  if (!reason && requireNonnegativeK && state[k] < 0.0)
    { reason = "turbulent kinetic energy must be nonnegative"; component = k; }
  const CFreal lower = std::log(std::numeric_limits<CFreal>::min());
  const CFreal upper = std::log(std::numeric_limits<CFreal>::max());
  if (!reason && !(state[logOmega] >= lower && state[logOmega] < upper))
    { reason = "logomega outside normal finite exponential range"; component = logOmega; }
  if (reason) {
    std::ostringstream message;
    message << "Transition source inadmissible state: local ID " << stateID
            << ", component " << component << ", value " << state[component]
            << ": " << reason;
    throw Common::BadValueException(FromHere(), message.str());
  }
}

template<class Gradients>
inline void validateTransitionSourceEnvironment(const Gradients& gradients,
                                                const CFreal wallDistance,
                                                const CFuint dim, const CFuint equations)
{
  if (!(std::isfinite(wallDistance) && wallDistance >= 0.))
    throw Common::BadValueException(FromHere(), "Transition source requires finite nonnegative wall distance");
  if (gradients.size() != equations)
    throw Common::BadValueException(FromHere(), "Transition source gradient equation count mismatch");
  for (CFuint e = 0; e < equations; ++e) {
    if (gradients[e].size() != dim)
      throw Common::BadValueException(FromHere(), "Transition source gradient dimension mismatch");
    for (CFuint d = 0; d < dim; ++d)
      if (!std::isfinite(gradients[e][d]))
        throw Common::BadValueException(FromHere(), "Transition source requires finite gradients");
  }
}

// Residual-domain check, before any source evaluation. Pressure-gradient
// correlations require a streamline direction; none is selected at rest.
inline CFreal transitionResidualSpeed(const CFreal u, const CFreal v,
                                      const CFreal w, const bool pressureGradient)
{
  const CFreal speed = std::hypot(std::hypot(u, v), w);
  if (!std::isfinite(speed))
    throw Common::BadValueException(FromHere(), "Transition source requires finite velocity");
  if (speed == 0.0 && pressureGradient)
    throw Common::BadValueException(FromHere(),
        "Transition pressure-gradient source requires positive speed (streamline direction undefined at rest)");
  return speed;
}

// No zero-speed analytical derivative policy is implemented. In GA the
// speed-norm source is generally nondifferentiable there; in GR this is a
// conservative support restriction, not a claim of universal nondifferentiability.
inline CFreal transitionAnalyticalSpeed(const CFreal u, const CFreal v, const CFreal w)
{
  const CFreal speed = transitionResidualSpeed(u, v, w, false);
  if (speed == 0.0)
    throw Common::BadValueException(FromHere(),
        "Transition analytical source Jacobian requires positive speed; zero-speed derivative policy is not implemented");
  return speed;
}

// Saturate before dividing, so a bounded intensity never requires an
// overflowing intermediate. k=V=0 has no unique intensity limit.
inline CFreal transitionIntensity(const CFreal k, const CFreal speed, const CFreal cap)
{
  const CFreal fluctuation = std::sqrt((2.0/3.0)*k);
  if (speed == 0.0 && fluctuation == 0.0)
    throw Common::BadValueException(FromHere(), "Transition intensity undefined at zero k and zero speed");
  if (speed <= (100.0/cap)*fluctuation) return cap;
  return std::max(0.027, 100.0*(fluctuation/speed));
}

// Use 1/t directly: t = 500 mu/(rho V^2) is infinite at rest.
// Preconditions: finite, positive density/viscosity and finite speed >= 0.
inline CFreal transitionRelaxationRate(const CFreal rho, const CFreal mu,
                                      const CFreal speed)
{
  return (rho/(500.0*mu)*speed)*speed;
}

// In both models delta = 375 Omega y Re_theta mu/(rho V^2).
// Cancel y before evaluating exp(-(y/delta)^4); neither theta_BL nor delta
// needs to be formed. This also supplies the y -> 0 continuation.
// Preconditions: finite positive rho, mu, Re_theta; finite Omega, V >= 0.
inline CFreal transitionWakeExponential(const CFreal rho, const CFreal mu,
                                       const CFreal retheta, const CFreal vorticity,
                                       const CFreal speed)
{
  if (vorticity == 0.0) {
    if (speed == 0.0)
      throw Common::BadValueException(FromHere(),
          "Transition wake factor has no unique joint limit at zero speed and zero vorticity");
    return 0.0;
  }
  if (speed == 0.0) return 1.0;
  const CFreal ratio = (rho/(375.0*retheta*mu)*speed)*(speed/vorticity);
  // exp(-8^4) rounds to zero in CFreal; avoid overflowing the fourth power.
  if (ratio >= 8.0) return 0.0;
  const CFreal ratio2 = ratio*ratio;
  return std::exp(-ratio2*ratio2);
}

// V*d(alpha)/dn = v*d(alpha)/dx - u*d(alpha)/dy in the 2D GA source.
// The product has a unique zero-speed limit although the normal does not.
inline CFreal gammaAlphaCrossGradient(const CFreal u, const CFreal v,
                                     const CFreal dax, const CFreal day)
{
  return v*dax-u*day;
}

} // namespace detail
} // namespace FluxReconstructionMethod
} // namespace COOLFluiD
#endif
