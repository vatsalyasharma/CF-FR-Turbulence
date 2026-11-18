#ifndef COOLFluiD_Physics_NavierStokes_TransitionTransport_hh
#define COOLFluiD_Physics_NavierStokes_TransitionTransport_hh

#include "NavierStokes/TransitionArithmetic.hh"
#include "NavierStokes/NSTerm.hh"
#include "NavierStokes/SutherlandDynViscosity.hh"
#include "NavierStokes/SimplerSutherlandDynViscosity.hh"
#include "NavierStokes/FixedDynViscosity.hh"
#include "NavierStokes/ConstantDynViscosity.hh"
#include "NavierStokes/FixedKinViscosity.hh"
#include <typeinfo>

namespace COOLFluiD { namespace Physics { namespace NavierStokes {
namespace transition_detail {

// Normal dimensional and nondimensional viscosity are the supported transport
// domain. These tests select an evaluation path; they do not clip viscosity.
inline bool transportProduct(const CFreal a, const CFreal b, CFreal& out)
{
  if ((b > 1. && a >= std::numeric_limits<CFreal>::max()/b) ||
      (b < 1. && a <= std::numeric_limits<CFreal>::min()/b)) return false;
  out=a*b; return std::isnormal(out) && out>0.;
}
inline bool transportRatio(const CFreal a, const CFreal b, CFreal& out)
{
  if ((b < 1. && a >= std::numeric_limits<CFreal>::max()*b) ||
      (b > 1. && a <= std::numeric_limits<CFreal>::min()*b)) return false;
  out=a/b; return std::isnormal(out) && out>0.;
}
inline void requireTransportShift(const CFreal shift)
{
  if (!(std::isfinite(shift) && shift>=0.))
    throw Common::BadValueException(FromHere(), "Transition transport requires a finite nonnegative temperature shift");
}
inline CFreal transportLogSum(const CFreal a, const CFreal b)
{
  const CFreal hi=std::max(a,b), lo=std::min(a,b);
  return std::log(hi)+std::log1p(lo/hi);
}
inline CFreal transportExp(const CFreal logarithm)
{
  // Keep a one-log-unit margin at libm's rounded normal-range endpoints.
  if (!(std::isfinite(logarithm) && logarithm>std::log(std::numeric_limits<CFreal>::min())+1. &&
        logarithm<std::log(std::numeric_limits<CFreal>::max())-1.))
    throw Common::BadValueException(FromHere(), "Transition dimensional viscosity outside supported normal range");
  return std::exp(logarithm);
}
inline bool transportPower(const CFreal T, CFreal& out)
{
  const CFreal exponent=1.5*std::log(T);
  if (exponent<=std::log(std::numeric_limits<CFreal>::min())+1. ||
      exponent>=std::log(std::numeric_limits<CFreal>::max())-1.) return false;
  out=std::pow(T,1.5);return true;
}

inline CFreal sutherlandViscosity(const CFreal dimensionalT, const CFreal referenceT,
                                  const CFreal referenceMu, const CFreal shift)
{
  requirePositiveNormal(dimensionalT,"transport temperature");
  requirePositiveNormal(referenceT,"Sutherland reference temperature");
  requirePositiveNormal(referenceMu,"Sutherland reference viscosity");
  requireTransportShift(shift);
  // Retain the selected law's existing positive-temperature floor.
  const CFreal T=std::max(dimensionalT,1.e-10), limit=std::numeric_limits<CFreal>::max();
  CFreal numerator,quotient,ratio,power,result;
  if (referenceT<limit-shift && T<limit-shift &&
      transportProduct(referenceMu,referenceT+shift,numerator) &&
      transportRatio(numerator,T+shift,quotient) && transportRatio(T,referenceT,ratio) &&
      transportPower(ratio,power) && transportProduct(quotient,power,result)) return result;
  // A power/product may overflow even when the complete law is representable.
  return transportExp(std::log(referenceMu)+transportLogSum(referenceT,shift)-
      transportLogSum(T,shift)+1.5*(std::log(T)-std::log(referenceT)));
}
inline CFreal simplerSutherlandViscosity(const CFreal dimensionalT,
                                         const CFreal coefficient, const CFreal shift)
{
  requirePositiveNormal(dimensionalT,"transport temperature");
  requirePositiveNormal(coefficient,"simpler Sutherland coefficient");
  requireTransportShift(shift);
  const CFreal T=std::max(dimensionalT,1.e-10);
  CFreal power,numerator,result;
  if (T<std::numeric_limits<CFreal>::max()-shift && transportPower(T,power) &&
      transportProduct(coefficient,power,numerator) && transportRatio(numerator,shift+T,result)) return result;
  return transportExp(std::log(coefficient)+1.5*std::log(T)-transportLogSum(shift,T));
}

// Exact-type dispatch keeps unrelated KOmega/NavierStokes users on their
// existing law implementation and refuses unknown overrides before evaluation.
// Fixed/constant laws do not consume p/T scales; Sutherland laws do not consume
// pressure. Do not form or reject unused dimensional products.
inline CFreal laminarViscosity(NSTerm& model, const CFreal p, const CFreal T,
                               const CFreal pressureScale, const CFreal temperatureScale)
{
  const CFreal muScale=model.getReferencePhysicalData()[NSTerm::MU];
  requirePositiveNormal(muScale,"reference viscosity scale");
  const DynamicViscosityLaw& law=model.getDynamicViscosityLaw();
  CFreal mu;
  if (typeid(law)==typeid(SutherlandDynViscosity)) {
    const auto& s=static_cast<const SutherlandDynViscosity&>(law);
    mu=sutherlandViscosity(positiveProduct(temperatureScale,T,"dimensional transport temperature"),
        s.getReferenceTemperature(),s.getReferenceViscosity(),s.getSutherlandConstant());
  } else if (typeid(law)==typeid(SimplerSutherlandDynViscosity)) {
    const auto& s=static_cast<const SimplerSutherlandDynViscosity&>(law);
    mu=simplerSutherlandViscosity(positiveProduct(temperatureScale,T,"dimensional transport temperature"),
        s.getViscosityCoefficient(),s.getTemperatureShift());
  } else if (typeid(law)==typeid(FixedDynViscosity)) {
    mu=static_cast<const FixedDynViscosity&>(law).getFixedViscosity();
  } else if (typeid(law)==typeid(ConstantDynViscosity)) {
    mu=static_cast<const ConstantDynViscosity&>(law).getConstantViscosity();
  } else if (typeid(law)==typeid(FixedKinViscosity)) {
    const CFreal pdim=positiveProduct(pressureScale,p,"dimensional transport pressure");
    const CFreal Tdim=positiveProduct(temperatureScale,T,"dimensional transport temperature");
    mu=positiveProduct(positiveRatio(pdim,positiveProduct(287.046,Tdim,"fixed-kinematic R*T"),
        "fixed-kinematic density"),static_cast<const FixedKinViscosity&>(law).getKinematicViscosity(),
        "fixed-kinematic dynamic viscosity");
  } else {
    throw Common::BadValueException(FromHere(), "Checked transition transport is not implemented for this viscosity law");
  }
  return positiveRatio(mu,muScale,"nondimensional laminar viscosity");
}

} } } }
#endif
