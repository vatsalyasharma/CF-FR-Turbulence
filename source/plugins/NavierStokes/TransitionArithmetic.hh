#ifndef COOLFluiD_Physics_NavierStokes_TransitionArithmetic_hh
#define COOLFluiD_Physics_NavierStokes_TransitionArithmetic_hh

#include <cmath>
#include <algorithm>
#include <limits>
#include <string>
#include <initializer_list>
#include "Common/BadValueException.hh"

namespace COOLFluiD { namespace Physics { namespace NavierStokes {
namespace transition_detail {

inline void requirePositiveNormal(const CFreal x, const char* quantity)
{
  if (!(std::isnormal(x) && x > 0.))
    throw Common::BadValueException(FromHere(), std::string("Transition thermodynamic domain: positive normal ") + quantity + " required");
}

// Check before multiplication/division, including underflow of denominators.
// Equality at rounded range boundaries is conservatively rejected.
// These checks describe the implemented arithmetic domain, not a physical floor.
inline CFreal positiveProduct(const CFreal a, const CFreal b, const char* quantity)
{
  requirePositiveNormal(a, quantity); requirePositiveNormal(b, quantity);
  if ((b > 1. && a >= std::numeric_limits<CFreal>::max()/b) ||
      (b < 1. && a <= std::numeric_limits<CFreal>::min()/b))
    throw Common::BadValueException(FromHere(), std::string("Transition thermodynamic product out of range: ") + quantity);
  return a*b;
}

inline CFreal positiveRatio(const CFreal a, const CFreal b, const char* quantity)
{
  requirePositiveNormal(a, quantity); requirePositiveNormal(b, quantity);
  if ((b < 1. && a >= std::numeric_limits<CFreal>::max()*b) ||
      (b > 1. && a <= std::numeric_limits<CFreal>::min()*b))
    throw Common::BadValueException(FromHere(), std::string("Transition thermodynamic ratio out of range: ") + quantity);
  return a/b;
}

// Guard the actual EOS/enthalpy/sound-speed intermediates and the source's
// temperature denominator. Leave ordinary evaluation order unchanged in callers.
// This deliberately does not impose volume-source k bounds on reflected ghosts.
inline void validateThermodynamics(const CFreal p, const CFreal T,
                                   const CFreal R, const CFreal gamma)
{
  requirePositiveNormal(p, "pressure"); requirePositiveNormal(T, "temperature");
  requirePositiveNormal(R, "gas constant");
  if (!(std::isfinite(gamma) && gamma > 1.))
    throw Common::BadValueException(FromHere(), "Transition thermodynamic domain: finite gamma > 1 required");
  const CFreal RT = positiveProduct(R, T, "R*T");
  const CFreal rho = positiveRatio(p, RT, "density");
  positiveRatio(1., rho, "inverse density");
  positiveProduct(RT, T, "R*T*T");
  const CFreal gammaRatio = positiveRatio(gamma, gamma-1., "gamma/(gamma-1)");
  const CFreal hPressure = positiveProduct(gammaRatio, p, "gamma/(gamma-1)*p");
  positiveRatio(hPressure, rho, "static enthalpy");
  positiveRatio(positiveProduct(gamma, p, "gamma*p"), rho, "sound speed squared");
}

// Signed arithmetic for gradient intermediates. Zero/subnormal results are
// allowed. Reject before overflow; rounded upper endpoints are conservative.
inline CFreal finiteProduct(const CFreal a, const CFreal b)
{
  const CFreal x = std::abs(a), y = std::abs(b);
  if (!std::isfinite(a) || !std::isfinite(b) ||
      (y > 1. && x != 1. && x >= std::numeric_limits<CFreal>::max()/y))
    throw Common::BadValueException(FromHere(), "Transition gradient product out of range");
  return a*b;
}
inline CFreal finiteSum(const CFreal a, const CFreal b)
{
  const CFreal limit = std::numeric_limits<CFreal>::max();
  if (!std::isfinite(a) || !std::isfinite(b) ||
      (b > 0. && a > limit-b) || (b < 0. && a < -limit-b))
    throw Common::BadValueException(FromHere(), "Transition gradient sum out of range");
  return a+b;
}
inline CFreal finiteRatio(const CFreal a, const CFreal b)
{
  if (!std::isfinite(a) || !std::isfinite(b) || b == 0. ||
      (std::abs(b) < 1. && std::abs(a) >= std::numeric_limits<CFreal>::max()*std::abs(b)))
    throw Common::BadValueException(FromHere(), "Transition gradient ratio out of range");
  return a/b;
}

// Preserve the historical operation order for ordinary gradients. At extreme
// magnitudes check tensor sums, normalize before squares, and check rescaling.
// The fast-path threshold bounds the largest 3D sum by 18*scale^2 < max.
inline CFreal velocityGradientNorm(const CFreal g[3][3], const unsigned dim,
                                   const bool vorticity, const CFreal axis = 0.)
{
  CFreal scale = std::abs(axis);
  if (!std::isfinite(axis))
    throw Common::BadValueException(FromHere(), "Transition requires finite axisymmetric velocity/radius");
  for (unsigned i=0; i<dim; ++i) for (unsigned j=0; j<dim; ++j) {
    if (!std::isfinite(g[i][j]))
      throw Common::BadValueException(FromHere(), "Transition requires finite velocity gradients");
    scale = std::max(scale,std::abs(g[i][j]));
  }
  if (scale == 0.) return 0.;
  const bool scaled = scale > std::sqrt(std::numeric_limits<CFreal>::max()/32.) ||
                      scale < 4.*std::sqrt(std::numeric_limits<CFreal>::min());
  if (scaled) {
    // Form each tensor component before choosing the norm scale. Scaling raw
    // gradients first can erase small nonzero components when large opposing
    // gradients cancel (or when an unused diagonal dominates vorticity).
    CFreal terms[7] = {}, weights[7] = {};
    unsigned count = 0;
    if (!vorticity) {
      for (unsigned i=0; i<dim; ++i) { terms[count]=g[i][i]; weights[count++]=2.; }
      if (dim == 2) { terms[count]=axis; weights[count++]=2.; }
    }
    for (unsigned i=0; i<dim; ++i) for (unsigned j=i+1; j<dim; ++j) {
      terms[count]=finiteSum(g[i][j],vorticity ? -g[j][i] : g[j][i]);
      weights[count++]=1.;
    }
    scale = 0.;
    for (unsigned i=0; i<count; ++i) scale=std::max(scale,std::abs(terms[i]));
    if (scale == 0.) return 0.;
    CFreal sum=0.;
    for (unsigned i=0; i<count; ++i) {
      const CFreal normalized=terms[i]/scale;
      sum += weights[i]*normalized*normalized;
    }
    return finiteProduct(std::sqrt(sum),scale);
  }
  const CFreal (*h)[3] = g;
  const CFreal a = axis;
  CFreal result;
  if (vorticity) {
    const CFreal xy = h[0][1]-h[1][0];
    if (dim == 2) result = std::abs(xy);
    else {
      const CFreal xz = h[0][2]-h[2][0], yz = h[2][1]-h[1][2];
      result = std::sqrt(xy*xy+xz*xz+yz*yz);
    }
  } else {
    const CFreal xy = h[0][1]+h[1][0];
    if (dim == 2)
      result = std::sqrt(2.*(h[0][0]*h[0][0]+0.5*xy*xy+h[1][1]*h[1][1]+a*a));
    else {
      const CFreal xz = h[0][2]+h[2][0], yz = h[1][2]+h[2][1];
      result = std::sqrt(2.*(h[0][0]*h[0][0]+h[1][1]*h[1][1]+h[2][2]*h[2][2])+xy*xy+xz*xz+yz*yz);
    }
  }
  return result;
}

// Preflight the uncancelled velocity-gradient quadratics used by the volume
// sources and their legacy analytical matrices. This deliberately requires
// both strain/vorticity squares and both legacy production forms to fit even
// if a selected option/cap would later discard one. No gradient is clipped.
// Coefficient/state products and turbulence-gradient products need separate
// checks; passing this guard is NOT a guarantee of the complete source domain.
inline void validateVelocityGradientQuadratics(const CFreal g[3][3],
                                                const unsigned dim, const CFreal axis)
{
  const CFreal strain = velocityGradientNorm(g,dim,false,axis);
  const CFreal rotation = velocityGradientNorm(g,dim,true);
  finiteProduct(strain,strain); finiteProduct(rotation,rotation);
  const CFreal difference = finiteSum(g[0][0],-g[1][1]);
  const CFreal shear = finiteSum(g[0][1],g[1][0]);
  const CFreal diagonal = finiteSum(g[0][0],g[1][1]);
  const CFreal base = finiteSum(finiteProduct(difference,difference),finiteProduct(g[0][0],g[1][1]));
  const CFreal cylindrical = finiteProduct(finiteSum(diagonal,-axis),axis);
  const CFreal shear2 = finiteProduct(shear,shear);
  finiteSum(finiteProduct(4./3.,finiteSum(base,-cylindrical)),shear2);
  finiteSum(finiteSum(finiteProduct(4./3.,base),shear2),-cylindrical);
  finiteSum(diagonal,dim == 3 ? g[2][2] : axis);
}

// Positive powers used by the GA freestream reconstruction. Keep one natural-
// log unit inside the normal range to avoid rounded libm endpoint ambiguity.
// This is an explicit arithmetic support restriction, not physical clipping.
inline CFreal positivePower(const CFreal base, const CFreal exponent, const char* quantity)
{
  requirePositiveNormal(base, quantity); requirePositiveNormal(exponent, quantity);
  const CFreal logarithm = std::log(base);
  if ((logarithm > 0. && exponent >= (std::log(std::numeric_limits<CFreal>::max())-1.)/logarithm) ||
      (logarithm < 0. && exponent >= (std::log(std::numeric_limits<CFreal>::min())+1.)/logarithm))
    throw Common::BadValueException(FromHere(), std::string("Gamma-alpha freestream power out of range: ") + quantity);
  return std::pow(base, exponent);
}

struct GammaAlphaFreestream {
  CFreal density, mach, speed, temperature;
};

// Preserve the existing isentropic reconstruction and its operation order.
// Reject unsupported inputs/intermediates before division, powers or output
// mutation. This does not repair the freestream relation itself (F14), nor
// validate the later viscosity/correlation products (remaining F08 work).
inline GammaAlphaFreestream gammaAlphaFreestream(const CFreal p, const CFreal R,
                                                 const CFreal gamma, const CFreal machInf,
                                                 const CFreal uInf, const CFreal pInf)
{
  requirePositiveNormal(p, "GA local pressure");
  requirePositiveNormal(R, "GA gas constant");
  requirePositiveNormal(machInf, "GA machInf");
  requirePositiveNormal(uInf, "GA uInf");
  requirePositiveNormal(pInf, "GA pInf");
  if (!(std::isfinite(gamma) && gamma > 1.))
    throw Common::BadValueException(FromHere(), "Gamma-alpha freestream requires finite gamma > 1");
  const CFreal gm = gamma-1.;
  const CFreal rhoInf = positiveProduct(positiveProduct(positiveRatio(
      positiveProduct(machInf,machInf,"GA machInf squared"),
      positiveProduct(uInf,uInf,"GA uInf squared"),"GA Mach/velocity squares"),
      gamma,"GA freestream density factor"),pInf,"GA freestream density");
  const CFreal term = positiveProduct(positiveProduct(positiveRatio(gm,2.,"GA isentropic factor"),
      machInf,"GA stagnation factor"),machInf,"GA stagnation factor");
  if (term >= std::numeric_limits<CFreal>::max()-1.)
    throw Common::BadValueException(FromHere(), "Gamma-alpha freestream stagnation factor out of range");
  const CFreal pTotTerm = 1.+term;
  const CFreal exponent = positiveRatio(gamma,gm,"GA stagnation exponent");
  const CFreal pTotal = positiveProduct(positivePower(pTotTerm,exponent,"GA stagnation pressure power"),p,"GA stagnation pressure");
  GammaAlphaFreestream result;
  result.density = positiveProduct(rhoInf,positivePower(positiveRatio(p,pInf,"GA local/freestream pressure"),
      positiveRatio(1.,gamma,"GA inverse gamma"),"GA density power"),"GA local freestream density");
  const CFreal machTerm = positivePower(positiveRatio(pTotal,p,"GA total/static pressure"),
      positiveRatio(1.,exponent,"GA inverse stagnation exponent"),"GA local Mach power")-1.;
  // Small positive machInf can round the reconstructed Mach term to zero.
  // The downstream source assumes positive freestream speed; do not invent it.
  result.mach = std::sqrt(positiveProduct(positiveRatio(2.,gm,"GA Mach factor"),machTerm,"GA local Mach squared"));
  result.speed = positiveProduct(std::sqrt(positiveRatio(positiveProduct(gamma,p,"GA gamma*p"),
      result.density,"GA local sound speed squared")),result.mach,"GA local freestream speed");
  result.temperature = positiveRatio(p,positiveProduct(result.density,R,"GA density*R"),"GA local freestream temperature");
  validateThermodynamics(p,result.temperature,R,gamma);
  return result;
}

// Gamma-alpha local freestream after Jacobs, Van den Eynde & Steelant (2025), Appendix B: the edge Mach number from
// the local static pressure and the INLET total pressure (B1), the edge temperature from the inlet total temperature
// (B2), c_e = M_e sqrt(kappa R T_e) (B3). Density and viscosity use the inlet reference (Sect. 2.1). The legacy
// gammaAlphaFreestream() forms the total pressure from the local pressure, so its Mach is always machInf.
inline GammaAlphaFreestream gammaAlphaFreestreamInlet(const CFreal p, const CFreal R,
                                                      const CFreal gamma, const CFreal machInf,
                                                      const CFreal uInf, const CFreal pInf)
{
  requirePositiveNormal(p, "GA local pressure");
  requirePositiveNormal(R, "GA gas constant");
  requirePositiveNormal(machInf, "GA machInf");
  requirePositiveNormal(uInf, "GA uInf");
  requirePositiveNormal(pInf, "GA pInf");
  if (!(std::isfinite(gamma) && gamma > 1.))
    throw Common::BadValueException(FromHere(), "Gamma-alpha freestream requires finite gamma > 1");
  const CFreal gm = gamma-1.;
  const CFreal rhoInf = gamma*pInf*machInf*machInf/(uInf*uInf);
  const CFreal TInf = pInf/(rhoInf*R);
  const CFreal stag = 1. + 0.5*gm*machInf*machInf;
  const CFreal exponent = gamma/gm;
  const CFreal p0 = pInf*std::pow(stag, exponent);
  const CFreal T0 = TInf*stag;
  // p >= p0 (stagnation) gives no positive edge Mach; keep a tiny positive value as the legacy routine requires.
  const CFreal machTerm = std::max(std::pow(p0/p, 1./exponent) - 1., 1.e-24);
  GammaAlphaFreestream result;
  result.mach = std::sqrt(2./gm*machTerm);
  result.temperature = T0/(1. + 0.5*gm*result.mach*result.mach);
  result.speed = result.mach*std::sqrt(gamma*R*result.temperature);
  result.density = rhoInf;
  validateThermodynamics(p,result.temperature,R,gamma);
  return result;
}

// Evaluate min(product(numerator)/product(denominator), cap) without
// forming an overflowing ratio or an underflowed denominator. All factors
// must be finite; numerator factors may be zero, denominators must be positive.
// The fast path preserves the original left-to-right arithmetic. Only extreme
// products use normalized binary mantissas/exponents (no wider type required).
// Use only where the caller is already insensitive to values above cap.
inline CFreal cappedPositiveProductRatio(std::initializer_list<CFreal> numerator,
                                         std::initializer_list<CFreal> denominator,
                                         const CFreal cap)
{
  if (!(std::isfinite(cap) && cap > 0.))
    throw Common::BadValueException(FromHere(), "Transition bounded ratio requires finite positive cap");
  bool zero = false;
  for (const CFreal x : numerator) {
    if (!(std::isfinite(x) && x >= 0.))
      throw Common::BadValueException(FromHere(), "Transition bounded ratio requires finite nonnegative numerator factors");
    zero = zero || x == 0.;
  }
  for (const CFreal x : denominator)
    if (!(std::isfinite(x) && x > 0.))
      throw Common::BadValueException(FromHere(), "Transition bounded ratio requires finite positive denominator factors");
  if (zero) return 0.;

  const CFreal largest = std::numeric_limits<CFreal>::max();
  const CFreal smallest = std::numeric_limits<CFreal>::min();
  const auto product = [largest, smallest](std::initializer_list<CFreal> factors, CFreal& result) {
    result = 1.;
    for (const CFreal x : factors) {
      if ((x > 1. && result >= largest/x) ||
          (x < 1. && result <= smallest/x)) return false;
      result *= x;
    }
    return true;
  };
  CFreal n, d;
  if (product(numerator, n) && product(denominator, d) &&
      !(d < 1. && n >= largest*d) && !(d > 1. && n <= smallest*d))
    return std::min(n/d, cap);

  const auto scaledProduct = [](std::initializer_list<CFreal> factors, int& exponent) {
    CFreal mantissa = 1.;
    exponent = 0;
    for (const CFreal x : factors) {
      int e, shift;
      const CFreal part = std::frexp(x, &e);
      mantissa = std::frexp(mantissa*part, &shift);
      exponent += e + shift;
    }
    return mantissa;
  };
  int ne, de, shift, ce;
  n = scaledProduct(numerator, ne);
  d = scaledProduct(denominator, de);
  const CFreal m = std::frexp(n/d, &shift);
  const int exponent = ne-de+shift;
  const CFreal cm = std::frexp(cap, &ce);
  if (exponent > ce || (exponent == ce && m >= cm)) return cap;
  return std::min(std::ldexp(m, exponent), cap);
}

inline CFreal expNegativeSquare(const CFreal x)
{
  static const CFreal cutoff = std::sqrt(-std::log((std::numeric_limits<CFreal>::denorm_min() > 0. ? std::numeric_limits<CFreal>::denorm_min() : std::numeric_limits<CFreal>::min()))) + 1.;
  return std::abs(x) >= cutoff ? 0. : std::exp(-x*x);
}
inline CFreal expNegativeFourth(const CFreal x)
{
  static const CFreal cutoff = std::sqrt(std::sqrt(-std::log((std::numeric_limits<CFreal>::denorm_min() > 0. ? std::numeric_limits<CFreal>::denorm_min() : std::numeric_limits<CFreal>::min())))) + 1.;
  return std::abs(x) >= cutoff ? 0. : std::exp(-std::pow(x,4));
}
inline CFreal tanhSquare(const CFreal x)
{
  return std::abs(x) >= 8. ? 1. : std::tanh(x*x);
}
inline CFreal tanhFourth(const CFreal x)
{
  return std::abs(x) >= 4. ? 1. : std::tanh(x*x*x*x);
}

} } } }
#endif
