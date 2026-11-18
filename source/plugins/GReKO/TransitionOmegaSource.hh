#ifndef COOLFluiD_GReKO_TransitionOmegaSource_hh
#define COOLFluiD_GReKO_TransitionOmegaSource_hh

#include <cmath>
#include "NavierStokes/TransitionArithmetic.hh"

namespace COOLFluiD { namespace Physics { namespace GReKO { namespace omega_detail {

// Local forward chain rule in k and log(omega), with fixed spatial gradients.
// This is confined to the GR log-omega source, not a numerical perturbation.
struct Partial {
  CFreal value, k, logOmega;
  Partial(CFreal v=0., CFreal dk=0., CFreal dw=0.):value(v),k(dk),logOmega(dw) {}
};
inline Partial operator+(Partial a, Partial b) {return {a.value+b.value,a.k+b.k,a.logOmega+b.logOmega};}
inline Partial operator-(Partial a, Partial b) {return {a.value-b.value,a.k-b.k,a.logOmega-b.logOmega};}
inline Partial operator-(Partial a) {return {-a.value,-a.k,-a.logOmega};}
inline Partial operator*(Partial a, Partial b) {return {a.value*b.value,a.k*b.value+a.value*b.k,a.logOmega*b.value+a.value*b.logOmega};}
inline Partial operator/(Partial a, Partial b) {
  const CFreal v=a.value/b.value;
  return {v,(a.k-v*b.k)/b.value,(a.logOmega-v*b.logOmega)/b.value};
}
// Match std::min/max's first-argument choice on a branch tie.
inline Partial minimum(Partial a,Partial b) {return b.value<a.value?b:a;}
inline Partial maximum(Partial a,Partial b) {return a.value<b.value?b:a;}
inline Partial root(Partial a) {
  const CFreal v=std::sqrt(a.value);
  return v>0.?Partial(v,a.k/(2.*v),a.logOmega/(2.*v)):Partial(0.);
}
inline Partial tanhPower(Partial a,unsigned power) {
  if(a.value>=(power==2?8.:4.))return Partial(1.);
  const CFreal v=std::tanh(std::pow(a.value,power));
  const CFreal slope=power*std::pow(a.value,power-1)*(1.-v*v);
  return {v,slope*a.k,slope*a.logOmega};
}
// The native bounded-ratio evaluator receives separate factors so a tiny
// wall distance or omega cannot underflow a product before saturation.
inline Partial ratioPartial(CFreal value,CFreal cap,CFreal k,CFreal kPower,CFreal omegaPower) {
  if(value>=cap || value==0.)return Partial(value);
  const CFreal dk=kPower!=0. && k>0.?(value/k)*kPower:0.;
  return {value,dk,value*omegaPower};
}
inline Partial expNegativeEighth(Partial a) {
  // Above three the exponential and its derivatives round to zero in CFreal.
  if(a.value>=3.)return Partial(0.);
  const CFreal p=std::pow(a.value,8),v=std::exp(-p);
  const CFreal slope=-8.*std::pow(a.value,7)*v;
  return {v,slope*a.k,slope*a.logOmega};
}

struct Inputs {
  enum Model {Base,BSL,SST} model;
  CFreal rho,mu,k,omega,distance,strain,soundSpeedSquared;
  CFreal gradKDotGradLogOmega,gradLogOmegaSquared;
  CFreal sigma1,sigma2,beta1,beta2,betaStar,gamma1,gamma2,coeffTau;
  CFreal shearProduction,dilatation;
  Partial destructionFactor;
  bool limitProduction;
};

inline Partial source(const Inputs& p) {
  const Partial k(p.k,1.,0.),omega(p.omega,0.,p.omega);
  using NavierStokes::transition_detail::cappedPositiveProductRatio;
  Partial F1(1.),q=omega,mut;
  if(p.model!=Inputs::SST)mut=p.rho*k/omega;
  if(p.model!=Inputs::Base) {
    const CFreal d=std::max(p.distance,1.e-12);
    const CFreal cd=std::max(2.*p.rho*p.sigma2*p.gradKDotGradLogOmega,1.e-10);
    const Partial a=ratioPartial(cappedPositiveProductRatio({std::sqrt(p.k)},{0.09,p.omega,d},4.),4.,p.k,0.5,-1.);
    const Partial b=ratioPartial(cappedPositiveProductRatio({500.,p.mu},{p.rho,d,d,p.omega},4.),4.,p.k,0.,-1.);
    const Partial c=ratioPartial(cappedPositiveProductRatio({4.,p.rho,p.sigma2,p.k},{cd,d,d},4.),4.,p.k,1.,0.);
    const Partial org=tanhPower(minimum(maximum(a,b),c),4);
    const Partial ry=p.rho*d*root(k)/p.mu;
    F1=maximum(org,expNegativeEighth(ry/120.));
  }
  if(p.model==Inputs::SST) {
    const Partial a=ratioPartial(cappedPositiveProductRatio({2.,std::sqrt(p.k)},{0.09,p.omega,p.distance},8.),8.,p.k,0.5,-1.);
    const Partial b=ratioPartial(cappedPositiveProductRatio({500.,p.mu},{p.rho,p.distance,p.distance,p.omega},8.),8.,p.k,0.,-1.);
    const Partial F2=tanhPower(maximum(a,b),2);
    const Partial denominator=maximum(0.31*omega,p.strain*F2);
    q=denominator/0.31;
    mut=0.31*p.rho*k/denominator;
  }
  const Partial sigma=F1*p.sigma1+(1.-F1)*p.sigma2;
  const Partial gamma=F1*p.gamma1+(1.-F1)*p.gamma2;
  const Partial mt=2.*k/p.soundSpeedSquared;
  const Partial fmt=mt.value>0.25?mt*mt-0.25*0.25:Partial(0.);
  const Partial betaStar=p.betaStar*(1.+1.5*fmt);
  const Partial beta=F1*p.beta1+(1.-F1)*p.beta2-p.betaStar*1.5*fmt;
  Partial production=gamma*p.rho*(p.shearProduction-q*p.dilatation)/omega;
  if(p.limitProduction)
    production=minimum(10.*p.rho*q*betaStar*gamma*p.destructionFactor,production);
  production=production+(1.-F1)*2.*p.rho/omega*p.sigma2*p.gradKDotGradLogOmega;
  production=production+p.coeffTau*(p.mu+mut*sigma)*p.gradLogOmegaSquared;
  return maximum(0.,production)+minimum(0.,-p.rho*omega*beta);
}

} } } }
#endif
