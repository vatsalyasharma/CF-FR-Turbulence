// Copyright (C) 2026 COOLFluiD contributors; LGPLv3, see doc/lgpl.txt.
#ifndef COOLFluiD_GReKO_LM2009Source_hh
#define COOLFluiD_GReKO_LM2009Source_hh
#include <array>
#include <algorithm>
#include <cmath>
#include "Common/COOLFluiD.hh"
#include "Common/BadValueException.hh"
#include "NavierStokes/TransitionArithmetic.hh"

namespace COOLFluiD { namespace Physics { namespace GReKO { namespace lm2009 {
// State-local forward differentiation. Spatial gradients are held fixed, as in
// StdSourceTerm's numerical Jacobian. Nine entries cover both primitive layouts.
struct Scalar {
  CFreal v; std::array<CFreal,9> d;
  Scalar(CFreal value=0.):v(value),d{} {}
  static Scalar variable(CFreal value,unsigned i) {Scalar x(value);x.d[i]=1;return x;}
};
inline Scalar operator+(const Scalar& a,const Scalar& b) {Scalar c(a.v+b.v);for(unsigned i=0;i<9;++i)c.d[i]=a.d[i]+b.d[i];return c;}
inline Scalar operator-(const Scalar& a,const Scalar& b) {Scalar c(a.v-b.v);for(unsigned i=0;i<9;++i)c.d[i]=a.d[i]-b.d[i];return c;}
inline Scalar operator-(const Scalar& a) {return Scalar(0)-a;}
inline Scalar operator*(const Scalar& a,const Scalar& b) {Scalar c(a.v*b.v);for(unsigned i=0;i<9;++i)c.d[i]=a.d[i]*b.v+a.v*b.d[i];return c;}
inline Scalar operator/(const Scalar& a,const Scalar& b) {Scalar c(a.v/b.v);for(unsigned i=0;i<9;++i)c.d[i]=(a.d[i]-c.v*b.d[i])/b.v;return c;}
inline Scalar minimum(const Scalar& a,const Scalar& b) {return b.v<a.v?b:a;}
inline Scalar maximum(const Scalar& a,const Scalar& b) {return a.v<b.v?b:a;}
inline Scalar power(const Scalar& a,CFreal p) {Scalar b(std::pow(a.v,p));CFreal slope=(a.v!=0 || p>=1)?p*std::pow(a.v,p-1):0;for(unsigned i=0;i<9;++i)b.d[i]=slope*a.d[i];return b;}
inline Scalar exponential(const Scalar& a) {Scalar b(std::exp(a.v));for(unsigned i=0;i<9;++i)b.d[i]=b.v*a.d[i];return b;}
inline Scalar expNegativePower(const Scalar& a,unsigned p) {
  // Saturate only where exp(-a^p) and its derivative round to zero.
  if(a.v>=std::pow(800.,1./p))return Scalar(0);
  return exponential(-power(a,p));
}
inline Scalar tanhPower(const Scalar& a,unsigned p) {
  if(a.v>=(p==2?8.:4.))return Scalar(1);
  Scalar z=power(a,p),b(std::tanh(z.v));for(unsigned i=0;i<9;++i)b.d[i]=(1-b.v*b.v)*z.d[i];return b;
}
inline Scalar cappedRatio(std::initializer_list<Scalar> n,std::initializer_list<Scalar> d,CFreal cap) {
  // Logarithms avoid overflowing intermediate products in saturated ratios.
  CFreal l=0;unsigned zeros=0;
  for(const auto& x:n){if(x.v==0)++zeros;else l+=std::log(x.v);}
  for(const auto& x:d){if(x.v==0)return Scalar(zeros?0:cap);l-=std::log(x.v);}
  if(zeros)return Scalar(0); // active bound convention at zero k
  if(l>=std::log(cap))return Scalar(cap);
  Scalar r(std::exp(l));
  for(unsigned i=0;i<9;++i){CFreal g=0;for(const auto& x:n)g+=x.d[i]/x.v;for(const auto& x:d)g-=x.d[i]/x.v;r.d[i]=r.v*g;}
  return r;
}
inline Scalar speed(const Scalar* state,unsigned dim) {
  CFreal scale=0;for(unsigned i=1;i<=dim;++i)scale=std::max(scale,std::abs(state[i].v));
  if(scale==0)return Scalar(0);
  Scalar sum;for(unsigned i=1;i<=dim;++i)sum=sum+(state[i]/scale)*(state[i]/scale);
  return scale*power(sum,.5);
}
inline Scalar rawEquilibrium(const Scalar& tu) {
  return tu.v<=1.3?1173.51-589.428*tu+0.2196/(tu*tu):331.50*power(tu-.5658,-.671);
}
inline Scalar pressureFactor(const Scalar& lambda,const Scalar& tu) {
  if(lambda.v<=0)return 1+(12.986*lambda+123.66*lambda*lambda+405.689*lambda*lambda*lambda)*exponential(-power(tu/1.5,1.5));
  return 1+.275*(1-exponential(-35.*lambda))*exponential(-tu/.5);
}
inline Scalar equilibrium(const Scalar& tu,const Scalar& lambdaCoef,bool pressureGradient) {
  const Scalar raw=rawEquilibrium(tu);
  if(!pressureGradient)return maximum(20,raw);
  // Solve Re = max(20, raw(Tu)*F(clamp(lambdaCoef*Re^2))). This has a
  // bracket because 0 < F <= 1.275 over the specified lambda interval.
  CFreal lo=20,hi=std::max(20.,1.275*raw.v);
  for(unsigned j=0;j<64;++j){
    CFreal mid=(lo+hi)/2;
    Scalar lambda=minimum(.1,maximum(-.1,lambdaCoef.v*mid*mid));
    CFreal target=std::max(20.,(raw.v*pressureFactor(lambda,tu.v)).v);
    if(mid<target)lo=mid;else hi=mid;
  }
  Scalar re((lo+hi)/2);
  Scalar lambda=minimum(.1,maximum(-.1,lambdaCoef*re.v*re.v));
  Scalar target=maximum(20,raw*pressureFactor(lambda,tu));
  // Implicit-function derivative: dRe/dq = target_q/(1-target_Re).
  Scalar rvar=Scalar::variable(re.v,0);
  Scalar lr=minimum(.1,maximum(-.1,lambdaCoef.v*rvar*rvar));
  Scalar tr=maximum(20,raw.v*pressureFactor(lr,tu.v));
  const CFreal den=1-tr.d[0];
  if(!(std::isfinite(den) && den>0))throw Common::BadValueException(FromHere(),"LM2009 equilibrium root has invalid slope");
  for(unsigned i=0;i<9;++i)re.d[i]=target.d[i]/den;
  return re;
}
struct Inputs {
  unsigned dim; enum Model {Base,BSL,SST} model;
  CFreal state[9],grad[9][3],distance,R,isentropicGamma;
  CFreal mu,muP,muT,sigmaK1,sigmaK2,sigmaW1,sigmaW2,beta1,beta2,betaStar,gamma1,gamma2,coeffTau;
  bool limitProduction,pressureGradient,decouple,limitReProduction,vorticityProduction,neglectDilatation;
  Inputs():dim(2),model(SST),state{},grad{},distance(0),R(0),isentropicGamma(0),mu(0),muP(0),muT(0),sigmaK1(.85),sigmaK2(1),sigmaW1(.5),sigmaW2(.856),beta1(.075),beta2(.0828),betaStar(.09),gamma1(5./9.),gamma2(.44),coeffTau(1),limitProduction(true),pressureGradient(true),decouple(false),limitReProduction(false),vorticityProduction(false),neglectDilatation(false){}
};
struct Result {std::array<Scalar,9> source;Scalar gammaEff,gammaSep,onset,length,equilibriumRe,mut;};
inline Result evaluate(const Inputs& p) {
  const unsigned kID=p.dim+2,wID=kID+1,gID=kID+2,rID=kID+3,tID=p.dim+1;
  Scalar q[9];for(unsigned i=0;i<p.dim+6;++i)q[i]=Scalar::variable(p.state[i],i);
  const Scalar rho=q[0]/(p.R*q[tID]),k=maximum(0,q[kID]),omega=exponential(q[wID]);
  // Do not cap gamma at one: its signed source restores small overshoots.
  const Scalar ga=maximum(0,q[gID]),re=maximum(20,q[rID]);
  Scalar mu(p.mu);mu.d[0]=p.muP;mu.d[tID]=p.muT;
  CFreal s2=0,v2=0,div=0,dotKW=0,gradW2=0;
  for(unsigned i=0;i<p.dim;++i){
    div+=p.grad[1+i][i];dotKW+=p.grad[kID][i]*p.grad[wID][i];gradW2+=p.grad[wID][i]*p.grad[wID][i];
    for(unsigned j=0;j<p.dim;++j){CFreal s=.5*(p.grad[1+i][j]+p.grad[1+j][i]),v=.5*(p.grad[1+i][j]-p.grad[1+j][i]);s2+=2*s*s;v2+=2*v*v;}
  }
  const CFreal strain=std::sqrt(s2),vorticity=std::sqrt(v2),d=std::max(p.distance,1.e-12);
  Scalar F1(1),viscosityRatio=omega;
  if(p.model!=Inputs::Base){
    Scalar cd=maximum(1.e-10,2*rho*p.sigmaW2*dotKW);
    Scalar a=cappedRatio({power(k,.5)},{.09,omega,d},4),b=cappedRatio({500.,mu},{rho,d,d,omega},4),c=cappedRatio({4.,rho,p.sigmaW2,k},{cd,d,d},4);
    Scalar ry=cappedRatio({rho,d,power(k,.5)},{mu,120.},3);
    F1=maximum(tanhPower(minimum(maximum(a,b),c),4),expNegativePower(ry,8));
  }
  if(p.model==Inputs::SST){
    Scalar a=cappedRatio({2.,power(k,.5)},{.09,omega,d},8),b=cappedRatio({500.,mu},{rho,d,d,omega},8);
    viscosityRatio=maximum(.31*omega,strain*tanhPower(maximum(a,b),2))/.31;
  }
  const Scalar mut=rho*k/viscosityRatio;
  const Scalar beta=F1*p.beta1+(1-F1)*p.beta2,gamma=F1*p.gamma1+(1-F1)*p.gamma2,sigmaW=F1*p.sigmaW1+(1-F1)*p.sigmaW2;
  const Scalar U=speed(q,p.dim);
  if(U.v==0 && vorticity==0)throw Common::BadValueException(FromHere(),"LM2009 wake limit at joint zero speed/vorticity is path dependent");
  Scalar tu=U.v>0?maximum(.027,100*power((2./3.)*k,.5)/U):Scalar(k.v>0?1.e100:.027);
  Scalar lambdaCoef;
  if(p.pressureGradient && U.v>0){
    Scalar dUds;for(unsigned i=0;i<p.dim;++i)for(unsigned j=0;j<p.dim;++j)dUds=dUds+(q[i+1]/U)*(q[j+1]/U)*p.grad[i+1][j];
    lambdaCoef=(mu/rho/U/U)*dUds;
  }
  const Scalar eq=equilibrium(tu,lambdaCoef,p.pressureGradient);
  Scalar rec=re.v<=1870?-3.96035+1.0120656*re-.00086823*re*re+6.96506e-7*re*re*re-1.74105e-10*re*re*re*re:re-(593.11+.482*(re-1870));
  Scalar length;
  if(re.v<400)length=39.8189-.011927*re-.000132567*re*re;
  else if(re.v<596)length=263.404-1.23939*re+.00194548*re*re-1.01695e-6*re*re*re;
  else if(re.v<1200)length=.5-.0003*(re-596);
  else length=.3188;
  Scalar rew=rho*d*d*omega/mu,rt=rho*k/mu/omega,rev=rho*d*d*strain/mu;
  Scalar wakeRatio;
  if(vorticity>0)wakeRatio=cappedRatio({rho,U,U},{375.,vorticity,mu,re},6.);
  else wakeRatio=Scalar(U.v>0?6.:0.);
  const Scalar ftheta=minimum(1,maximum(expNegativePower(rew/1.e5,2)*expNegativePower(wakeRatio,4),1-power((50*ga-1)/49,2)));
  Scalar sep=minimum(2,2*maximum(0,rev/(3.235*rec)-1)*expNegativePower(rt/20,4))*ftheta;
  Scalar eff=p.decouple?Scalar(1):maximum(ga,sep),dkFactor=p.decouple?Scalar(1):minimum(1,maximum(.1,eff));
  Scalar fon1=rev/(2.193*rec),fon2=fon1.v>=2?Scalar(2):minimum(2,maximum(fon1,power(fon1,4)));
  Scalar fon3=rt.v>=2.5?Scalar(0):maximum(0,1-power(rt/2.5,3)),onset=maximum(0,fon2-fon3);
  Scalar sub=expNegativePower(rew/200,2);length=length*(1-sub)+40*sub;
  const CFreal shear=p.coeffTau*(p.vorticityProduction?v2:s2-(2./3.)*div*div);
  const CFreal dilatation=p.vorticityProduction && p.neglectDilatation?0:(2./3.)*div;
  Scalar production=mut*shear-rho*k*dilatation;
  Scalar productionW=gamma*rho*(shear-viscosityRatio*dilatation)/omega;
  // Apply the underlying SST limiter before transition modifies k alone.
  const Scalar destruction=rho*omega*k*p.betaStar;
  if(p.limitProduction){production=minimum(10*destruction,production);productionW=minimum(10*rho*viscosityRatio*p.betaStar*gamma,productionW);}
  Result out;out.mut=mut;out.gammaEff=eff;out.gammaSep=sep;out.onset=onset;out.length=length;out.equilibriumRe=eq;
  out.source[kID]=eff*maximum(0,production)-dkFactor*destruction;
  // Cross diffusion is signed. The log transform adds D*|grad(log omega)|^2.
  out.source[wID]=productionW-rho*omega*beta+(1-F1)*2*rho/omega*p.sigmaW2*dotKW+p.coeffTau*(mu+mut*sigmaW)*gradW2;
  out.source[gID]=length*2*rho*strain*power(ga*onset,.5)*(1-ga)-.06*rho*vorticity*ga*expNegativePower(rt/4,4)*(50*ga-1);
  out.source[rID]=.03*rho*(rho*U*U/(500*mu))*(eq-q[rID])*(1-ftheta);
  if(p.limitReProduction)out.source[rID]=maximum(0,out.source[rID]);
  out.source[tID]=-out.source[kID];
  return out;
}
}}}}
#endif
