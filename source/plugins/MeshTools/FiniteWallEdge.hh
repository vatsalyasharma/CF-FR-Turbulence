// Copyright (C) 2026 COOLFluiD contributors
// Distributed under LGPLv3; see doc/lgpl.txt.
#ifndef COOLFluiD_MeshTools_FiniteWallEdge_hh
#define COOLFluiD_MeshTools_FiniteWallEdge_hh
#include <algorithm>
#include <cmath>

namespace COOLFluiD { namespace MeshTools { namespace detail {

// Closest distance to a finite P1/P2 edge. P2 ordering is endpoints, midpoint
// (LagrangeShapeFunctionLineP2). For a curved edge all stationary points of
// squared distance are considered, together with both endpoints.
template<typename Real>
Real finiteWallEdgeDistance(const Real* p, const Real* n0, const Real* n1,
                            const Real* mid = 0)
{
  Real a[2], b[2], c[2];
  for (unsigned d=0; d<2; ++d) {
    b[d]=(n1[d]-n0[d])/2;
    a[d]=mid ? (n0[d]-mid[d])/2+(n1[d]-mid[d])/2 : 0;
    c[d]=mid ? mid[d]-p[d] : n0[d]-p[d]+b[d];
  }
  const auto distance = [&](Real t) {
    return std::hypot((a[0]*t+b[0])*t+c[0],(a[1]*t+b[1])*t+c[1]);
  };
  Real best=std::min(distance(-1),distance(1));
  Real scale=0;
  for (unsigned d=0; d<2; ++d)
    scale=std::max(scale,std::max(std::abs(a[d]),std::max(std::abs(b[d]),std::abs(c[d]))));
  if (scale==0) return 0;
  // Normalize before forming the cubic coefficients (roots are unchanged).
  for (unsigned d=0; d<2; ++d) { a[d]/=scale; b[d]/=scale; c[d]/=scale; }
  const Real A=2*(a[0]*a[0]+a[1]*a[1]);
  const Real B=3*(a[0]*b[0]+a[1]*b[1]);
  const Real C=b[0]*b[0]+b[1]*b[1]+2*(a[0]*c[0]+a[1]*c[1]);
  const Real D=b[0]*c[0]+b[1]*c[1];
  const auto derivative = [&](Real t) { return ((A*t+B)*t+C)*t+D; };
  const auto scaledDistance = [&](Real t) { return scale*distance(t); };
  if (A==0) {
    if (C>0) best=std::min(best,scaledDistance(std::max(Real(-1),std::min(Real(1),-D/C))));
    return best;
  }
  // Split at critical points of the cubic. Each interval is monotone, so
  // bisection cannot miss a second local minimum as a single Newton solve can.
  Real cuts[4]={-1,1,0,0}; unsigned count=2;
  const Real disc=B*B-3*A*C;
  if (disc>=0) {
    const Real q=-B-std::copysign(std::sqrt(disc),B);
    const Real r0=q/(3*A), r1=q==0 ? r0 : C/q;
    if (r0>-1 && r0<1) cuts[count++]=r0;
    if (r1>-1 && r1<1 && r1!=r0) cuts[count++]=r1;
  }
  std::sort(cuts,cuts+count);
  for (unsigned j=0; j<count; ++j) best=std::min(best,scaledDistance(cuts[j]));
  for (unsigned j=1; j<count; ++j) {
    Real lo=cuts[j-1],hi=cuts[j],fl=derivative(lo),fh=derivative(hi);
    if ((fl<0 && fh>0) || (fl>0 && fh<0)) {
      for (unsigned k=0; k<64; ++k) {
        const Real t=(lo+hi)/2, f=derivative(t);
        if (t==lo || t==hi || f==0) { lo=hi=t; break; }
        if ((fl<0 && f<0) || (fl>0 && f>0)) { lo=t; fl=f; }
        else hi=t;
      }
      best=std::min(best,scaledDistance((lo+hi)/2));
    }
  }
  return best;
}
}}}
#endif
