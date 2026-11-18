// Copyright (C) 2026 COOLFluiD contributors; LGPLv3.
#ifndef COOLFluiD_FR_GReKODiffusiveUpdate_hh
#define COOLFluiD_FR_GReKODiffusiveUpdate_hh
#include "Framework/PhysicalModel.hh"
#include "NavierStokes/NavierStokesVarSet.hh"
#include "NavierStokes/EulerTerm.hh"
namespace COOLFluiD { namespace FluxReconstructionMethod { namespace detail {
inline bool isGReKLogO() {
  const std::string name=Framework::PhysicalModelStack::getActive()->getDiffusiveName();
  return name=="NavierStokes2DGReKLogO" || name=="NavierStokes3DGReKLogO";
}
inline CFreal greDiffusivity(Physics::NavierStokes::NavierStokesVarSet& var,
                            const RealVector& state,const std::vector<RealVector*>& gradients) {
  const CFreal rho=var.getDensity(state);
  const CFreal totalMu=var.getDynViscosity(state,gradients);
  const auto euler=Framework::PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm().d_castTo<Physics::NavierStokes::EulerTerm>();
  const CFreal cv=euler->getR()/(euler->getGamma()-1.);
  const CFreal thermal=var.getModel().getCoeffQ()*var.getThermConductivity(state,totalMu*std::max(1.,var.getModel().getPrOverPrT()))/cv;
  // ReTheta diffusion dominates the k, logOmega, gamma and 4/3 stress
  // coefficients for the supported base/BSL/SST gamma–ReTheta models.
  return std::max(2*var.getModel().getCoeffTau()*totalMu,thermal)/rho;
}
}}}
#endif
