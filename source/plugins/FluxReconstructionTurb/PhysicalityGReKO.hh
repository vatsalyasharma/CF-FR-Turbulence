// Copyright (C) 2026 COOLFluiD contributors; LGPLv3, see doc/lgpl.txt.
#ifndef COOLFluiD_FluxReconstructionMethod_PhysicalityGReKO_hh
#define COOLFluiD_FluxReconstructionMethod_PhysicalityGReKO_hh
#include "FluxReconstructionMethod/BasePhysicality.hh"
namespace COOLFluiD { namespace FluxReconstructionMethod {
/// Scale scalar polynomials toward their density-weighted cell means. Keeps
/// the quadrature integrals of rho*k, rho*gamma and rho*ReTheta unchanged.
/// Requires positive mean thermodynamics and positive quadrature weights;
/// it cannot repair an inadmissible scalar cell mean.
class PhysicalityGReKO : public BasePhysicality {
public:
 explicit PhysicalityGReKO(const std::string& name):BasePhysicality(name),m_gasConstant(0){}
 void setup();
protected:
 bool checkPhysicality();
 void enforcePhysicality();
 std::vector<CFreal> densityWeights();
 Common::SafePtr<RealVector> m_averageWeights;
 CFreal m_gasConstant;
};
}}
#endif
