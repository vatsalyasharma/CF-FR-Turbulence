// Copyright (C) 2026 COOLFluiD contributors; LGPLv3, see doc/lgpl.txt.
#include <algorithm>
#include <cmath>
#include <limits>
#include "Common/BadValueException.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/GeometricEntity.hh"
#include "NavierStokes/EulerTerm.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "FluxReconstructionTurb/FluxReconstructionKOmega.hh"
#include "FluxReconstructionTurb/PhysicalityGReKO.hh"
using namespace COOLFluiD::Framework;
namespace COOLFluiD { namespace FluxReconstructionMethod {
MethodCommandProvider<PhysicalityGReKO,FluxReconstructionSolverData,FluxReconstructionKOmegaModule> grScalarPositivity("PhysicalityGReKO");

void PhysicalityGReKO::setup() {
 BasePhysicality::setup();
 const std::string model=PhysicalModelStack::getActive()->getConvectiveName();
 const std::string vars=getMethodData().getUpdateVarStr();
 if(!((m_dim==2 && model=="Euler2DGReKLogO" && vars=="Puvt") ||
      (m_dim==3 && model=="Euler3DGReKLogO" && vars=="Pvt")) || m_nbrEqs!=m_dim+6)
   throw Common::BadValueException(FromHere(),"PhysicalityGReKO requires primitive 2D/3D gamma-ReTheta log-omega variables");
 auto& elements=getMethodData().getFRLocalData();
 if(elements.size()!=1)throw Common::BadValueException(FromHere(),"PhysicalityGReKO requires a single FR element type");
 m_averageWeights=elements[0]->getCellAvgSolCoefs();
 m_solPntsLocalCoords=elements[0]->getSolPntsLocalCoords();
 m_gasConstant=PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm().d_castTo<Physics::NavierStokes::EulerTerm>()->getR();
 for(CFuint i=0;i<m_nbrSolPnts;++i)
   if(!((*m_averageWeights)[i]>0 && std::isfinite((*m_averageWeights)[i])))
     throw Common::BadValueException(FromHere(),"PhysicalityGReKO requires positive solution-point quadrature weights");
}

bool PhysicalityGReKO::checkPhysicality() {
 bool valid=true;
 const CFuint k=m_dim+2;
 for(CFuint i=0;i<m_nbrSolPnts;++i){
   const State& q=*(*m_cellStates)[i];
   for(CFuint e=0;e<m_nbrEqs;++e)if(!std::isfinite(q[e]))
     throw Common::BadValueException(FromHere(),"PhysicalityGReKO cannot limit a nonfinite state");
   if(!(q[0]>0 && q[m_dim+1]>0))throw Common::BadValueException(FromHere(),"PhysicalityGReKO requires positive nodal pressure and temperature");
   valid=valid && q[k]>=0 && q[k+2]>=0 && q[k+3]>0;
 }
 for(CFuint i=0;i<m_maxNbrFlxPnts;++i){
   const RealVector& q=m_cellStatesFlxPnt[i];
   for(CFuint e=0;e<m_nbrEqs;++e)if(!std::isfinite(q[e]))
     throw Common::BadValueException(FromHere(),"PhysicalityGReKO cannot limit a nonfinite face state");
   if(!(q[0]>0 && q[m_dim+1]>0))throw Common::BadValueException(FromHere(),"PhysicalityGReKO requires positive face pressure and temperature");
   valid=valid && q[k]>=0 && q[k+2]>=0 && q[k+3]>0;
 }
 return valid;
}

std::vector<CFreal> PhysicalityGReKO::densityWeights() {
 const auto jac=m_cell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);
 std::vector<CFreal> weights(m_nbrSolPnts);
 for(CFuint i=0;i<m_nbrSolPnts;++i){
   const State& q=*(*m_cellStates)[i];
   weights[i]=(*m_averageWeights)[i]*jac[i]*q[0]/(m_gasConstant*q[m_dim+1]);
   if(!(weights[i]>0 && std::isfinite(weights[i])))
     throw Common::BadValueException(FromHere(),"PhysicalityGReKO requires positive finite density-weighted cell measures");
 }
 return weights;
}

void PhysicalityGReKO::enforcePhysicality() {
 const auto weights=densityWeights();CFreal mass=0;
 for(CFreal w:weights)mass+=w;
 if(!(mass>0 && std::isfinite(mass)))throw Common::BadValueException(FromHere(),"PhysicalityGReKO requires finite positive cell mass");
 const CFuint ids[3]={m_dim+2,m_dim+4,m_dim+5};
 CFreal average[3],theta[3];
 // Preflight all three means before changing any state.
 for(unsigned scalar=0;scalar<3;++scalar){
   const CFuint e=ids[scalar];CFreal integral=0,low=std::numeric_limits<CFreal>::max();
   for(CFuint i=0;i<m_nbrSolPnts;++i){const CFreal q=(*(*m_cellStates)[i])[e];integral+=weights[i]*q;low=std::min(low,q);}
   for(CFuint i=0;i<m_maxNbrFlxPnts;++i)low=std::min(low,m_cellStatesFlxPnt[i][e]);
   average[scalar]=integral/mass;theta[scalar]=1;
   if(!std::isfinite(average[scalar]) || average[scalar]<0 || (scalar==2 && average[scalar]==0))
     throw Common::BadValueException(FromHere(),"PhysicalityGReKO cannot preserve an inadmissible scalar cell mean; reduce the time step");
   if(low<0 || (scalar==2 && low==0)){
     if(average[scalar]==0)theta[scalar]=0;
     else theta[scalar]=(1-1.e-12)*average[scalar]/(average[scalar]-low);
   }
 }
 for(unsigned scalar=0;scalar<3;++scalar)if(theta[scalar]<1){
   const CFuint e=ids[scalar];
   for(CFuint i=0;i<m_nbrSolPnts;++i){
     State& q=*(*m_cellStates)[i];q[e]=(1-theta[scalar])*average[scalar]+theta[scalar]*q[e];
   }
 }
 computeFlxPntStates(m_cellStatesFlxPnt);
 if(!checkPhysicality())throw Common::BadValueException(FromHere(),"PhysicalityGReKO failed to restore scalar positivity");
}
}}
