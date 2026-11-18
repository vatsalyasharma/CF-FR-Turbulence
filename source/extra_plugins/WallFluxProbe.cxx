// Read-only measurement of the actual compact BR2 wall flux, before Jacobian perturbations.
#include <fstream>
#include <iomanip>
#include <set>
#include <mpi.h>
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "FluxReconstructionTurb/FluxReconstructionKOmega.hh"
#include "FluxReconstructionTurb/DiffBndCorrectionsRHSJacobFluxReconstructionTurb.hh"
#include "FluxReconstructionTurb/DiffBndCorrectionsRHSFluxReconstructionTurb.hh"
#include "FluxReconstructionNavierStokes/NSJacobBndGradientComputer.hh"
#include "FluxReconstructionNavierStokes/NSBndGradientComputer.hh"
using namespace COOLFluiD;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::FluxReconstructionMethod;
template<class Base> class WallFluxProbe : public Base {
 CFuint lastIter=~CFuint(0);std::set<CFuint> seen;
public:
 explicit WallFluxProbe(const std::string& n):Base(n){}
 void computeInterfaceFlxCorrection() {
  Base::computeInterfaceFlxCorrection();
  if(this->getCurrentTRS()->getName()!="plate" || !(*this->m_cellStates)[0]->isParUpdatable())return;
  const CFuint iteration=SubSystemStatusStack::getActive()->getNbIter();
  if(iteration!=lastIter){lastIter=iteration;seen.clear();}
  if(!seen.insert(this->m_face->getID()).second)return;
  int rank;MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  const std::string file="wall-flux-P"+std::to_string(rank)+".csv";
  const bool exists=std::ifstream(file).good();
  std::ofstream out(file,std::ios::app);out<<std::setprecision(17);
  if(!exists)out<<"iteration,face,point,x,y,nx,ny,weight,flux_u,flux_energy,flux_gamma,flux_retheta,u_wall,v_wall,k_wall,gamma_wall,retheta_wall\n";
  for(CFuint i=0;i<this->m_nbrFaceFlxPnts;++i) {
   const auto& x=this->m_flxPntCoords[i];const auto& n=this->m_unitNormalFlxPnts[i];
   const auto& f=this->m_flxPntRiemannFlux[i];const auto& q=*this->m_bndStateFlux[i];
   out<<iteration<<","<<this->m_face->getID()<<","<<i<<","<<x[0]<<","<<x[1]<<","<<n[0]<<","<<n[1]<<","<<this->m_faceJacobVecAbsSizeFlxPnts[i]*(*this->m_faceIntegrationCoefs)[i]<<","<<f[1]<<","<<f[3]<<","<<f[6]<<","<<f[7]<<","<<q[1]<<","<<q[2]<<","<<q[4]<<","<<q[6]<<","<<q[7]<<"\n";
  }
 }
};
MethodCommandProvider<WallFluxProbe<DiffBndCorrectionsRHSJacobFluxReconstructionTurb>,FluxReconstructionSolverData,FluxReconstructionKOmegaModule> wj("DiffBndCorrectionsRHSJacobGRWall");
MethodCommandProvider<NSJacobBndGradientComputer,FluxReconstructionSolverData,FluxReconstructionKOmegaModule> cj("ConvBndCorrectionsRHSJacobGRWall");
MethodCommandProvider<WallFluxProbe<DiffBndCorrectionsRHSFluxReconstructionTurb>,FluxReconstructionSolverData,FluxReconstructionKOmegaModule> we("DiffBndCorrectionsRHSGRWall");
MethodCommandProvider<NSBndGradientComputer,FluxReconstructionSolverData,FluxReconstructionKOmegaModule> ce("ConvBndCorrectionsRHSGRWall");
