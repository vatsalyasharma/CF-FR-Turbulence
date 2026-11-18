// Verify the actual PETSc reason; StdParSolveSys's message alone is insufficient.
#include "Petsc/PetscHeaders.hh"
#include "Petsc/StdParSolveSys.hh"
#include "Petsc/Petsc.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Common/BadValueException.hh"
#include <fstream>
#include <iomanip>
using namespace COOLFluiD;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Petsc;
class CheckedParSolveSys : public StdParSolveSys {
public:
 explicit CheckedParSolveSys(const std::string& n):StdParSolveSys(n){}
 void execute(){
  StdParSolveSys::execute();
  KSP ksp=getMethodData().getKSP();
  KSPConvergedReason reason;PetscInt iterations;PetscReal norm;
  const PetscErrorCode a=KSPGetConvergedReason(ksp,&reason);
  const PetscErrorCode b=KSPGetIterationNumber(ksp,&iterations);
  const PetscErrorCode c=KSPGetResidualNorm(ksp,&norm);
  if(a || b || c)throw Common::BadValueException(FromHere(),"PETSc convergence query failed");
  int rank;MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if(rank==0){
   const bool exists=std::ifstream("linear-convergence.csv").good();
   std::ofstream out("linear-convergence.csv",std::ios::app);out<<std::setprecision(17);
   if(!exists)out<<"iteration,ksp_iterations,reason,reported_residual_norm\n";
   out<<SubSystemStatusStack::getActive()->getNbIter()<<","<<iterations<<","<<int(reason)<<","<<norm<<"\n";
  }
  if(reason<=0)throw Common::BadValueException(FromHere(),"PETSc linear solve did not converge; update not accepted");
 }
};
MethodCommandProvider<CheckedParSolveSys,PetscLSSData,PetscModule> checked("CheckedParSolveSys");
