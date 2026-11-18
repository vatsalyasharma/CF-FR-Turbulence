// Read-only spatial RHS measurement before pseudo-time assembly and linear solve.
#include <fstream>
#include <iomanip>
#include <cmath>
#include <mpi.h>
#include "NewtonMethod/NewtonIteratorData.hh"
#include "NewtonMethod/NewtonMethod.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/State.hh"
#include "Framework/SubSystemStatus.hh"
using namespace COOLFluiD;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Numerics::NewtonMethod;
class SteadyResidualProbe : public NewtonIteratorCom {
 DataSocketSink<CFreal> rhs;
 DataSocketSink<State*, GLOBAL> states;
public:
 explicit SteadyResidualProbe(const std::string& n):NewtonIteratorCom(n),rhs("rhs"),states("states"){}
 std::vector<SafePtr<BaseDataSocketSink> > needsSockets(){return {&rhs,&states};}
 void execute(){
  const auto r=rhs.getDataHandle(); const auto s=states.getDataHandle();
  const CFuint neq=PhysicalModelStack::getActive()->getNbEq();
  int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  const std::string file="steady-residual-P"+std::to_string(rank)+".csv";
  const bool exists=std::ifstream(file).good();
  std::ofstream out(file,std::ios::app);out<<std::setprecision(17);
  if(!exists)out<<"iteration,region,equation,count,sum_squares,max_abs\n";
  for(CFuint region=0;region<2;++region)for(CFuint e=0;e<neq;++e){
   CFreal sum=0.,peak=0.;CFuint count=0;
   for(CFuint i=0;i<s.size();++i)if(s[i]->isParUpdatable()){
    const auto& x=s[i]->getCoordinates();
    if(region && !(x[0]>=0. && x[0]<=1. && x[1]<=0.1))continue;
    const CFreal v=r[i*neq+e];sum+=v*v;peak=std::max(peak,std::abs(v));++count;
   }
   out<<SubSystemStatusStack::getActive()->getNbIter()<<","<<region<<","<<e<<","<<count<<","<<sum<<","<<peak<<"\n";
  }
 }
};
MethodCommandProvider<SteadyResidualProbe,NewtonIteratorData,NewtonMethodModule> sr("SteadyResidualProbe");
