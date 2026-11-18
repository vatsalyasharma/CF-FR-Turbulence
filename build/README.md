# Building the solver

**Tested with:**
- GCC 13.3;
- Open MPI 4.1.6;
- PETSc 3.24.3;
- Boost 1.88;
- ParMETIS / METIS;
- CMake 3.28;
- Ubuntu 24.04.

## 1. Configure and build

```
export MPI_HOME=/path/to/openmpi           # contains bin/mpirun, include/, lib/libmpi.so
export DEPS=/path/to/deps                  # Boost, ParMETIS and METIS (include/, lib/); PETSc in $DEPS/petsc
mkdir build-fr && cd build-fr
CC=$MPI_HOME/bin/mpicc CXX=$MPI_HOME/bin/mpic++ \
cmake -C ../build/fr_turbulence.cmake \
      -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DBOOST_ROOT=$DEPS -DPARMETIS_HOME=$DEPS -DPETSC_DIR=$DEPS/petsc \
      ../source
make -j8
```

- The cache file turns on only what these cases need: the FR space method, Navier–Stokes, SST k–ω, γ–Re_θ (GReKO) and γ-α, the Newton iterator, PETSc, and the CFmesh/Gmsh readers and writers. Everything else is off.
- The libraries go to `build-fr/dso/` and the solver to `build-fr/apps/Solver/coolfluid-solver`.

## 2. Build the four extra plugins

Every case file loads these plugins from `source/extra_plugins/`:

| Plugin | Role |
|---|---|
| `WallFluxProbe` | wall output for γ–Re_θ |
| `WallFluxProbeGA` | wall output for γ-α |
| `SteadyResidualProbe` | per-equation residual files |
| `CheckedParSolveSys` | PETSc solve that stops on a failed linear solve |

They are single files built against the libraries from step 1:

```
B=$PWD                       # the build-fr folder
S=$B/../source
FLAGS="-DCF_HAVE_CONFIG_H -DCF_HAVE_MPI -I$DEPS/include -I$B -I$S/src -I$S/plugins -O2 -g -fPIC -std=c++14 -shared"
LINK="-L$B/dso -Wl,-rpath,$B/dso -lFramework -lEnvironment -lCommon -lMathTools"
$MPI_HOME/bin/mpic++ $FLAGS $S/extra_plugins/WallFluxProbe.cxx $LINK -lFluxReconstructionNavierStokes -lFluxReconstructionGReKO -lFluxReconstructionKOmega -lGReKO -o $B/dso/libWallFluxProbe.so
$MPI_HOME/bin/mpic++ $FLAGS $S/extra_plugins/SteadyResidualProbe.cxx $LINK -lFluxReconstructionNavierStokes -lFluxReconstructionGReKO -lFluxReconstructionKOmega -lGReKO -lNewtonMethod -o $B/dso/libSteadyResidualProbe.so
$MPI_HOME/bin/mpic++ $FLAGS -I$DEPS/petsc/include $S/extra_plugins/CheckedParSolveSys.cxx $LINK -L$DEPS/petsc/lib -Wl,-rpath,$DEPS/petsc/lib -lpetsc -lFluxReconstructionNavierStokes -lFluxReconstructionGReKO -lFluxReconstructionKOmega -lGReKO -lNewtonMethod -lPetscI -o $B/dso/libCheckedParSolveSys.so
$MPI_HOME/bin/mpic++ $FLAGS $S/extra_plugins/WallFluxProbeGA.cxx $LINK -lFluxReconstructionNavierStokes -lFluxReconstructionGammaAlpha -lFluxReconstructionKOmega -lFluxReconstructionMethod -lGammaAlpha -o $B/dso/libWallFluxProbeGA.so
```

Add the Boost define that your configure step uses (for Boost 1.88, `-DCF_HAVE_BOOST_1_88`) to FLAGS. It appears in `build-fr/plugins/FluxReconstructionTurb/CMakeFiles/*/flags.make`.
