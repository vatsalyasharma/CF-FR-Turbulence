# CF-FR-Turbulence

COOLFluiD flux reconstruction (FR) solver with RANS turbulence and laminar–turbulent transition models, plus validated test cases.

- **Turbulence model:** SST k–ω.
- **Transition model γ–Re_θ:** Langtry–Menter 2009, the same model NASA's OVERFLOW uses.
- **Transition model γ-α:** SST γ-α.

## Main result

**T3A flat plate (ERCOFTAC), γ–Re_θ at the NASA TMR set-up.** At P2 and P5 it reproduces NASA OVERFLOW with no tuning or ad hoc model change:

| Order | Mesh | C_f = 3e-3 crossing (Re_x) | C_f minimum (Re_x) | Case |
|---|---|---|---|---|
| P0 | t3a_small | none (turbulent from the leading edge; first order) | — | `testcases/gamma_retheta_LM2009/T3A_TMR_P0_t3a_small` |
| P1 | t3a_small | 1.22e5 (early) | 2.89e-3 (8.6e4) | `testcases/gamma_retheta_LM2009/T3A_TMR_P1_t3a_small` |
| P2 | t3a_small | **1.65e5** | 2.45e-3 (1.14e5) | `testcases/gamma_retheta_LM2009/T3A_TMR_P2_t3a_small` |
| P5 | t3a_p5r | **1.70e5** | 2.43e-3 (1.16e5) | `testcases/gamma_retheta_LM2009/T3A_TMR_P5_t3a_p5r` |
| — | OVERFLOW, LM2009, 2817 × 1537 points | 1.69e5 | 2.40e-3 (1.12e5) | NASA TMR |

- The transition location agrees with OVERFLOW within 2.4 % at P2 and 0.2 % at P5.
- From the start of the rise on, C_f at the ERCOFTAC stations agrees with OVERFLOW within 0.6–3.8 points.
- Both codes transition slightly ahead of the experiment (about 2.1e5).

The γ-α cases are in `testcases/gamma_alpha`. Their best agreement with the experiment uses LaminarFlowShare = 0.25, a calibrated value that is not in the published model. A case with the model as published (share 1) is included for comparison. See `testcases/README.md`.

## Layout

| Folder | Content |
|---|---|
| `source/` | The complete COOLFluiD source used for every result here: `src/`, `plugins/`, `apps/`, `cmake/`, `config/`, `tools/`, `doc/`, `CMakeLists.txt`. `source/extra_plugins/` holds four small plugins that the case files load. |
| `build/` | `fr_turbulence.cmake`, the initial CMake cache that selects the plugins and options of the build, and `README.md` with the build steps. |
| `meshes/` | The Gmsh meshes (second-order, Q2 geometry), one folder per mesh. See `meshes/README.md`. |
| `testcases/` | One folder per case, each with its own case file, start solution, reference solution and README. See `testcases/README.md`. |

No binaries, build trees or third-party dependencies are included.

## Quick start

1. Build the solver: see `build/README.md`.
2. Run a case:
   ```
   export CF_BUILD=/path/to/build
   cd testcases/gamma_retheta_LM2009/T3A_TMR_P2_t3a_small
   mpirun -np 4 $CF_BUILD/apps/Solver/coolfluid-solver --conf solver.xml \
       --bdir $PWD --scase $PWD/case.CFcase --ldir $CF_BUILD/dso
   ```
3. Compare the wall output with `reference_wall.dat` in the case folder.
