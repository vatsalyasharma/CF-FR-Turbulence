# Test cases: ERCOFTAC T3A flat plate

Each case has its own folder with:

| File | Content |
|---|---|
| `case.CFcase` | the full solver input |
| `solver.xml` | an empty solver configuration that the command line needs |
| `restart.CFmesh` | the start solution (mesh and solution together); the P0 case links the Gmsh mesh instead |
| `reference_solution.CFmesh`, `reference_wall.dat` | the reference result |
| `README.md` | the case, its mesh, how to run it, and the expected result |

## γ–Re_θ (Langtry–Menter 2009), NASA TMR set-up: the main result

A P0 → P1 → P2 → P5 sequence. Each case starts from the previous one.

| Case | Mesh | Order | Transition, C_f = 3e-3 at Re_x |
|---|---|---|---|
| `gamma_retheta_LM2009/T3A_TMR_P0_t3a_small` | `t3a_small` | P0 | none (first order, turbulent) |
| `gamma_retheta_LM2009/T3A_TMR_P1_t3a_small` | `t3a_small` | P1 | 1.22e5 |
| `gamma_retheta_LM2009/T3A_TMR_P2_t3a_small` | `t3a_small` | P2 | **1.65e5** |
| `gamma_retheta_LM2009/T3A_TMR_P5_t3a_p5r` | `t3a_p5r` | P5 | **1.70e5** |
| NASA OVERFLOW (TMR), same model | 2817 × 1537 points | — | 1.69e5 |

No tuning and no ad hoc change to the published model.

## SST γ-α, experiment's conditions

| Case | Mesh | Order | LaminarFlowShare | Transition at Re_x | ERCOFTAC within 10 % |
|---|---|---|---|---|---|
| `gamma_alpha/T3A_EXP_P1_t3a_small_exp_share025` | `t3a_small_exp` | P1 | 0.25 | 2.23e5 | 16/16 |
| `gamma_alpha/T3A_EXP_P2_t3a_small_exp_share025` | `t3a_small_exp` | P2 | 0.25 | 2.16e5 | 15/16 |
| `gamma_alpha/T3A_EXP_P5_t3a_p5r_exp_share025` | `t3a_p5r_exp` | P5 | 0.25 | 2.17e5 | 15/16 |
| `gamma_alpha/T3A_EXP_P1_jatayu_L3_share1` | `t3a_exp_jl3` | P1 | 1 (as published) | 2.75e5 | 6/16 |

- The experiment's transition is about 2.1e5.
- LaminarFlowShare = 0.25 is a calibrated value, not part of the published model. The share-1 case shows the model as published.

## Notes

- **Build of the reference runs:**
  - The γ–Re_θ P5 case and the γ-α P2 and Jatayu-grid cases were run with exactly this source.
  - The other references came from earlier builds of the same model code.
  - Re-running two of them on this source reproduced them: the γ–Re_θ P5 to within 0.02 %, and the γ-α P1 as well.
- **Ranks:** any number of MPI ranks works. The reference rank counts are listed in each README.
