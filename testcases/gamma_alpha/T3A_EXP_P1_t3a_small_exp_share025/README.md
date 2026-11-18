# T3A_EXP_P1_t3a_small_exp_share025

| | |
|---|---|
| Model | SST γ-α, model version 2025 (`NavierStokes2DGammaAlpha`) |
| Conditions | ERCOFTAC T3A experiment's conditions: U 5.4 m/s, T 290.3 K, Re about 3.6e5 per metre, inlet k 0.0976 m²/s², ln ω 6.084, α 17.79 (inlet at x = −0.139 m; Tu 3.35 % at the leading edge) |
| Mesh | `meshes/t3a_small_exp` |
| Order | P1 |
| Ranks used for the reference | 2 |
| Pseudo-time stepping | implicit Newton, PETSc GMRES; CFL 30 ×1.02 → 600 (the reference was stopped after 1266 updates, once settled) |

## Start

`restart.CFmesh` = a settled P1 γ-α state with the same settings, except for the γ diffusion weighting. This case gives the same C_f to within 0.011 %.

## Run

From this folder:

```
export CF_BUILD=/path/to/build-fr
mpirun -np 2 $CF_BUILD/apps/Solver/coolfluid-solver --conf solver.xml \
    --bdir $PWD --scase $PWD/case.CFcase --ldir $CF_BUILD/dso
```

The solver writes into this folder:
- `final-iter_*.CFmesh`: the solution, every 100 updates;
- `convergence.plt*`: the residual history;
- `wall-flux-*.csv`: the wall output;
- `steady-residual-*.csv`: the per-equation residuals.

## Reference result

Transition (C_f = 3e-3) at Re_x 2.23e5 (experiment about 2.1e5); C_f minimum 2.01e-3 at Re_x 1.50e5; 16/16 ERCOFTAC stations within 10 %, RMS 5.5 %. **Uses LaminarFlowShare = 0.25**, a calibrated value that is not in the published model, plus a γ diffusion weighted by the wall-law γ.

- `reference_solution.CFmesh`: the final solution of the reference run.
- `reference_wall.dat`: its wall distribution, in Tecplot ASCII. The columns are x, Re_x, C_f, the wall values of k, γ and Re_θt (γ–Re_θ) or α (γ-α), and the quadrature weight.
