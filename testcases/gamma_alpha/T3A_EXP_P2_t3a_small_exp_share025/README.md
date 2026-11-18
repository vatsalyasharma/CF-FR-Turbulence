# T3A_EXP_P2_t3a_small_exp_share025

| | |
|---|---|
| Model | SST γ-α, model version 2025 (`NavierStokes2DGammaAlpha`) |
| Conditions | ERCOFTAC T3A experiment's conditions: U 5.4 m/s, T 290.3 K, Re about 3.6e5 per metre, inlet k 0.0976 m²/s², ln ω 6.084, α 17.79 (inlet at x = −0.139 m; Tu 3.35 % at the leading edge) |
| Mesh | `meshes/t3a_small_exp` |
| Order | P2 |
| Ranks used for the reference | 4 |
| Pseudo-time stepping | implicit Newton, PETSc GMRES; CFL 10 ×1.02 → 600, 3000 updates |

## Start

`restart.CFmesh` = the P1 reference solution of `T3A_EXP_P1_t3a_small_exp_share025`. The case raises it to P2.

## Run

From this folder:

```
export CF_BUILD=/path/to/build-fr
mpirun -np 4 $CF_BUILD/apps/Solver/coolfluid-solver --conf solver.xml \
    --bdir $PWD --scase $PWD/case.CFcase --ldir $CF_BUILD/dso
```

The solver writes into this folder:
- `final-iter_*.CFmesh`: the solution, every 100 updates;
- `convergence.plt*`: the residual history;
- `wall-flux-*.csv`: the wall output;
- `steady-residual-*.csv`: the per-equation residuals.

## Reference result

Transition at Re_x 2.16e5; C_f minimum 1.94e-3 at Re_x 1.57e5; 15/16 within 10 %, RMS 7.3 %. Resolved along the plate (C_f jumps ≤ 0.9 % between elements). Not fully converged at 3000 updates: C_f still changes 0.13 % per 100 updates near the start of the rise, falling by about 10 % per 100 updates. Same settings as the P1 case (share 0.25).

- `reference_solution.CFmesh`: the final solution of the reference run.
- `reference_wall.dat`: its wall distribution, in Tecplot ASCII. The columns are x, Re_x, C_f, the wall values of k, γ and Re_θt (γ–Re_θ) or α (γ-α), and the quadrature weight.
