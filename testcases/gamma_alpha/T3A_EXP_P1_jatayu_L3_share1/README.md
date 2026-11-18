# T3A_EXP_P1_jatayu_L3_share1

| | |
|---|---|
| Model | SST γ-α, model version 2025 (`NavierStokes2DGammaAlpha`) |
| Conditions | ERCOFTAC T3A experiment's conditions: U 5.4 m/s, T 290.3 K, Re about 3.6e5 per metre, inlet k 0.0976 m²/s², ln ω 6.084, α 17.79 (inlet at x = −0.139 m; Tu 3.35 % at the leading edge) |
| Mesh | `meshes/t3a_exp_jl3` |
| Order | P1 |
| Ranks used for the reference | 8 |
| Pseudo-time stepping | implicit Newton, PETSc GMRES; CFL 10 ×1.02 → 600, 3000 updates |

## Start

`restart.CFmesh` = a P1 γ-α solution on another mesh, interpolated onto this one.

## Run

From this folder:

```
export CF_BUILD=/path/to/build-fr
mpirun -np 8 $CF_BUILD/apps/Solver/coolfluid-solver --conf solver.xml \
    --bdir $PWD --scase $PWD/case.CFcase --ldir $CF_BUILD/dso
```

The solver writes into this folder:
- `final-iter_*.CFmesh`: the solution, every 100 updates;
- `convergence.plt*`: the residual history;
- `wall-flux-*.csv`: the wall output;
- `steady-residual-*.csv`: the per-equation residuals.

## Reference result

**The γ-α model as published (LaminarFlowShare = 1, inlet γ 0.99)** on the L3 grid of the Jatayu finite-volume code. Transition is late: Re_x 2.75e5 against about 2.1e5 measured; 6/16 within 10 %. Kept for comparison with Jatayu, which places transition earlier on the same grid with the same settings. Not fully converged: C_f changes 0.12 % per 100 updates.

- `reference_solution.CFmesh`: the final solution of the reference run.
- `reference_wall.dat`: its wall distribution, in Tecplot ASCII. The columns are x, Re_x, C_f, the wall values of k, γ and Re_θt (γ–Re_θ) or α (γ-α), and the quadrature weight.
