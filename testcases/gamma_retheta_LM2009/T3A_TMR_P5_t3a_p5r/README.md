# T3A_TMR_P5_t3a_p5r

| | |
|---|---|
| Model | γ–Re_θ (Langtry–Menter 2009) on SST k–ω, log-ω form (`NavierStokes2DGReKLogO`) |
| Conditions | NASA TMR T3A set-up: M 0.2, T 300 K, Re 2e5 per metre, inlet Tu and μt/μ as in TMR (k 24.80 m²/s², ln ω 8.700, inlet at x = −0.25 m) |
| Mesh | `meshes/t3a_p5r` |
| Order | P5 |
| Ranks used for the reference | 8 |
| Pseudo-time stepping | implicit Newton, PETSc GMRES; CFL 100 ×1.05 → 1000, at most 800 updates (the reference was stopped after 201 updates, once converged in practice) |

## Start

`restart.CFmesh` = a P5 state on `t3a_p5r` after 500 updates. That run started from the P2 solution interpolated onto this mesh.

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

**Matches NASA OVERFLOW.** Transition (C_f = 3e-3) at Re_x 1.70e5 (OVERFLOW 1.69e5, 0.2 % apart); C_f minimum 2.43e-3 at Re_x 1.16e5; 9/16 ERCOFTAC stations within 10 %. Converged in practice: C_f changes by at most 1.3e-4 per 100 updates.

- `reference_solution.CFmesh`: the final solution of the reference run.
- `reference_wall.dat`: its wall distribution, in Tecplot ASCII. The columns are x, Re_x, C_f, the wall values of k, γ and Re_θt (γ–Re_θ) or α (γ-α), and the quadrature weight.
