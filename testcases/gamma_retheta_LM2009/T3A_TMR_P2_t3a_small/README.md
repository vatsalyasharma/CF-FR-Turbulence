# T3A_TMR_P2_t3a_small

| | |
|---|---|
| Model | γ–Re_θ (Langtry–Menter 2009) on SST k–ω, log-ω form (`NavierStokes2DGReKLogO`) |
| Conditions | NASA TMR T3A set-up: M 0.2, T 300 K, Re 2e5 per metre, inlet Tu and μt/μ as in TMR (k 24.80 m²/s², ln ω 8.700, inlet at x = −0.25 m) |
| Mesh | `meshes/t3a_small` |
| Order | P2 |
| Ranks used for the reference | 4 |
| Pseudo-time stepping | implicit Newton, PETSc GMRES; CFL 1 ×1.02 → 1000, 3500 updates |

## Start

`restart.CFmesh` = the converged P1 solution (an equivalent P1 run from another start; it is the same solution as the P1 case to plotting accuracy). The case raises it to P2.

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

**Matches NASA OVERFLOW.** Transition (C_f = 3e-3) at Re_x 1.65e5 (OVERFLOW 1.69e5); C_f minimum 2.45e-3 at Re_x 1.14e5 (OVERFLOW 2.40e-3 at 1.12e5); 9/16 ERCOFTAC stations within 10 % (OVERFLOW 12/16; the gap is the first three laminar stations, which lie just either side of the 10 % line). Steady: from update 1500 on the residuals sit on a plateau at 4e-7 to 5.5e-6 of their peak, and C_f changes by at most 2.8e-5 per 100 updates.

- `reference_solution.CFmesh`: the final solution of the reference run.
- `reference_wall.dat`: its wall distribution, in Tecplot ASCII. The columns are x, Re_x, C_f, the wall values of k, γ and Re_θt (γ–Re_θ) or α (γ-α), and the quadrature weight.
