# T3A_TMR_P0_t3a_small

| | |
|---|---|
| Model | γ–Re_θ (Langtry–Menter 2009) on SST k–ω, log-ω form (`NavierStokes2DGReKLogO`) |
| Conditions | NASA TMR T3A set-up: M 0.2, T 300 K, Re 2e5 per metre, inlet Tu and μt/μ as in TMR (k 24.80 m²/s², ln ω 8.700, inlet at x = −0.25 m) |
| Mesh | `meshes/t3a_small` |
| Order | P0 |
| Ranks used for the reference | 4 |
| Pseudo-time stepping | implicit Newton, PETSc GMRES; CFL 1 ×1.02 → 1000, at most 3000 updates |

## Start

No start file. The case reads the Gmsh mesh and starts from the uniform free stream. Link the mesh first:
```
ln -s ../../../meshes/t3a_small/t3a_small_q2.msh input.msh
```

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

Converged after 849 updates. First order: turbulent from the leading edge, with no transition (C_f at Re_x 5e5 = 5.41e-3). It is the start of the P1 case.

- `reference_solution.CFmesh`: the final solution of the reference run.
- `reference_wall.dat`: its wall distribution, in Tecplot ASCII. The columns are x, Re_x, C_f, the wall values of k, γ and Re_θt (γ–Re_θ) or α (γ-α), and the quadrature weight.
