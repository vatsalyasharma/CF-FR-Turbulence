# Meshes

All meshes are 2D Gmsh files with quadrilateral elements of second-order geometry (Q2).
- **Boundary tags:** `inlet`, `upstream_slip` (the symmetry line ahead of the plate), `plate` (the no-slip wall from x = 0), `top` and `outlet`.
- **Lengths:** in metres. The leading edge is at x = 0.

| Mesh | Elements | Domain (x; height) | First wall element | Used by |
|---|---|---|---|---|
| `t3a_small` | 1,536 (48 streamwise × 32 normal: 8 ahead of the plate, 24 over the first metre, 16 to the outlet) | −0.25 … 20 m; 5 m | 1.0e-4 m high | γ–Re_θ P0, P1, P2 (TMR set-up) |
| `t3a_p5r` | 640 (40 × 16; the leading edge and the first metre as in `t3a_small`, coarser beyond) | −0.25 … 20 m; 5 m | 2.36e-4 m high | γ–Re_θ P5 (TMR set-up) |
| `t3a_small_exp` | 1,536, `t3a_small` scaled by 2e5/3.6e5 = 0.5556 | −0.139 … 11.1 m; 2.78 m | 5.6e-5 m high | γ-α P1, P2 at the experiment's conditions |
| `t3a_p5r_exp` | 640, `t3a_p5r` scaled by 0.5556 | −0.139 … 11.1 m; 2.78 m | 1.31e-4 m high | γ-α P5 at the experiment's conditions |
| `t3a_exp_jl3` | 12,288 (192 × 64), the L3 grid of the Jatayu finite-volume code, node for node | −0.139 … 11.1 m; 2.78 m | 5.65e-5 m high, 0.5 mm long at the leading edge | γ-α P1, code comparison with Jatayu |

- **TMR set-up:** the TMR meshes use NASA TMR's unit Reynolds number (Re = 2e5 per metre).
- **Experiment set-up:** the `_exp` meshes are the same geometry at the experiment's lengths (Re ≈ 3.6e5 per metre). Each case README names its mesh.
- **Restart cases:** these read the mesh from their `restart.CFmesh`, which already holds the geometry and the solution. Only the P0 case reads a Gmsh file directly.
