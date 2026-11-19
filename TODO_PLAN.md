# RANS-LTT FR: action plan (FWO project + Andrea Lani's case list)

2026-11-19. This is the short, decision-oriented list. The evidence and detail are in `STATUS.md`, `TODO.md` and `MANIFEST.md` in this folder.

**Focus: hypersonic cones with order blending (OB) and/or matrix-free (JFNK).** Ray did neither, so any converged OB or JFNK result compared with his solutions is reportable progress.

---

## 1. What is available now

### 1.1 Code (current build: `~/Desktop/COOLFluiD/OPENMPI/optim`, branch `fr-turbulence`)

| Item | State | Evidence |
|---|---|---|
| γ-α (SST) and LM γ-Reθ source terms, 2D and 3D, analytical source Jacobian | ✅ in code | `plugins/FluxReconstructionTurb/{GammaAlpha,NavierStokesGReKO}{2D,3D}SourceTerm*` |
| **Order blending for turbulence, implicit** | ✅ in code (RD, 2026-09-21, committed), **never used in a turbulence case** | `ConvRHSJacobTurbBlending`, `ConvBndCorrectionsRHSJacobTurbBlending`, `ConvBndCorrectionsRHSJacobGammaAlphaBlending`, `GammaAlphaJacobBndGradientComputerBlending`. Selected by `UseBlending = true` + `PreProcessCom = OrderBlending` (the solver appends "Blending" to `ConvSolveCom` and to the boundary command) |
| Order blending for turbulence, explicit | ⚠️ only generic fallbacks | No `ConvBndCorrectionsRHSGammaAlphaBlending` exists; the explicit path falls back to the generic `ConvBndCorrectionsRHSBlending`. Not tested with γ-α |
| **Matrix-free (JFNK)** | ✅ infrastructure, **never used with turbulence** | PETSc `ParJFSetup` / `ParJFSolveSys`, `PCSHELL` + `FRP0Precond` (p-multigrid), `DoComputeJacobian = false` |
| **Templates** | ✅ | JFNK+OB: `plugins/NavierStokes/testcases/Wedge/wedge2DFR-JFNK-OB-PhysCFL.CFcase`. OB, hypersonic: `NavierStokes/testcases/Cylinder_FR/gnoffoFR_P{1,3}_OB.CFcase`, `NEQ/testcases/TCNEQ/{HEG,Hornung}_cyl_FR/*`. JFNK: `NavierStokes/testcases/Cylinder/gnoffoFR2DP2-JFNK.CFcase` |
| Wall probes, residual probe, checked linear solve (from CF-FR-Turbulence) | ✅ built locally and working | used by the T3A runs (§3) |
| Ray's code (`~/Desktop/CF_RAY_Version/COOLFluiD_RAY`) | source only, not built | needed only for a strict reproduction of Ray (optional) |

### 1.2 Hypersonic cases and Ray's reference solutions

Folders are in `02_validation_baseline/` unless stated.

| Andrea item | What exists | Runnable now? |
|---|---|---|
| **Kimmel 2D, P1–P3** (M 7.93, axisymmetric, 1764 quads) | Ray's original CFcases `RAY_kimmel2D_GA_P1/P2/P3_expl`. Each reads **Ray's own Pn solution** (`kimmelP1/P2New/P3.CFmesh`), so it doubles as the reference. `RAYPROD_kimmel2D_GA_P2-CTR_expl` is a production variant (Tw 320 K). Ray's plots: `05_reference_results_RAY/PLT_cone/kimmelP1..P5.plt` | ✅ but **MaxIter = 1 in P1/P3** (set it). The `.inter` files are placeholders; copy settings from `RAYPROD_*/transition*.inter` |
| Kimmel 2D P0 | `RAY_kimmel2D_GA_P0_impl` (Newton, starts from mesh) | ✅ |
| **DiCristina 2D, P1–P2** (M 10.2) | `RAY_dicristina2D_GA_P1_expl` (Ray's P1 solution `dicristinaP1Tu075`), `RAY_dicristina2D_GA_P2_expl` (Ray's P2 `…P2Tu075ns8Cmg2-1`). Tu and spot-growth sweeps: `05_reference_results_RAY/Res_cone/Dicristina/` | ✅ (P2 MaxIter = 1) |
| **Kimmel 3D, P1–P3** (Q2 hexa) | **A complete Ray ladder on the 3240-hex tip sector:**<ul><li>P1 solution `05…/Res_cone/3D/kimmel3DP1.CFmesh`</li><li>`RAY_kimmel3D_GA_P2_expl` (upgrades P1 → P2)</li><li>`RAYPROD_kimmel3D_GA_P3bis_expl` (P2 `kimmel3DP2` + P3 restart `run2`, production `.inter`)</li></ul>Alternative: the 19,440-hex half cone at 1° AoA (`RAYPROD_kimmel3D_GA_P1_expl`, P1 only) | ✅ P2 and P3bis. P1 from Ray's solution |
| **DiCristina 3D, P1–P2** | `RAY_dicristina3D_GA_P1_impl` (P1 restart, 59,904 hexa, Newton). `RAY_dicristina3D_GA_P2_expl` (upgrades Ray's P1 `dicristinaNewP1Run4`, 79,248 hexa). Extra P1 solutions: `05…/Res_cone/dicriP1BeforeCrash*.CFmesh` (59,904) | ✅ (**two different meshes:** choose one) |
| HIFiRE (M 2.64, P1–P4), HEXAFLY (M 7.03) | `RAYPROD_hifire3D_*` (restarts missing; P4 solution available), `RAY_hexafly3D_GA_P1_expl` (setup only) | ⚠️ later |

### 1.3 Low-speed reference (build and model check, not hypersonic)

`~/Desktop/Vincent-scripts/CF-cases/Validation_TURB_26`: T3A flat plate, LM 2009 and γ-α, with references (see §3).

---

## 2. What is missing (blocking or important)

| # | Missing | Impact | Fix |
|---|---|---|---|
| M1 | **Any turbulence case with OB or JFNK** | Andrea's list cannot start without one | Build from the templates in §1.1 (first tasks below) |
| M2 | Experimental reference data for **Kimmel** and **DiCristina** (wall heat flux / Stanton versus x), and HIFiRE flight data | Can only compare with Ray, not with experiment | Get them from Ray's thesis or papers (Kimmel 1997; DiCristina 1970) and store as `.dat` next to the cases |
| M3 | Ray's interactive files for the `RAY_*` 2D cases | Runtime LLAV and laminar-wall settings unknown | Copy from `RAYPROD_*` (same keys); record the values used |
| M4 | Lost production restarts: Kimmel 3D P2/P3 on 19,440 hexa, DiCristina 3D P1/P3, HIFiRE P1–P3, Kimmel 2D P0 (iter 70000) | Exact continuation of Ray's production runs is impossible | Use the 3240-hex Kimmel ladder and the DiCristina P1 solutions that do exist; restart others from their upgrade inputs |
| M5 | Explicit OB for γ-α boundary terms | "Explicit + OB" may be inconsistent at the wall | Test the generic fallback first; implement `ConvBndCorrectionsRHSGammaAlphaBlending` if needed |
| M6 | JFNK preconditioner with the turbulence/transition source Jacobian | JFNK may stall on stiff k-ω-γ-α sources | Check what `FRP0Precond` includes; add the analytical source Jacobian to the P0 / block-Jacobi preconditioner |
| M7 | Uncommitted turbulence work on `fr-turbulence` (2025 γ-α, Karsch, probes…) | Risk of loss; results not reproducible from git | Commit or merge (merge kit in `~/Desktop/MERGE_CF-FR-Turbulence`) |

---

## 3. Validation_TURB_26 (T3A flat plate): what is already fulfilled

Checked 2026-11-19 16:35 against each case's `reference_wall.dat`:
- **C_f:** shear from the wall probe, normalised on the laminar stretch.
- **Wall coverage:** both rank files merged (`wall-flux-P1.csv` and `-P3.csv` are ranks 1 and 3, not polynomial orders).

| Case | Status | Result |
|---|---|---|
| LM 2009 T3A **P0** (`T3A_TMR_FR/t3aFR_P0`, local chain from mesh, current build) | ✅ **reproduced** | 829 Newton updates to residual about −8; C_f within **0.8 %** of the reference at all 38 wall points |
| LM 2009 T3A **P1** (`t3aFR_P1`) | ✅ **reproduced** | 1499 updates, residual about −8. C_f within **1.2 %** at 79 points; C_f minimum 2.888e-3 at Re_x 8.6e4 (ref. 2.885e-3); C_f = 3e-3 crossing at Re_x 1.18e5 (ref. 1.22e5, on the steep rise) |
| LM 2009 T3A **P2** (`t3aFR_P2`) | ⏳ **running** | started 16:18; iteration 107, CFL 8, about 8 s/iteration at 16:33. Reference: transition 1.65e5 |
| LM 2009 T3A **P5** (`T3A_TMR_P5_t3a_p5r`) | ❌ not run | reference: 1.70e5 (OVERFLOW 1.69e5) |
| γ-α T3A **P1 / P2 / P5**, LaminarFlowShare 0.25 (`gamma_alpha/T3A_EXP_*_share025`) | ❌ not run | references: 2.23e5 / 2.16e5 / 2.17e5 (experiment about 2.1e5) |
| γ-α T3A P1 Jatayu grid, share 1 (`T3A_EXP_P1_jatayu_L3_share1`) | ❌ not run | reference 2.75e5 |

**Meaning:**
- The current local build reproduces the LM 2009 references at P0 and P1, so the LM FR path and the new probes are verified.
- **The γ-α path, which the hypersonic cases use, is not yet verified on this build.** Run `T3A_EXP_P1_t3a_small_exp_share025` early.

---

## 4. What to run first (in this order)

Core counts follow the VSC rule: 100–200 solution points per core minimum, 1000 maximum.

| Case | Solution points | Cores |
|---|---|---|
| Kimmel 2D P2 | about 16k | 16–100 |
| Kimmel 2D P3 | about 28k | 28–200 |
| Kimmel 3D P3 (3240 hexa) | about 207k | 210–1000 |
| DiCristina 3D P2 (79,248 hexa) | about 2.1M | 2100+ (Ray used 288) |

### Step 0: finish the build check (today and tomorrow, workstation)
1. Let LM T3A P2 finish; compare it with its reference (1.65e5).
2. Run **γ-α T3A P1 share 0.25** (`Validation_TURB_26/testcases/gamma_alpha/T3A_EXP_P1_t3a_small_exp_share025`). It verifies the γ-α path of the current build. *Gate: transition within a few % of 2.23e5.*

### Step 1: Kimmel 2D control, no OB (baseline with the current build)
3. `RAY_kimmel2D_GA_P1_expl`: restart from Ray's `kimmelP1.CFmesh`.
   - Set MaxIter; copy the production `.inter` values.
   - Run with Ray's settings (LLAV, explicit).
   - *Gate: wall heat flux matches `kimmelP1.plt`.* This separates code regressions from the new numerics.

### Step 2: Kimmel 2D with order blending (first new result)
4. **Kimmel 2D P1, implicit + OB.** Start from Ray's `kimmelP1`.
   - Settings: `SpaceRHSJacobCom = RHSJacobGammaAlpha`, `ConvSolveCom = ConvRHSJacobTurb` (becomes `…TurbBlending`), `DiffSolveCom = DiffRHSJacobTurb`, `UseBlending = true`, `PreProcessCom = OrderBlending`.
   - Take `OrderBlending.*` from the wedge or Gnoffo templates: S0, Kappa, monitored `rho` or `p`.
   - Turn off LLAV and keep the analytical source Jacobian.
5. Same at **P2**, then **P3**: upgrade from the P1 / P2 OB results, or start from Ray's `kimmelP2New` / `kimmelP3`.
6. **Explicit + OB** at P1–P3 (`ConvSolveCom = ConvRHS` + `UseBlending`). Note M5.
7. **Record per run:** iterations, wall time, memory, residual drop, Stanton versus Ray (LLAV). This is the comparison for the FWO report.

### Step 3: Kimmel 2D with JFNK
8. Same as step 4 with `NewtonIteratorLSS.SetupCom = ParJFSetup`, `SysSolver = ParJFSolveSys`, `PCType = PCSHELL`, `ShellPreconditioner = FRP0Precond`, `Data.DoComputeJacobian = false` (from the wedge JFNK-OB template).
   - P1, then P2–P3, with and without OB.
   - If it stalls, see M6.

### Step 4: DiCristina 2D (P1–P2)
9. Repeat steps 1–3 with `RAY_dicristina2D_GA_P1_expl` / `P2_expl` (Ray's solutions inside), with Tu 0.75 % as in Ray's production runs.

### Step 5: Kimmel 3D (P1–P3), best 2D settings (VSC)
10. On the **3240-hex tip sector** (a complete Ray ladder exists):
    - P1 from `Res_cone/3D/kimmel3DP1.CFmesh`;
    - P2 = `RAY_kimmel3D_GA_P2_expl`;
    - P3 = `RAYPROD_kimmel3D_GA_P3bis_expl`.
    - Compare with Ray's `kimmel3DP2/run2` and `05_reference_results_RAY/Res_cone/Conv/3D30DegGAP2.plt`.

### Step 6: DiCristina 3D (P1–P2), best 2D settings (VSC)
11. Choose a mesh: 59,904 hexa (`RAY_dicristina3D_GA_P1_impl`, `dicriP1BeforeCrash`) or 79,248 (`RAY_dicristina3D_GA_P2_expl`). Then run P1 → P2.

### In parallel, low effort
- Finish the T3A references: LM P5, γ-α P2/P5, Jatayu.
- Fetch the experimental data (M2).

---

## 5. TODO to fulfil the project proposal

Status:
- ✅ done
- 🟡 partial (what exists)
- ⬜ not started

**→** marks the next concrete action.

### WP3: models
- 🟡 **A1. Tw-dependent spot growth.** `KarschSpotGrowth` exists (2D, uncommitted, constant Tw, numerical Jacobian only).
  → Use the local Tw/T0, port it to 3D, add the Jacobian, commit.
- ⬜ **A2. Verify A1 at several Tw/T0.**
  → Kimmel (0.42), Kimmel CTR (Tw 320 K), DiCristina (0.27), HEXAFLY (Tw 850 K).
- ⬜ **A3. Low-Tu onset correlation.** Tu is only floored at 0.027 %.
  → Implement it behind an option, in 2D and 3D.
- ⬜ **A4. Tu sweep.**
  → Ray's DiCristina Tu 0.35–1 % series is the baseline.
- 🟡 **A5. SA-γ-α.** SA FR exists in 2D only, with a numerical Jacobian.
  → Variable sets, source term (2D and 3D) with an analytical Jacobian, SA 3D.
- ⬜ **A6. Compare SA-γ-α, SST-γ-α and improved SST-γ-α.**
- 🟡 **A7. LM as the reference model.** LM T3A P0/P1 reproduced on the current build (§3).
  → Hypersonic LM Kimmel (`RAY_kimmel2D_LM_*`, 3D `RAY_kimmel3D_LM_*`).

### WP3: 3D numerics
- 🟡 **B1. RANS-LTT on tets, prisms and hybrid meshes.** FR Tetra (≤ P5) and Prism (≤ P4) exist in the current code. All RANS-LTT runs so far are on hexahedra.
  → Exact 3D wall distance on curved faces; 3D subsonic turbulent BCs; pyramids or pyramid-free meshes; Kimmel 3D on a prism-layer + tet mesh against Ray's hexa results.
- 🟡 **B2. RANS-LTT up to P9.** Highest so far: 3D P4 (HIFiRE, hexa), 2D P5 (Kimmel), 2D P7 LM flat plate (case only).
  → Extend Tetra/Prism to P9 and geometric shape functions above P2; **OB as the high-P shock capturing** (Andrea's list is the first step).
- ⬜ **B3. P-refinement P1 → P9 on a hybrid mesh.**

### WP3: validation
- 🟡 **C1. HIFiRE-5 versus flight.** Ray: M 2.64, P1–P4, P4 solution stored, CFcases recovered (`RAYPROD_hifire3D_*`), restarts lost, no flight data on disk.
  → Restart from the P4 solution or the Tu 0.3 % P1 solution; get the flight data.
- 🟡 **C2. DiCristina versus experiment.** Ray: 2D P1–P3, 3D P1–P3 (not residual-converged).
  → Andrea items 2 and 4, plus the experimental data.
- 🟡 **C3. Complex vehicle (HEXAFLY-INT).** Setup only (`RAY_hexafly3D_GA_P1_expl`, Q1 hexa).
  → Run P1 once the 3D settings are fixed.
- ⬜ **C4. Datasets with measured free-stream disturbance level.** Literature task.
- 🟡 **C5. Coarse FR versus FV reference.** FV LM Kimmel cone (36,400 cells) against FR (1764).
  → Compare Stanton; there is no FV γ-α.

### WP4: optimisation
- ⬜ **D1. 3D r-AMR coupled to RANS-LTT.** MeshFitting FR is 2D only and has no turbulence coupling.
- 🟡 **D2. Matrix-free implicit FR.** The infrastructure exists.
  → **Andrea items 1–4 with JFNK** (step 3), then a preconditioner with the source Jacobian (M6).
- ⬜ **D3. Lightweight Jacobian fallback.** `BlockDecoupledJacob` and `SinkOnlyDiagonal` are first steps.
- ⬜ **D4. Memory and time versus the full Jacobian above P3–P4.**
  → Comes free from logging steps 2–3 (step 2, item 7).
- ⬜ **D5. Trajectory ATD analysis.** HEXAFLY is the candidate.

### Project level
- ⬜ **E1. TCNEQ + RANS-LTT:** decide.
- 🟡 **E2. Meshes.** Gmsh hexa meshes and sources exist for Kimmel, DiCristina, HIFiRE and HEXAFLY. Nothing for tets/prisms.
- ⬜ **E3. Papers.**
  → The first natural paper is "Order blending and JFNK for hypersonic RANS-LTT FR" (Andrea's list + D2/D4).
- 🟡 **E4. Open-source merge.** The 2021–22 work is in master; the 2025–26 work is uncommitted (M7).

---

## 6. Decisions needed

1. **Kimmel 3D mesh for item 3:** the 3240-hex tip sector (complete Ray ladder, recommended) or the 19,440-hex 1° AoA half cone (P1 only)?
2. **DiCristina 3D mesh for item 4:** 59,904 or 79,248 hexa?
3. **Explicit + OB:** accept the generic boundary fallback for a first look, or implement the γ-α-specific explicit blending command first?
4. **Reference Tu:** Kimmel and DiCristina (Ray's production runs use 0.75 % for DiCristina); HIFiRE 0.3 % or 0.4 %?
