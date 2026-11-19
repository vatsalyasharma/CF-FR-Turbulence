# TODO: RANS-LTT FR development (status as of 2026-11-19, revised after the WORK-HD audit)

Legend:
- `[x]` done
- `[~]` partially done (what exists is stated)
- `[ ]` not started

Evidence is in STATUS.md (§1.8 covers Ray's archive). Code paths are relative to `COOLFluiD/plugins/`. Case folders are relative to `Validation_TURB/`.

## ★ Near-term cases to try (user list, 2026-11-19)

The "best settings" are to be picked from the 2D experiments. Order blending exists in the current FR (`PreProcessCom = OrderBlending`, `UseBlending = true`; RD 2026-10). The only RANS-LTT case already using it is `02_validation_baseline/SALV_kimmel_GA-SST_P2_implOB` (unreliable origin).

- [ ] **N1. Kimmel 2D (up to P3), with order blending: explicit, implicit, and implicit with or without matrix-free.**
  - Start from Ray's chain: `02_validation_baseline/RAY_kimmel2D_GA_P0_impl` → `P1_expl` → `P2_expl` → `P3_expl`. Ray's P1–P3 solutions are inside the folders; his P0 Newton run reached a ρ drop to −4.85.
  - [ ] Add the `OrderBlending` settings to each order (replacing or alongside LLAV); keep the LLAV runs as the reference.
  - [ ] Explicit: FwdEuler. Implicit: Newton, `RHSJacobGammaAlpha` + `ConvDiffLLAVRHSJacobTurb`, analytical source Jacobian. Matrix-free: PETSc `ParJF*`/`ParMF*` with `FRP0Preconditioner` (see M1).
  - [ ] Compare the wall Stanton/heat flux with Ray (`05_reference_results_RAY/PLT_cone/kimmelP*.plt`), iterations, wall time and memory.
- [ ] **N2. DiCristina 2D (up to P2), same matrix as N1.**
  - Start from `02_validation_baseline/RAY_dicristina2D_GA_P1_expl` → `P2_expl`.
  - Ray's Tu and spot-growth sweeps are in `05_reference_results_RAY/Res_cone/Dicristina/`.
- [ ] **N3. Kimmel 3D (up to P3), with the best settings from N1.**
  - Inputs: `02_validation_baseline/RAYPROD_kimmel3D_GA_P1_expl` (Ray's production case, runnable: 19,440 hexa, P1 restart `3DP1AoA1RunRestart4`) and `RAYPROD_kimmel3D_GA_P3bis_expl` (runnable: 3240 hexa, P3 restart `run2`). P2 and P3 have no restart. The older variants are the `RAY_kimmel3D_GA_P*` folders.
  - Reference solutions: `05_reference_results_RAY/Res_cone/3D/`.
  - [ ] Decide the mesh: the 3240-hex tip sector (stored CFcases) or the 19,440-hex half cone at 1° AoA (as-run).
  - [ ] Port the configuration from `AS_RUN_on_VSC_*.txt`.
- [ ] **N4. DiCristina 3D (up to P2), with the best settings from N2.**
  - Inputs: `02_validation_baseline/RAY_dicristina3D_GA_P1_impl` (P1 restart on 59,904 hexa) and `RAY_dicristina3D_GA_P2_expl` (P1 solution on 79,248 hexa, `dicristinaNewP1Run4`).
  - Production configurations and `.inter` files: `RAYPROD_dicristina3D_GA_P1..P3_expl__MISSING_RESTART`.
  - Extra P1 solutions: `05_reference_results_RAY/Res_cone/dicriP1BeforeCrash*.CFmesh`.
  - [ ] Choose one 3D mesh: `Bis` (79,248) or the 59,904 mesh.

## ★ M. Matrix-free implicit FR (priority)

This is linked to D2/D3/D4 below. It is the main lever on cost: Ray's explicit 3D runs needed 0.2–3 M iterations on 72–288 ranks and never converged by residual.

- [~] **M1. Get JFNK running on FR turbulence (2D first).**
  - **Exists:** `Petsc/ParJFSetup`, `ParJFSolveSys(GMRESR)`, `ParMFSetup`, `ParMFSolveSys`, Seq variants, and `FluxReconstructionPetsc/FRP0Preconditioner` (p-multigrid shell preconditioner "for JFNK in FR"), `FRBlockJacobiPreconditioner`, `FRJacobAssembler` (RD 2026-03 to 05).
  - **Never used with** `*Turb` / γ-α commands; no case exists.
  - [ ] Build a 2D Kimmel γ-α P1–P3 JFNK case (part of N1).
- [ ] **M2. Preconditioner keeping the turbulence/transition source-term Jacobian.**
  - Use the analytical block-decoupled Jacobian of `GammaAlpha{2D,3D}SourceTerm` and the AD Jacobian of `NavierStokesGReKO*SourceTerm_Lang`.
  - Put it inside the FR block-Jacobi or P0 preconditioner (`FluxReconstructionPetsc/FRBlockJacobiPreconditioner.*`).
- [ ] **M3. Positivity and robustness of k, ω, γ, α under JFNK.** Includes the perturbation size for the log ω variable and clipping (`PhysicalityGReKO`, `TurbUpdateSol`).
- [ ] **M4. Memory and time comparison** (matrix-free versus assembled Jacobian versus explicit) at P2–P4, in 2D and 3D. This is the same as D4.

## 0. Prerequisites found during the audit

- [ ] **P0. Commit or merge the uncommitted turbulence work** on `fr-turbulence`.
  - The work covers: γ-α ModelVersion 2025, Karsch spot growth, onset table, `PhysicalityGReKO`, `FRTurbProbes`, `LM2009Source.hh`, `TransitionSourceScales.hh`, etc.
  - Reconcile with CF-FR-Turbulence (`~/Desktop/MERGE_CF-FR-Turbulence/REPORT.md`), `COOLFluiD_VIN_26` and `CF_RAY_Version`.
- [ ] **P1. Re-establish Ray's 2D baseline with the current build.**
  - Run Ray's Kimmel 2D P0 → P4 (not Salvatore's) and compare with Ray's stored `kimmelP1..P5` solutions and plots.
- [x] **P2. Use the as-run configurations.** Done by finding Ray's production CFcases and `.inter` files in `~/Desktop/CF_RAY_Version` (matched to the VSC job logs). They are copied as `02_validation_baseline/RAYPROD_*`.
  - [ ] Decide how to restart the cases whose production restart is lost: start from the upgrade input, or use a stored alternative (HIFiRE P4 → `hifireTu04P4Run4`; DiCristina 3D P1 → `dicriP1BeforeCrash`).
- [ ] **P2b. Build Ray's code** (`CF_RAY_Version/COOLFluiD_RAY`, not built; `Install_coolfluid_RAY.sh` is there) and rerun `RAYPROD_kimmel3D_GA_P1_expl` and `RAYPROD_kimmel2D_GA_P2-CTR_expl` with both Ray's build and the current build. This separates code regressions from setup differences.
- [ ] **P3. Implicit robustness for hypersonic γ-α and LM.**
  - All implicit runs (Ray and Salvatore) stall, except Ray's Kimmel 2D P0.
  - LM Kimmel failed early in Salvatore's runs; check it with Ray's originals (`01_models/RAY_kimmel2D_LM_*`).
- [ ] **P4. One free-stream Tu per validation case.** Ray swept Tu (DiCristina 0.35–1 %, HIFiRE `Tu04`); fix the reference values and their source.
- [ ] **P5. Recover the lost restarts if possible:** HIFiRE P1–P3, DiCristina 3D P1–P3, Kimmel 3D P2/P3 and Kimmel 2D P0. Check VSC `$VSC_DATA/CPU/COOLFluiD_VSC/plugins/GReKO/testcases/Cone`, if the account still exists. Otherwise regenerate them.
- [ ] **P6. Regression testcases** in `FluxReconstructionTurb/testcases/` (currently none): γ-α T3A P1 (meshes now located in Ray's FlatPlate folder), Kimmel 2D P1, Kimmel 3D P1.

## A. WP3: γ-α model improvements

- [~] **A1. Wall-temperature-dependent spot growth.**
  - **Exists:** `KarschSpotGrowth` / `KarschWallTemperature` / `KarschRecoveryFactor` in `FluxReconstructionTurb/GammaAlpha2DSourceTerm.{hh,cxx}` (uncommitted, 2D only, `ModelVersion = 2025`, numerical Jacobian only, constant Tw).
  - **Related (2021):** Ray's spot-growth sweeps (`Res_cone/GGrow`, DiCristina `ns`/`CEg`/`Cmg` series), which tuned constants rather than adding a Tw dependence.
  - [ ] Use the local Tw (or Tw/T0) instead of the constant.
  - [ ] Port to `GammaAlpha3DSourceTerm`.
  - [ ] Add the analytical Jacobian.
  - [ ] Commit.
- [ ] **A2. Verify A1 at several Tw/T0.**
  - Available: Kimmel (Tw/T0 ≈ 0.42), Kimmel CTR variant (Tw = 320 K), DiCristina (≈ 0.27), HIFiRE (Tw = 332 K, M 2.64), HEXAFLY (Tw = 850 K).
- [ ] **A3. Low-Tu onset correlation.** Currently Tu is only floored at 0.027 % (`TransitionSourceScales.hh`, `getTuInfLocal`; the onset correlation is in `GammaAlpha2DSourceTerm.cxx` about lines 1365–1393 and the 3D counterpart).
- [ ] **A4. Verify A3** with a Tu sweep. Ray's DiCristina Tu 0.35–1 % series is a ready baseline.
- [~] **A5. SA-γ-α in HO FR.**
  - **Exists:** SA FR, 2D only (`FluxReconstructionTurb/SA2DSourceTerm`, numerical Jacobian).
  - [ ] SA-γ-α variable sets.
  - [ ] `SAGammaAlpha{2D,3D}SourceTerm`, with the analytical Jacobian.
  - [ ] `SA3DSourceTerm` in FR.
  - [ ] BCs.
- [ ] **A6. Compare SA-γ-α, SST-γ-α and improved SST-γ-α.**
- [~] **A7. Keep LM as the reference.**
  - **Exists:** LM FR 2D and 3D. Ray ran LM in 3D on Kimmel (P1/P2; solutions in `Res_cone/3D`, plots `3D30DegLMP*`). T3A flat plates P0–P7.
  - [ ] Make LM run with the current build.
  - [ ] Add LM to the regression set.

## B. WP3: extension to the 3D numerics developed for TCNEQ

- [ ] **B1. RANS-LTT on tets, prisms and hybrid meshes.**
  - **Exists:** FR Tetra (≤ P5) and Prism (≤ P4) element data in the current main tree; not in Ray's code. All RANS-LTT 3D runs so far are hexahedral.
  - [ ] Exact 3D wall distance to curved Q2 faces (`MeshTools/ComputeWallDistanceVectorFRMPI.cxx`).
  - [ ] 3D subsonic turbulent BCs.
  - [ ] Pyramid support, or meshes without pyramids.
  - [ ] Kimmel 3D on prism layer plus tets, compared with Ray's hexa results.
- [ ] **B2. RANS-LTT up to P9.**
  - Extend Tetra/Prism to P9; add shape functions above P2.
  - Positivity at high P; interaction with LLAV and order blending (see N1).
  - **Highest RANS-LTT orders so far:** 3D P4 (HIFiRE, hexa) and 2D P5 (Kimmel).
- [ ] **B3. P-refinement P1 → P9 on a hybrid mesh.**
  - Available hexa ladders: Kimmel 3D P0–P3, HIFiRE P1–P4 (configs), Kimmel 2D P0–P5.

## C. WP3: validation

- [~] **C1. HIFiRE-5 versus flight.**
  - **Exists (Ray, 2022):** a HIFiRE elliptic-cone case at **M 2.64**, Tw 332 K, P1–P4 on 98,496 hexa. A P4 solution is stored; there is a P1 CFcase only.
  - **Missing:** residual convergence, the flight comparison data, and the P2–P4 CFcases (configs only).
  - [ ] Confirm the trajectory point and Tu.
  - [ ] Collect the flight data.
- [~] **C2. DiCristina cone (Mach 10.2).**
  - **Exists:** Ray's 2D P1–P3 and 3D P1–P3 runs (3D not residual-converged), with Tu and spot-growth sweeps.
  - **Missing:** experimental data on the drive; most 3D restarts.
- [~] **C3. Complex 3D vehicle (HEXAFLY-INT).**
  - **Exists:** Ray's setup case `hexaflyFRP1` (M 7.03, Tw 850 K, 161k Q1 hexa, P1) and the Gmsh mesh `gmsh4test2vol`. Never run [I].
- [ ] **C4. Validation datasets with measured free-stream disturbance levels.** Nothing on the drive.
- [~] **C5. Mesh-coarsening comparison versus an FV reference.**
  - **Exists:** FV LM Kimmel cone (36,400 cells) against the FR Kimmel mesh (1764 cells). There is no FV γ-α.

## D. WP4: optimization (LTT share)

- [ ] **D1. r-AMR in 3D, coupled to RANS-LTT.**
  - `FluxReconstructionMethod/MeshFittingAlgorithmFR*` is 2D only, also in Ray's code. Wall-distance recomputation is needed after mesh motion.
- [~] **D2. Matrix-free implicit FR.** → see **★ M** (priority).
- [ ] **D3. Fallback: lightweight Jacobian.** `BlockDecoupledJacob` and `SinkOnlyDiagonal` are first steps.
- [ ] **D4. Memory and runtime versus full Jacobian above P3–P4.** → M4.
- [ ] **D5. Flight-trajectory ATD of a full 3D vehicle.** The HEXAFLY setup (C3) could be the starting point.

## E. Project-level

- [~] **E1. TCNEQ + RANS-LTT.** Not implemented; a decision is required.
- [~] **E2. Mesh generation.**
  - **Exists:** Ray's Gmsh hexa meshes and sources for Kimmel 3D, DiCristina 3D, HIFiRE (`.geo` incl. a hexa version) and HEXAFLY (`05_reference_results_RAY/MSH_sources`, `HiFire_geometry_Firas`). Nothing for tets/prisms.
- [ ] **E3. Publications.** None found.
- [~] **E4. Merged into open-source COOLFluiD.**
  - The 2021–22 FR turbulence work is in master.
  - The 2025–26 work is uncommitted or external.
