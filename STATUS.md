# STATUS: COOLFluiD HO FR RANS-LTT solver (audit of 2026-11-19)

Scope:
- **First pass:** `/home/u0142014/Desktop/Cases-Salvatore` (S1) and `/home/u0142014/Desktop/COOLFluiD/plugins` (S2).
- **Second pass (same day):** `/media/u0127536/WORK-HD/RAY-FIRAS` (S3) and `/media/u0127536/WORK-HD/RAY_VSC_FILES` (S4), which hold Ray Vandenhoeck's original cases, results and VSC job logs.
- Read-only inspection only; nothing was run or built.

> **Revision note (second pass).**
> - S1 (Salvatore) is now treated as **unreliable**. Its cases are modified copies of Ray's 2D cases (S3): 20–57 active lines differ per case.
> - S3/S4 contain the 3D work missing from S1/S2: Kimmel 3D, DiCristina 3D, HIFiRE and HEXAFLY, all on Q2/Q1 hexahedra.
> - The baseline table (§2) has been revised accordingly; §1.8 and §3 items 14 onwards are new.

Each claim is tagged:
- **[V]** verified in code, files or logs.
- **[I]** inferred.

---

## 1. Inventory

### 1.1 Repositories and git state

**Main tree: `/home/u0142014/Desktop/COOLFluiD`** [V]
- **Branch:** `fr-turbulence`, which has no commits beyond `master`.
- **HEAD:** `f92e6297` (2026-11-04, A. Lani, MHD testcase).
- **Last commit touching FR turbulence:** `f6c7239a` (2026-09-21, RD), "FR diffusive path: compact BR2 … fixed Jacobian assembly".
- **Commits before that:** 2021–2022 (Ray Vandenhoeck). Examples:
  - `556a545b` "FR 3D gamma-alpha ST Jacobian"
  - `56c0a4d8` "3D GReKLogO FR"
  - `7d3b975c` "3D source term, BC, transition criterion for FR gamma-alpha"
  - `1aa27669` "SST-V and SST-Vm"
- **Uncommitted work:** 84 files modified, 826 untracked, 747 deleted. The deletions are almost all MHD/COCONUT `.xz` testcase files.
  - **Modified:** most of `FluxReconstructionTurb`, `GammaAlpha` and `GReKO`. For example, `GammaAlpha2DSourceTerm.cxx` has +1028 lines.
  - **Untracked turbulence pieces:** the `FRTurbProbes/` plugin and `GammaAlphaOnsetTable.hh`, `TransitionSourceScales.hh`, `AUSMPlusUpFluxGammaAlpha2D`, `PhysicalityGReKO`, `GReKO/LM2009Source.hh`, `GammaAlpha{ArtificialDiffusion,EGammaFlux,MuTFactor}.hh`, `NewtonMethod/PhysicalityCFLUpdateSolTurb.*`, `MeshTools/FiniteWallEdge.hh`.
  - **⚠ All of the "ModelVersion 2025" γ-α work, including Karsch spot growth and the onset table, exists only as uncommitted or untracked files.**
- **Build:** `OPENMPI/optim` was built 2026-11-19 14:16. It contains every library the copied cases need [V].

**Other trees on the Desktop** [V]
- `COOLFluiD_VIN_26` (master, 2026-09-08): its FluxReconstructionTurb has 46 differences from the main tree.
- `CF_RAY_Version/COOLFluiD_RAY` (not a git repo): also 46 differences.
- `MERGE_CF-FR-Turbulence/`: a merge kit for the external repo `github.com/vatsalyasharma/CF-FR-Turbulence`. It contains `REPORT.md`, patches, and `T3A_GA_chain.tgz`. The external repo has T3A flat-plate γ-Reθ (LM2009) and SST γ-α cases at P0–P5 *with meshes* [V from REPORT.md].
- These trees are outside the audit scope; I only listed them.

**Trees named in Salvatore's notes:** `COOLFluiD_VSC`, `COOLFluiD-merged-clean` and `COOLFluid` are not on this machine. His runs used `/home/r1132394/COOLFluid` (the symlinks in S1 are broken) [V].

### 1.2 Source code: turbulence and transition models

| Model | FR implementation (`plugins/FluxReconstructionTurb`, lib in brackets) | FV implementation | Physics plugin |
|---|---|---|---|
| SST/BSL k-logω (+SST-V, SST-Vm) | `KLogOmega2DSourceTerm`, `KLogOmega3DSourceTerm` (libFluxReconstructionKOmega) | none for k-logω | `KOmega` |
| k-ω Wilcox | `KOmega2DSourceTerm` (2D only) | `FiniteVolumeTurb/NavierStokesKOmega{2D,3D}SourceTerm` | `KOmega` |
| Langtry-Menter γ-Reθ | `NavierStokesGReKO2DSourceTerm_Lang`, `NavierStokesGReKO3DSourceTerm_Lang` (libFluxReconstructionGReKO) | `FiniteVolumeGReKO/NavierStokesGReKO2DSourceTerm_Lang*` (2D only) | `GReKO` (`NavierStokes{2D,3D}GReKLogO{,BSL,SST}…`) |
| γ-α (VdE & Steelant) | `GammaAlpha2DSourceTerm`, `GammaAlpha3DSourceTerm`, `GammaAlphaTransitionCriterion` (libFluxReconstructionGammaAlpha) | **none** | `GammaAlpha` (`NavierStokes{2D,3D}GammaAlpha{,BSL,SST}…`) |
| Spalart-Allmaras | `SA2DSourceTerm` (2D only; libFluxReconstructionSA) | `FiniteVolumeTurb/SA{2D,3D}SourceTerm`, DES variants | `SA` |
| SA-γ-α | **none** | none | none |

**FR commands, BCs, fluxes and tools** [V]
- **Volume and diffusive terms:** `ConvDiffRHSJacobTurb`, `ConvDiffLLAVRHSTurb`, `ConvDiffLLAVRHSJacobTurb`, `DiffRHS(Jacob)Turb`.
- **Boundary corrections:** `Diff/ConvBndCorrectionsRHS(Jacob)GammaAlpha`, including `…Blending` and `…SubcellBlending`.
- **Wall BCs:** `BCNoSlipWallGammaAlpha{2D,3D}`, `BCNoSlipWallTurb{2D,3D}`.
- **Inflow/outflow BCs:** `BCSubInletTurbUVT2D` and `BCSubOutletTurb2D`. **There are no 3D turbulent subsonic inlet/outlet BCs.**
- **Riemann flux:** `AUSMPlusUpFluxGammaAlpha2D`.
- **Wall distance:** `ComputeWallDistanceVectorFRMPI` (MeshTools).

**γ-α onset correlation** (`GammaAlpha2DSourceTerm.cxx` ~1365-1393) [V]
- `Re_θt = 0.664·sqrt(400094·Tu^-1.38 − 105254·Tu^-7/8)`, with the local Tu recovered with a Mach factor and clipped to [0.027, 13] %.
- There is an optional pressure-gradient variant (`PGrad`).
- The 2025 version adds an onset function (wall criterion τw/√(ρμ)) and an onset table keyed on (x, y), which is 2D only.

**LM onset correlation** [V]
- The standard LM 2009 correlation (Tu ≤ 1.3 / > 1.3 branches), with F(λ) and Re_θt ≥ 20.
- In both models, low Tu is only floored at 0.027 %.

**γ-α spot growth** [V]
- **Default (2018):** `nσ ∝ Tu^(7/4)·f_k·f_Mnσ` with `f_Mnσ = 1/(1+0.58 M^0.6)²`. It depends on Mach only, **not on Tw**.
- **`KarschSpotGrowth` option (uncommitted, 2D only, `ModelVersion = 2025`):** Karsch f_σ (Jacobs et al. 2025, eqs. 10-16), which uses Tw/Te.
  - Tw is a **constant input** (`KarschWallTemperature`), not the local wall temperature.
  - `ModelVersion = 2025` requires `AnalyticalJacob = false` (the code checks this at about line 238).
  - **No case uses it.**

### 1.3 Element types and orders (FluxReconstructionMethod) [V]

**Element data and maximum solution order**

| Element | Max solution order | Source of the limit |
|---|---|---|
| Quad, Hexa | P10 | quadrature tables |
| Triag | P9 | quadrature tables |
| Tetra | P5 | quadrature tables 0–5 |
| Prism | P4 | element data reaches P5, but `SetPrismLagrange*StateCoord` stops at P4 |
| Pyramid | none | `MeshUpgradeBuilder` throws NotImplemented |

**Geometry**
- Q1/Q2 geometry is supported for all five element types.
- The Framework has P2 Lagrange shape functions for Hexa (20 and 27 nodes), Prism and Tetra.
- **There are no geometric shape functions above P2.**

**Turbulence code and element types**
- The FR turbulence code has no element-type-specific branches. It is split by dimension into separate 2D and 3D classes.
- **Missing in 3D:** SA, k-ω and the subsonic turbulent BCs are 2D only. Every γ-α 2025 option is also 2D only.

**Wall distance** (`ComputeWallDistanceVectorFRMPI`) [V]
- It loops over solution points, so it is element-generic.
- In 2D it computes the exact distance to finite P1/P2 edges.
- In 3D it **projects onto the plane of the nearest face centre**. The code comment says curved finite surfaces are not implemented. Near curved walls in 3D this is approximate [I: can matter for prism-layer meshes].

### 1.4 Jacobian and implicit [V]

- **Mechanism:** `StdSourceTerm` has `AddJacob` and `AnalyticalJacob`. The numerical Jacobian perturbs the full cell block. The Jacobian can be frozen with `FreezeJacob`/`FreezeJacobInterval`.
- **Analytical source Jacobian by model:**
  - **KLogOmega 2D/3D:** analytical, plus a gradient Jacobian.
  - **γ-α 2D/3D:** hand-coded analytical, all equations, with an option for a block-decoupled Jacobian.
  - **LM 2D/3D:** analytical by forward-mode AD (`lm2009::Scalar`), except that the axisymmetric path uses legacy hand-coded terms. Extra options: `BlockDecoupledJacob` and `SinkOnlyDiagonal`.
  - **SA2D and KOmega2D:** **numerical only.** With `AnalyticalJacob = true` they contribute nothing [I].
- **γ-α 2025 options:** these force the numerical Jacobian (`GradientJacobian`, `EGammaImplicitDiagonal`).
- **Linear solver:** PETSc Newton with GMRES, preconditioned by ASM or BJacobi.

### 1.5 r-AMR and matrix-free [V]

- **r-AMR:** `MeshFittingAlgorithmFR`, `…FRQ2`, `…FRQ2Quads` and `StdMeshFittingUpdateFR` (spring analogy, F. Ben Ameur, 2022).
  - It is **2D only**: the 3D Hexa/Tet assembly calls are commented out.
  - There are no turbulence references and no testcase uses it.
- **Matrix-free / JFNK:**
  - PETSc provides `ParJFSetup`/`ParJFSolveSys`, `ParMFSetup`/`ParMFSolveSys` and the Seq variants.
  - FluxReconstructionPetsc adds `FRP0Preconditioner` (p-multigrid shell preconditioner "for JFNK in FR"), `FRBlockJacobiPreconditioner` and `FRJacobAssembler` (RD, 2026-03 to 05).
  - **No testcase combines JFNK with FR turbulence.**

### 1.6 TCNEQ + turbulence [V]

- `FluxReconstructionNEQ` has no turbulence.
- `NEQKOmega` exists only as physics variable sets, with no FR or FV source term and no case.

### 1.7 Cases using a turbulence or transition model

**S1 (Cases-Salvatore): 29 CFcases, all FR, all 2D axisymmetric, all Q2 quads** [V]

All S1 cases load the same libraries. Each run upgrades a restart in a separate mesh subsystem before the FR solve.

| Group | Model | P | Conv. method | Mach | Mesh/restart present | Results / convergence (log10 residual of ρ) |
|---|---|---|---|---|---|---|
| Kimmel GA `Cases/kimmelFRP0` | GA-SST | P0 | Newton | 7.93 | yes | `Res_implicit` 6470 it, stalled 1.0 → 0.98 |
| Kimmel GA `kimmelFRP1`, `P1bis` | GA | P1 | FwdEuler | 7.93 | yes | 7.29M it in 14 h on 32 ranks, final 1.12; `bis` 1.9M it, final 0.34 |
| Kimmel GA `kimmelFRP2`, `P2bis`, `P2-implTry` | GA | P2 | FwdEuler / Newton | 7.93 | yes | 307k it, −0.65 (best); `bis` 0.53; `implTry` 5600 it, 3.16 (14 negative k) |
| Kimmel GA `kimmelFRP3`, `P3bis`, `P3-impl` | GA | P3 | FwdEuler / Newton | 7.93 | yes | 750k it (MaxIter), 0.46; `bis` no run phase; `impl` 166 it |
| Kimmel GA `kimmelFRP4`, `P4-impl` | GA | P4 | FwdEuler / Newton | 7.93 | yes | 150 it; `impl` 3526 it, 1.5 → 0.6 |
| Kimmel GA `kimmelFRP5` | GA | P5 | FwdEuler, MaxIter = 1 | 7.93 | yes | no residual history |
| Kimmel LM `kimmelFRGReKOP0..P3` | LM-SST | P0–P3 | Newton (P0), FwdEuler | 7.93 | yes | P0: 2 it; P1: 3276 it (3.2 → 3.1); P2: 50 it; P3: 41 it |
| DiCristina GA `dicristinaFRP1`, `P2` | GA | P1, P2 | FwdEuler / Newton | 10.2 | yes | P1 150k it, 4.1 → 0.5; P2 18100 it, ~3.5 |
| `Usable_Cases` (kimmel2D P0–P4, diCristina2D P1–P2) | GA | P0–P4 | mixed | 7.93 / 10.2 | several restarts missing (see MANIFEST); CFcases otherwise identical to `Cases/` except P0 and P2 | no results stored |

**S2 (repo testcases)** [V]

There are no `.CFcase` files in FluxReconstructionTurb, GammaAlpha, FiniteVolumeGReKO, FiniteVolumeTurb or FRTurbProbes.

| Path | Model | Method | Dim | P | Mach | Mesh |
|---|---|---|---|---|---|---|
| `GReKO/testcases/FlatPlate/flatPlateFRImplGReKLogO_T3A_Coarse_P0..P7` | LM | FR Newton | 2D | P0–P7 | 0.0158 | present except P0 and P5 inputs (and P4's direct `T3AP4New`) |
| `…_T3A_Tu35_P0..P6` | LM | FR | 2D | P0–P6 | 0.0158 | P1–P5 present |
| `…flatPlateFRImplGAlphaKLogO_T3A_P0..P2` | **GA** | FR | 2D | P0–P2 | 0.0158 | **all missing** |
| `GReKO/testcases/Cone/cone.CFcase` | LM (k-ω) | **FV** | 2D | – | 7.93 | present (36,400 cells) |
| `SA/testcases/FlatPlate/flatPlateFRImplSA_Re1e7_P0..P4`, `…Coarse P0..P6` | SA | FR | 2D | ≤ P6 | 0.2 | present (except P0Coarse) |
| `SA/testcases/{FlatPlate,Cylinder}` FV (2D and 3D) | SA | FV | 2D/3D | – | 0.1 | meshes absent |
| `KOmega/testcases/FlatPlate/flatPlateFRImplKOmega_*` | k-ω BSL | FR | 2D | P0 | 0.2–0.3 | present |

**No 3D RANS/RANS-LTT FR case exists in S1 or S2.** That includes HIFiRE-5, 3D cones, any tet/prism/hexa mesh, HEXAFLY-INT and Space Rider. **S3/S4 have them; see §1.8.**

### 1.8 WORK-HD: Ray's archive (S3 = RAY-FIRAS, S4 = RAY_VSC_FILES) [V unless marked]

**What the two folders contain**
- **S3** `RAY-FIRAS/CFTestCases/GReKO-GA/Cone` (13 GB) holds 39 CFcases dated 2020-05 to 2022-01, plus `CFmesh/`, `MSH/`, `Res/` and `PLT/`. `…/FlatPlate` (418 MB) holds 48 T3A/T3Am/T3B flat-plate cases (LM, LM-CTR, γ-α) with meshes.
- **S4** `RAY_VSC_FILES/vsc32241` holds:
  - `jobs/`: PBS scripts and **all stdout/stderr** of Ray's VSC Genius runs;
  - `COOLFluiD_Genius`: Ray's svn-era code. Compared with the main tree it has 15 differing files in FluxReconstructionTurb and 34 fewer files, and **no FR Tetra/Prism element data**. JFNK `ParJF*` already existed there.

**Cases found** (3D meshes are Q2 hexa unless noted)

| Case | Model | Dim / mesh | P | Mach | Status of what Ray ran (from job stdout) |
|---|---|---|---|---|---|
| Kimmel 2D | γ-α SST | 2D-axi quad (1764) | P0–P5 | 7.93 | P0 Newton: 1.03M it, ρ 1.8 → **−4.85** (best on the drive). P1–P5 solutions stored (`kimmelP1..P5.CFmesh`, `.plt`) |
| Kimmel 2D | LM | 2D-axi | P0–P3 | 7.93 | CFcases and some restarts only |
| Kimmel 3D | γ-α SST | 3240 hexa (tip sector) and 19,440 (`Dpdx-2Half`); "AoA1" restarts | P1–P3 | 7.93 | P1: 17 jobs, 3.2 → −0.75. P2: 2.8 → −0.06. P3: stuck around 1.5. All explicit, CFL ≤ 0.1, 72–108 ranks |
| Kimmel 3D | SST k-logω, LM, laminar NS (30°) | hexa | P0–P3 | 7.93 | solutions in `Res/3D`, plots `3D30Deg*` |
| DiCristina 2D | γ-α | quad | P1–P3 | 10.2 | Tu sweep 0.35–1 % and spot-growth sweeps (`ns`, `CEg`, `Cmg`) in `Res/Dicristina` |
| DiCristina 3D | γ-α | 79,248 hexa (`Bis`), 59,904 (`P1Run5`) | P1–P3 | 10.2 | P1: 22 jobs, 2.7 → −0.8. P2: stuck around 1.6. P3: ~2.0; 144–288 ranks |
| **HIFiRE** (elliptic cone, very likely HIFiRE-5) | γ-α | 98,496 hexa, full body, ~1.5° incidence | P1–**P4** | **2.64** | P1: ~1.5M it, ρ around 0.2. P2: around 0.6. P3: around 1.3–1.7. P4: around 1.6 (25k it). **P4 solution stored** (`hifireTu04P4Run4.CFmesh`, 2.7 GB). Only the P1 CFcase survives; P2–P4 exist only as echoed configs |
| HEXAFLY-INT | γ-α | 161,328 hexa, **Q1** | P1 | 7.03 | setup CFcase only (MaxIter = 1); no run output found [I: never run] |

**Reliability of S3/S4**
- The stored CFcases differ from what was run. The job stdout echoes the effective configuration, and it shows:
  - `MaxIter = 1e10` (stored files often have 1);
  - different CFL ramps;
  - different LLAV S0, κ and Péclet;
  - different interactive files and restarts.
- **Almost all `*.inter` files are lost** (only 3 survive).
- **Most production restarts are lost**: HIFiRE P1–P3, DiCristina 3D P1–P3, Kimmel 3D P2/P3 and Kimmel 2D P0. HIFiRE P4 and some Kimmel 3D solutions survive.

**Conclusions**
- **Convergence:** Ray's 3D γ-α runs were all explicit, and none converged by residual (best 3D drop about 4 orders). Validation was presumably judged on surface heat flux [I]. The `−1.79769e+308` "reached residual" also appears in Ray's runs, so it is the signature of a manual or interactive stop.
- **Thesis coverage:** this matches the thesis claims of "3D cones, HIFiRE, Mach 2.6–10.2, up to P4" in a weak sense. **3D at P4 was run only for HIFiRE**, on hexahedra, explicitly, not converged by residual.

### 1.9 `/home/u0142014/Desktop/CF_RAY_Version` (third pass) [V unless marked]

**What it is**
- `COOLFluiD_RAY/` is a source snapshot dated 2026-04-21, with install scripts (`Install_coolfluid_RAY.sh`, `COOLFluid_VSC_RAY.conf`, Salvatore's conf and `prepare.pl`).
- **It is not built** (`OPENMPI/` is empty).
- Its FR turbulence code differs from both Ray's GPU tree on WORK-HD (30 diffs in FluxReconstructionTurb) and the main tree (46 diffs). It has **no FR Tetra/Prism**.
- [I] It is most likely Ray's CPU tree `$VSC_DATA/CPU/COOLFluiD_VSC`, the one his PBS jobs used. It is therefore the code that produced his results.

**Cone testcases** (`plugins/GReKO/testcases/Cone`, 1 GB)
- **14 production CFcases:** Kimmel 2D P0 and P2-CTR, Kimmel 3D P1/P2/P3/P3bis, DiCristina 3D P1–P3, HIFiRE P1–P4 and P4Mesh, plus `cone.CFcase` (FV).
- **They match the run logs:** in every setting echoed by the VSC job logs, each CFcase is identical to what Ray ran (the residual differences are the CFEnv, Paths and mesh-upgrade blocks, which the logs don't print).
- **All 14 production `.inter` files are present.** They set LLAV κ, Péclet and S0, the wall laminar enforcement (`EnforceLaminar`, `GammaThreshold`, `LeadingEdgeLaminarEnforce`) and `FreezeLimiterIter`.
- **Solutions:**
  - `hifireTu03P1` (HIFiRE P1, Tu 0.3 %);
  - `dicriP1BeforeCrash{,Earlier}` (DiCristina 3D P1, 59,904 hexa);
  - `kimmel3DP2`, `run2` (Kimmel 3D P2/P3);
  - 2D Kimmel meshes, including `KimmelSingleConeRoundedBis` (1710 quads, used by the production 2D P0).
- **Still missing:** most production restarts (`kimmelGAKLogOCoarseP1-iter_*`, `3DP2AoA1Run3`, `3DP3Run6`, `dicristinaNew*`).
- A stray file literally named `\` is an orphan CFcase copy.

**FlatPlate**
- The γ-α T3A P0–P2 CFcases here are byte-identical to the main-repo testcases.
- Extra γ-α T3Am P2/P2bis and T3B P3 cases have their P2 inputs but not their restarts.

**Effect on §1.8 caveats**
- The caveat about lost interactive files is **resolved for the cone cases**.
- The caveat that stored CFcases differ from the run configuration is **resolved** by using `RAYPROD_*`.
- Lost restarts remain the main gap.

---

## 2. Baseline verification

| # | Item | Status | Evidence | Remarks |
|---|---|---|---|---|
| B1 | LM in HO FR, implicit | **Partially confirmed** | [V] `NavierStokesGReKO{2D,3D}SourceTerm_Lang` with AD Jacobian; `RHSJacobTurb`/`ConvDiffRHSJacobTurb`; T3A flat-plate FR Newton cases P0–P7 | No stored results for the flat plates. Hypersonic LM (Kimmel) runs failed or stopped early (2–3276 it). The 3D LM source exists but has no 3D case. |
| B2 | γ-α in HO FR, implicit | **Partially confirmed** | [V] `GammaAlpha{2D,3D}SourceTerm`; S1 Kimmel and DiCristina runs | The implicit runs stall or stop early (≤ 6470 it). Every long run in S1 was explicit, and none reached the stop criterion. The γ-α flat-plate repo cases lack meshes. |
| B3 | Source-term Jacobian for turbulence/transition | **Confirmed** (with gaps) | [V] analytical for KLogOmega, γ-α and LM (2D/3D); block-decoupled option; freezing | SA and KOmega have a numerical Jacobian only. The γ-α 2025 model forces the numerical Jacobian. Effectiveness is not demonstrated: the implicit runs converge poorly. |
| B4 | 3D elements: only Q2 hexa | **Confirmed for RANS-LTT runs; code has more** | [V] All of Ray's 3D RANS-LTT cases (S3) use Q2 hexa (HEXAFLY Q1 hexa). Ray's code (S4) has no FR Tetra/Prism. The current main tree has FR element data for Tetra and Prism | The baseline is correct for what was *run*. The tet/prism FR capability is newer and never used with turbulence. |
| B5 | Max order P4 in 3D | **Confirmed (HIFiRE only)** | [V] HIFiRE P4 jobs (about 25k it, ρ around 1.6) and the stored P4 solution `hifireTu04P4Run4.CFmesh`. Kimmel 3D and DiCristina 3D only to P3 (attempted) | P4 in 3D was reached only on HIFiRE, explicitly, not residual-converged. Code limits now: Hexa P10, Tetra P5, Prism P4. |
| B6 | Validated: HIFiRE-5, DiCristina, 3D cones, M 2.6–10.2 | **Partially confirmed** | [V] S3/S4: HIFiRE (M 2.64, P1–P4), DiCristina 3D (M 10.2, P1–P3), Kimmel 3D (M 7.93, P1–P3), with solutions and surface plots. Kimmel 2D P0 is the only run with a deep residual drop (−4.85) | **Not verifiable from the residuals:** no 3D run converged by residual. Experiment or flight comparison data was not found on the drive. HIFiRE is at **M 2.64**, not "Mach ~7" as in your table. Only the HIFiRE P1 CFcase survives (P2–P4 only as echoed configs). |
| B7 | SA-γ-α: not existing | **Confirmed** | [V] no SA-γ-α anywhere | SA itself exists in FR (2D only, numerical Jacobian, flat-plate cases); SA 3D FR is missing. |
| B8 | Tw-dependent spot growth: not existing | **Contradicts partially** | [V] `KarschSpotGrowth` + `KarschWallTemperature` in `GammaAlpha2DSourceTerm` (uncommitted) | 2D only, ModelVersion 2025 only, numerical Jacobian only. It takes a constant Tw rather than the local Tw/T0. It is untested (no case) and uncommitted. |
| B9 | Low-Tu onset correlation: not existing | **Confirmed** | [V] Tu floored at 0.027 %; onset correlation unchanged | The onset table and wall criterion (2025, 2D) are a different mechanism, not a low-Tu correlation. |
| B10 | Tets, prisms, hybrid with RANS-LTT: not verified | **Confirmed** | [V] no case | Code risks: 3D wall distance is approximate on curved faces; no pyramids, so hybrid hex/prism/tet meshes need pyramids or have to be avoided; there are no 3D subsonic turbulent BCs. |
| B11 | r-AMR coupled to RANS-LTT: not existing | **Confirmed** | [V] MeshFitting FR is 2D, with the 3D code commented out and no turbulence coupling | |
| B12 | Matrix-free FR: not existing | **Contradicts partially** | [V] PETSc JF/MF solvers (already present in Ray's code, S4) and `FRP0Preconditioner` "for JFNK in FR" (2026, main tree) | The infrastructure exists. It has never been used with FR turbulence or LTT, and there is no case of it in S1–S4. |
| B13 | Combined TCNEQ + RANS-LTT: not existing | **Confirmed** | [V] no turbulence in FluxReconstructionNEQ; NEQKOmega only as variable sets | |

---

## 3. Broken or suspicious items

### Convergence
1. **No RANS-LTT run in S1 is converged.**
   - Every run stopped at MaxIter, stopped by hand, or stopped after a few iterations with the sentinel residual −1.79769e+308.
   - The best log10 residual drop is about 4 orders (Kimmel P2 explicit: 1.05 → −0.65).
   - Every `error.txt` contains "terminate called without an active exception / Aborted". For runs whose log ends with "COOLFluiD Environment Terminated" (LM P0/P1/P3, GA P3 explicit, P3/P4 implicit, DiCristina P2), this is a teardown crash and not divergence [I].
   - Other logs end in mid-iteration (GA P1, P1bis, P2, P2bis, P4 explicit, DiCristina P1, LM P2). Those processes were killed or aborted while running, and the cause cannot be determined from the files [V].
2. **The implicit γ-α runs perform worse than the explicit ones** [V]. For example, P2-implTry recorded 14 negative k and 16 negative ω, and the P4 implicit run stalled. Implicit robustness for hypersonic RANS-LTT is the main open numerical problem.
3. **The Langtry-Menter Kimmel runs failed almost immediately** at P0 implicit (2 it), P2 (50 it) and P3 (41 it). One `error.txt` notes edits that were tried (`IsAxiSym` false, `InitComds`).

### Case setup
4. **Inconsistent free-stream Tu across "the same" case.** Kimmel uses 0.475 / 0.50 / 0.75 % and DiCristina 1.0 % (P1) / 0.75 % (P2). Restarts mix Tu values (file names `…Tu075`, `…Tu0475`).
5. **Missing restarts:**
   - `kimmelGAKLogOCoarseP1-iter_8930000.CFmesh`, needed by `Usable_Cases/kimmel2D/P2/kimmelFRP2_impl`, exists nowhere.
   - Several `Usable_Cases` restarts exist only in `Cases/`. They were fixed in the copies.
   - The γ-α flat-plate repo meshes (`GAT3AP0New`, `GA600T3AP1LamNew`, `GReKLogOFlatPlate143.CFmesh`) are missing.
   - LM T3A `T3AP4New` and P0 inputs are missing.
6. **Misleading names:**
   - `Usable_Cases/kimmel2D/P1/kimmelFRP0.CFcase` is a P1 case.
   - `kimmelFRP2_expl.CFcase` uses NewtonIterator.
   - Several `Cases/*-impl` / `*bis` folders contain FwdEuler cases.
7. **DiCristina P1:** the mesh-upgrade subsystem reads the 1128-cell mesh, but the FR solver reads Ray's 1540-cell P1 restart. The two meshes differ, and the upgrade output is unused.
8. **Kimmel P5:** `MaxIter = 1`, so P5 was never run.
9. **`kimmelFRP2_implOB`** uses the generic NS `ConvRHSJacob`/`DiffRHSJacob` instead of the `*Turb` commands. Check that the turbulent diffusion terms are actually treated in that path [I].

### Code
10. **Uncommitted code:**
    - All 2025 γ-α work, Karsch spot growth, onset tables, PhysicalityGReKO, FRTurbProbes and the LM2009 AD Jacobian are uncommitted or untracked on `fr-turbulence`. That is a high risk of loss.
    - There is a parallel external repo (CF-FR-Turbulence) with a merge kit prepared today.
    - Three trees (main, VIN_26, RAY) have diverged FluxReconstructionTurb versions (46 differences each).
11. **3D coverage gaps in FR turbulence code:**
    - no `SA3DSourceTerm` in FR;
    - no `KOmega3DSourceTerm` in FR;
    - no 3D subsonic turbulent BCs;
    - `GammaAlpha3DSourceTerm` has none of the 2025 options.
12. **Wall distance in 3D** is a planar projection to the nearest face centre: approximate on curved Q2 walls.
13. **Stale binaries:** Salvatore's solver symlinks point to `/home/r1132394/COOLFluid` and are broken. The build his results came from is not on this machine.

### Second pass (WORK-HD)

14. **Stored CFcases differ from the configs Ray ran.** Rerunning a stored CFcase does not reproduce Ray; use the `AS_RUN_on_VSC_*.txt` echoes saved in the case folders. Example (Kimmel 3D P1):
    - the stored file has a 3240-hex mesh and MaxIter 1;
    - the run used a 19,440-hex "AoA1" restart, MaxIter 1e10 and `transitionP1.inter`.
15. **Lost runtime controls:**
    - 15 of Ray's interactive files are missing. Placeholders were created in the copies so the cases can start.
    - Most production restarts are missing: HIFiRE P1–P3, DiCristina 3D, Kimmel 3D P2/P3 and Kimmel 2D P0.
16. **Possible Tu inconsistency in HIFiRE:** the P4 restart is labelled `Tu04` (0.4 %). The γ-α onset depends strongly on it; check it against the flight estimate (relevant to C4/A3) [I].
17. **HEXAFLY-INT was only set up** (MaxIter 1, CFL 1e-5), on a Q1 mesh. There is no evidence that it ran [I].
18. **`kimmelFR3DP0` (γ-α variable set) has the γ-α source term commented out.** It is an initialisation step, not a transition case.
19. **No experimental or flight reference data** (Stanton, transition onset) was found in S3/S4. Comparison data for C1/C2 has to come from the literature or from Ray's thesis.
20. **Ray's convergence:** the 3D explicit runs needed 0.2–3 M iterations on 72–288 ranks, chained in 24–72 h jobs, and never converged by residual. This strongly motivates the implicit and matrix-free work (TODO D2).
