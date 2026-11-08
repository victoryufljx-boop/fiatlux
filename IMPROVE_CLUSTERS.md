# Embedded-cluster tier: improvements (branch `improve/clusters`)

Goal: move `gainmat/cluster` towards standard AIMP-embedded cluster CASSCF/CASPT2 practice
(Seijo & Barandiarán), where Ce3+ 4f→5d errors are typically 1000–2000 cm-1.
All numbers below were computed in this session on **2 CPU cores** (1 thread per job, two jobs
in parallel) at the cheap `fast` / `fast+` levels. Nothing at `full` was run here.

## What changed

| Item | Status | Where |
|---|---|---|
| (a) Repulsive embedding of the next cation shell. No full-core Y ECP exists in PySCF or basis_set_exchange (every Y ECP there leaves 4s4p in valence), so a **frozen-ion AIMP-like potential** was implemented. Each Y3+ ion adds −Z_eff/r + its semilocal def2 ECP28 + exact 2J − K of the frozen free-ion Y3+ RHF density + a Huzinaga/AIMP core-projection shift Σ_c(−2ε_c)\|c⟩⟨c\| to the cluster h_core. There are **no basis functions on Y**. All 8 nearest Y3+ are included (`r_cap=5.3`). `cap_mode='pc'\|'partial'\|'aimp'`; the old `partial_cap` still works as `cap_mode='partial'`. | **implemented + run** (Ce `fast+`) | `embedding.py` (`aimp_hcore`, `attach_aimp`) |
| (b) Larger basis. Ce def2-SVP/def2-TZVP (ECP28MWB, 4s4p4d5s5p4f5d6s valence, with g functions) are bundled as files from basis_set_exchange, so Colab needs no network. F gets def2-TZVPD and O gets def2-SVPD/TZVPD (polarisation + diffuse). `'bundled:<file>'` basis keys. | **coded**; Ce/F big basis **not run** here (too expensive on 2 cores). Ti with O def2-SVPD run at `fast+`. | `basis_data/`, `embedding.resolve_basis`, `ce_5d.LEVELS['full']` |
| (c) Semicore 5s5p correlation in NEVPT2. PySCF NEVPT2 has no frozen-core option, so the frozen orbitals are folded into an effective h_core (h + 2J − K, plus their energy as a constant) and NEVPT2 is rerun. The baseline already correlated everything, Ce 5s5p included. | **implemented**; result below | `casci_lf.nevpt2_frozen`, `semicore_indices` |
| (d) Relaxation of the first anion shell inside the frozen embedding. Symmetry-adapted coordinates: breathing (1 coordinate, or 2 for the S4 short/long Ce–F sets), and breathing + e_g(θ) for TiO6. Ground state from UHF; lowest 5d / 2Eg from MOM ΔSCF (which breaks the e_g degeneracy, i.e. finds the static JT minimum). Quadratic fits give Q_gs, Q_ex and both relaxation energies; Stokes = λ_ex + λ_gs. Ce absorption is NEVPT2 at Q_gs, and emission is vertical NEVPT2 at Q_gs minus the SCF Stokes shift. Ti reruns NEVPT2 at the 2Eg minimum. | Ce: **implemented + run**. Ti: **coded, not converged here** (see below) | `relax.py`, `ce_5d.run(relax=...)`, `ti_3d.run(relax='jt')` |
| Notebook | `notebooks/build_03.py` (new) regenerates `03_embedded_clusters.ipynb`. `LEVEL='full'` uses everything above. Flags `CAP_MODE`, `R_CAP`, `CE_RELAX`, `TI_RELAX`, `SEMICORE` switch each piece off. | done |
| Tests | AIMP layout and repulsion, bundled basis parse, relaxation modes and quadratic fit (fast); slow Ti JT test | `tests/test_cluster.py` (11 fast tests pass) |

## Before / after (cm-1)

Before = original `fast` (SBKJC Ce / def2-SVP F, 4 Y as bare point charges, vertical only).
After = `fast+` (same basis, all 8 Y as AIMP ions, breathing-relaxed F shell, SCF Stokes shift).

### Ce3+:LiYF4

| quantity | before | AIMP only (host geometry) | after: AIMP + relaxed (`fast+`) | experiment |
|---|---|---|---|---|
| 4f→5d1 | 36,406 | 36,886 | 36,052 | 33,433 |
| 4f→5d2 | 46,527 | 47,534 | 47,401 | 41,101 |
| 4f→5d3 | 57,630 | 58,902 | 59,708 | 48,564 |
| 4f→5d4 | 58,474 | 59,751 | 60,572 | 50,499 |
| 4f→5d5 | 63,721 | 65,110 | 66,217 | 52,790 |
| 5d spread (5d5−5d1) | 27,315 | 28,224 | 30,165 | 19,357 |
| emission 5d1→2F5/2 | 36,406 (vertical) | – | 35,580 (relaxed) | 32,680 (306 nm) |
| emission 5d1→2F7/2 | 34,206 (vertical) | – | 33,379 (relaxed) | 30,770 (325 nm) |
| Stokes shift | – | – | 473 | ~750 (derived) |

The "before" Ce row is the given baseline (Ti baseline re-run here: 24,258 cm-1, matching).
The baseline Ce rerun was stopped to free cores.

Relaxation (UHF / MOM ΔSCF, 1-coordinate breathing): the ground-state minimum is at
ΔR = −0.034 Å relative to the host Y–F distances, and the lowest 5d minimum at −0.054 Å.
λ_ex = 234 cm-1, λ_gs = 239 cm-1. The 5d state contracting relative to 4f agrees with Seijo &
Barandiarán. The ground-state *contraction* is unphysical, since Ce3+ is larger than Y3+;
experiment and EXAFS show an expansion of ~+0.03–0.05 Å. It points to the small basis /
UHF level and to the AIMP Y3+ potential being too soft or too repulsive. The 5d minimum also
lies just outside the scanned grid (0, 0.05, 0.10 Å), so it is extrapolated.

Semicore (c) (AIMP, host geometry; 4 orbitals = Ce 5s + 5p frozen vs correlated in NEVPT2).
5d levels with 5s5p frozen: 35,923 / 46,968 / 58,700 / 59,558 / 65,198. Correlated minus frozen:
**+963 / +566 / +202 / +193 / −89 cm-1**. Correlating 5s5p raises 5d1 by ~1,000 cm-1 and
compresses the 5d spread by ~1,050 cm-1. The effect is at the level of the target accuracy, so it
must stay correlated (the default; the baseline already did this). It does not explain the
~8,000 cm-1 over-splitting.

### Ti3+:Al2O3

| quantity | before | after | experiment |
|---|---|---|---|
| 2T2g→2Eg (mean) | 24,258 | 23,938 (O def2-SVPD, vertical) | ~20,400 (490 nm) |
| emission | – | not obtained | ~12,800 (780 nm) |
| Stokes shift | – | not obtained | ~7,600 (derived) |

The Ti JT scan (3×3 grid in breathing × e_θ) is coded, but MOM ΔSCF with O def2-SVPD collapsed
or failed to converge at the e_θ = 0 (trigonal) points. A retry with damping is built in;
failed points are dropped from the fit. Each failed point cost ~4 min on 1 core, so the scan was
stopped. Suggested fix, not done: follow the excited state along the grid from the previous
converged excited orbitals (MOM reference = neighbour), and start from e_θ ≠ 0.

## Interpretation

* The AIMP Y3+ embedding together with the (contracted) relaxed geometry did **not** reduce
  the 5d over-splitting. It widened it by ~10% (spread 27.3k → 30.2k vs 19.4k exp). With
  frozen Y3+ ions in place of bare +3 charges, the 5d leakage onto Y is not the main cause of
  the ~40% too-large spread at this level. The more likely causes are the minimal SBKJC Ce basis
  (5d radial flexibility) and the ligand basis (no diffuse F, which (b) addresses at `full`),
  plus the spurious ground-state contraction, which raises the crystal field. The lowest 5d and
  the emission moved ~400–800 cm-1 *towards* experiment.
* The NEVPT2 5d1 is still ~2,600 cm-1 too high. The Stokes shift is the right order of
  magnitude (473 vs ~750 cm-1).
* The honest status: (a), (c) and (d-Ce) are implemented and demonstrated, and (b) is wired
  into `full` but untested at that size. The decisive test, a `full` Colab run with
  Ce def2-TZVP / F def2-TZVPD, has not been done.

## Expected Colab runtimes (12-core high-RAM CPU runtime, estimates)

Measured here on 1 thread: Ti `fast` 49 s (2 threads). Ce `fast+` 848 s total: relaxation scan
286 s (3 points × UHF + MOM, with AIMP J/K ≈ 20 s per SCF), ROHF 44 s, SA-CASSCF 49 s,
NEVPT2 233 s, semicore NEVPT2 ≈ 230 s.

| run | estimate |
|---|---|
| Ti `fast+` (vertical + JT 9-point scan + 2 NEVPT2) | 5–15 min |
| Ce `fast+` | 5–10 min |
| Ti `full` (def2-TZVP/TZVPD, X2C, 20 Å, JT scan) | 1–2 h |
| Ce `full` (def2-TZVP Ce ~110 AO + 8× def2-TZVPD F ~ 400 AO total, 9-point 2-coordinate scan, SA-CASSCF(1,12), 2× NEVPT2 over 12 roots) | 4–8 h (NEVPT2 dominates; set `SEMICORE=False` and `CE_RELAX='breathing'` to roughly halve it) |
| Ce `full` with `CE_FULL_BASIS='bundled:ce_def2-svp.nw'`, `F_FULL_BASIS='def2-svpd'` | 1–2 h |

## Files
`gainmat/cluster/{embedding,casci_lf,relax,ce_5d,ti_3d,run_cluster}.py`, `gainmat/cluster/basis_data/*.nw`,
`notebooks/build_03.py`, `notebooks/03_embedded_clusters.ipynb`, `tests/test_cluster.py`, raw logs/JSON in `runs/`.
