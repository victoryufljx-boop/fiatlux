# Band-gap and refractive-index bake-off (GPAW 26.7, PAW, plane waves)

**Winners**

| Quantity | Winner | Notes |
|---|---|---|
| Band gap | **GLLB-SC** | Now the default in `gainmat/host/gap.py`. `method="pbe" \| "tb09" \| "hse06"` remain selectable. |
| ε∞ / n | **RPA with local fields on PBE, no scissor** | Already the reported `eps_inf`. Now also gives per-direction values for uniaxial crystals. |

GLLB-SC is the only method tested that lands near the experimental *fundamental* gap of wide-gap fluorides and oxides. It costs about 1.2–1.3× PBE. TB09 and HSE06@PBE both come in 1–3 eV low, and TB09 is several times more expensive and fails to converge for α-quartz.

## Setup

- **Code:** GPAW 26.7.0 (libxc 5.2.3, FFTW, OpenMPI). PW mode, Fermi–Dirac smearing 0.01 eV.
- **Gaps** are read on Γ-centred meshes, which contain the band edges here:
  - LiF, MgO, CaF2: 6×6×6.
  - Al2O3, quartz: 4×4×4.
  - Pass 1: ecut 400 eV, run serially.
  - Pass 2 (this one): runs under `mpiexec -n 2`.
- **Structures:**
  - LiF a = 4.0263 Å, MgO a = 4.212 Å (experimental).
  - CaF2, Al2O3, quartz from `gainmat/structures.py`.
- **HSE06** is non-self-consistent on PBE wavefunctions, evaluated at the PBE VBM/CBM k-points plus Γ.
- **TB09** uses `{"name": "TB09", "stencil": 1}` with the c parameter self-consistent.
- **Raw data:**
  - `bench/results.jsonl` and the `bench/*.txt` logs.
  - `eps/results_eps.jsonl`.
  - Scripts: `bench/bench_gap.py`, `eps/bench_eps.py`, `bench/run{A,B,C,D}.sh`, `bench/yag/`.

## 1. Experimental references

We compare against the **fundamental (quasiparticle) gap**, not the first exciton peak:

- PBE, GLLB-SC, TB09 and HSE06 are all single-particle or quasiparticle gap methods. None of them includes the electron–hole interaction.
- In these ionic crystals the exciton binding energy is large: ~1.3–1.6 eV in LiF and ~1 eV in CaF2. Comparing with the absorption peak would reward methods that underestimate the gap.
- The exciton peak is still the practical UV edge, so `compare_to_experiment` also reports GLLB-SC against it.

| Material | Fundamental gap (eV) | First exciton peak (eV) | Sources |
|---|---|---|---|
| LiF | 14.2 ± 0.2 (Piacentini, Lynch & Olson, PRB 13, 5530 (1976)) | 12.9 (Fowler 1968) – 13.5 | https://scilit.com/publications/754a15725ef14e29e21924f7096660e7 , https://www.homepages.ucl.ac.uk/~ucapahh/research/crystal/lif.htm |
| MgO | 7.83 (Whited 1973); 7.7 (Roessler & Walker 1967) | ~7.7 (binding energy ≈ 80 meV) | https://ar5iv.labs.arxiv.org/html/2012.08960 |
| CaF2 | 12.1 direct Γ–Γ; 11.8 indirect X–Γ (Rubloff, PRB 5, 662 (1972)) | 11.1–11.2 | https://arxiv.org/pdf/1107.5962 (Table III), https://arxiv.org/pdf/cond-mat/0306510 |
| α-Al2O3 | 8.8 (French 1990, used by benchmarks such as JARVIS); 9.57 (French, Jones & Loughin 1994); 9.6 (Will 1992) | ~9.2 | https://arxiv.org/pdf/2301.04199 , https://arxiv.org/pdf/cond-mat/0610176 , https://arxiv.org/pdf/1804.01032 |
| α-quartz | 8.9 (Chelikowsky & Schlüter 1977; DiStefano & Eastman ~9); 9.65 (Garvie et al. 2000, EELS) | ~10.3 (first absorption peak, as quoted) | https://cyberleninka.ru/article/n/ab-initio-calculations-of-band-structures-and-dielectric-functions-of-quartz , https://scilit.com/publications/ad3eb9d2cd922c053354f678a136f828 |

Al2O3 and quartz each have two experimental values that differ by about 0.7 eV. Both are carried below as a range ("lo / hi").

Refractive index and ε∞:

- n is taken at 589 nm.
- ε∞ (the long-wavelength electronic limit) is taken from the UV terms of the Sellmeier fit on each page.

| Material | n (589 nm) | ε∞ | Source |
|---|---|---|---|
| LiF | 1.392 | 1.93 (1.90–1.93) | https://refractiveindex.info/?shelf=main&book=LiF&page=Li |
| MgO | 1.737 | 2.94–2.95 (Roessler & Walker; Henderson & Hughes) | https://refractiveindex.info/?book=MgO&page=Stephens&shelf=main , https://ar5iv.labs.arxiv.org/html/2012.08960 |
| CaF2 | 1.434 | 2.04 | https://refractiveindex.info/?book=CaF2&page=Li&shelf=main |
| Al2O3 o / e | 1.768 / 1.760 | 3.08 / 3.06 | https://refractiveindex.info/?shelf=3d&book=crystals&page=sapphire , https://refractiveindex.info/?shelf=main&book=Al2O3&page=Malitson-e |
| quartz o / e | 1.544 / 1.553 | 2.36 / 2.38 | https://refractiveindex.info/?book=SiO2&page=Ghosh-o&shelf=main , https://refractiveindex.info/?shelf=main&book=SiO2&page=Ghosh-e |

## 2. Gap comparison (eV)

Errors are computed against the fundamental gap. For Al2O3 and quartz they are given as lo / hi, using the two experimental values in section 1.

| Material | Exp. | PBE | GLLB-SC | TB09 | HSE06@PBE |
|---|---|---|---|---|---|
| LiF | 14.2 | 9.10 (−5.10) | **14.88 (+0.68)** | 11.90 (−2.30) | 11.28 (−2.92) |
| MgO | 7.83 | 4.58 (−3.25) | **8.01 (+0.18)** | 6.77 (−1.06) | 6.29 (−1.54) |
| CaF2 | 12.1 | 7.26 (−4.84) | **12.94 (+0.84)** | 10.32 (−1.78) | 9.35 (−2.75) |
| Al2O3 | 8.8 / 9.57 | 6.09 (−2.7 / −3.5) | **9.51 (+0.71 / −0.06)** | 8.22 (−0.58 / −1.35) | 7.91 (−0.89 / −1.66) |
| α-quartz | 8.9 / 9.65 | 5.87 (−3.0 / −3.8) | **9.81 (+0.91 / +0.16)** | did not converge* | 7.68 (−1.22 / −1.97) |
| **MAE** (5 materials) | | 3.8 / 4.1 | **0.66 / 0.38** | 1.4 / 1.6 (4 materials) | 1.9 / 2.2 |
| Cost vs PBE (same settings) | | 1× | ~1.2–1.3× | ~3–10× (LiF 59 s vs 6 s serial) | ~1.1–1.6× (nsc at 1–3 k-points) |

\* **Quartz TB09:** the SCF diverged twice: with Pulay mixing β = 0.5, and again with β = 0.05 (the energy runs away after about 100 iterations). Logs: `bench/*quartz_TB09*diverged*`. This is consistent with the TB09 potential being ill-behaved in the low-density channels of quartz. It is a further practical argument against TB09 as a default.

### Comparison with the literature

**TB09 is under-converged relative to published TB-mBJ.** Our values are low:

- All-electron WIEN2k values (Tran & Blaha, PRL 102, 226401 (2009), https://scispace.com/pdf/accurate-band-gaps-of-semiconductors-and-insulators-with-a-1nnqe1ie65.pdf):

  | Material | Ours | WIEN2k (TB-mBJ) | Difference |
  |---|---|---|---|
  | LiF | 11.90 | 12.94 | −1.0 |
  | MgO | 6.77 | 7.17 | −0.4 |

- Borlido et al. 2020 (npj Comput. Mater. 6, 96, Table 5, https://arxiv.org/pdf/2003.10177), mBJ:

  | Material | WIEN2k | VASP | Ours |
  |---|---|---|---|
  | LiF | 12.82 | 12.08 | 11.90 |
  | Al2O3 | 8.31 | 7.86 | 8.22 |

  So our values sit with the pseudopotential/PAW results, not the all-electron ones.

- The cause is known: with a pseudo-density, the self-consistent c parameter comes out too small. Borlido 2020 shows the effect grows with the gap. Our c values are LiF 1.571, MgO 1.437, CaF2 1.569 and Al2O3 1.438.

**Convergence is not the problem.** The k-mesh and cutoff are converged (section 3), so the TB09 deficit is a PAW/frozen-core limitation, not under-convergence. Even all-electron TB-mBJ is still about 1.3 eV low for LiF.

**Overall error statistics** (Tran, Ehsan & Blaha 2018, 76 solids, https://ar5iv.arxiv.org/html/1712.06508):

| | PBE | mBJ | HSE06 | GLLB-SC |
|---|---|---|---|---|
| MAE (eV) | 1.99 | 0.47 | 0.82 | 0.64 |
| Mean error (eV) | −1.99 | −0.30 | −0.68 | +0.20 |

That paper also finds GLLB-SC "overestimat[es] band gaps larger than 4 eV", which is exactly our pattern. Borlido 2020 finds mBJ and HSE06 best on average, but over a set dominated by smaller gaps. For the >7 eV ionic insulators relevant to UV gain hosts, our GPAW data favour GLLB-SC.

**Our GLLB-SC matches the published GPAW value.** LiF 14.88 eV, against 14.96 eV in Kuisma et al., PRB 82, 115106 (2010), Table I, https://arxiv.org/pdf/1003.0296.

**HSE06 is low because it is non-self-consistent on PBE.** Our MgO 6.29 compares with 6.67 for self-consistent HSE in Tran & Blaha 2009. Self-consistent HSE would add ~0.3–0.5 eV but cost 10²× more. It would still be ~1–2 eV low for LiF and CaF2.

**Not verified:** the Borlido 2019 JCTC per-material table could not be read (captcha / proxy), so it is not quoted.

## 3. Convergence of GLLB-SC (LiF, CaF2)

Gaps in eV, from `bench/log_conv_*` and `bench/results.jsonl`.

| ecut (eV) / k-mesh | LiF | CaF2 |
|---|---|---|
| 400 / 6³ | 14.883 | 12.942 |
| 500 / 6³ | 14.895 | 12.965 |
| 600 / 6³ | 14.874 | 12.943 |
| 700 / 6³ | 14.871 | 12.943 |
| 500 / 4³ | 14.903 | 12.963 |
| 500 / 8³ | 14.899 | 12.962 |
| 500 / 10³ | 14.896 | 12.963 |

- The gap is converged to **±0.03 eV** at ecut ≥ 400 eV and k ≥ 4³. The discontinuity Δxc changes by < 0.03 eV.
- The `full` level (ecut 500, path band edges, k = 2³ for CaF2) gives 12.98 eV, within 0.04 eV of the converged value (`bench/yag/calib_caf2.jsonl`).
- The default settings (`fast`: ecut 500; `full`: ecut 500 + path) are therefore converged well below the method error (~0.4–0.7 eV).

## 4. Refractive index (ε∞)

- **Method:** `gpaw.response.DielectricFunction` at ω = 0, η = 10⁻⁴, with local-field effects, on a PBE ground state.
- **Bands:** 6–8 × n_occ. LiF changes by < 0.2 % going from 6× to 8×.
- **Response cutoff:** ecut_resp 50–60 eV.
- **k-mesh:** 6–8³ for cubic, 4³ for Al2O3 and quartz.
- **Scissor:** the GLLB-SC − PBE gap, applied via `eshift`.

| Material | ε∞ exp. | RPA no LFE | **RPA + LFE** | RPA + LFE + scissor | n (RPA + LFE) | n exp. (589 nm) |
|---|---|---|---|---|---|---|
| LiF | 1.93 | 2.03 | **1.943 (+0.9 %)** | 1.704 (−11.5 %) | 1.394 | 1.392 |
| MgO | 2.94–2.95 | 3.28 | **3.122 (+6 %)** | 2.652 (−10 %) | 1.767 | 1.737 |
| CaF2 | 2.04 | 2.34 | **2.192 (+7.5 %)** | 1.868 (−8.4 %) | 1.481 | 1.434 |
| Al2O3 ⊥c (o) | 3.08 | 3.32 | **3.162 (+2.6 %)** | 2.723 (−11.6 %) | 1.778 | 1.768 |
| Al2O3 ∥c (e) | 3.06 | 3.27 | **3.122 (+2.2 %)** | 2.694 (−11.8 %) | 1.767 | 1.760 |
| quartz ⊥c (o) | 2.36 | 2.54 | **2.363 (+0.3 %)** | 2.067 (−12.3 %) | 1.537 | 1.544 |
| quartz ∥c (e) | 2.38 | 2.54 | **2.390 (+0.3 %)** | 2.086 (−12.5 %) | 1.546 | 1.553 |

**Unshifted RPA + LFE @PBE wins:**

- Mean |error| is 3 % in ε∞ and ~1 % in n (against n at 589 nm).
- Adding the scissor makes every material 8–12 % too low. A pure scissor without renormalising the momentum matrix elements removes the PBE error cancellation: the gap is too small, but there are no excitonic effects.
- `compute_eps_inf` therefore reports the unshifted value as `eps_inf` and keeps the scissored one only for reference.
- For non-cubic cells it now averages x, y and z, and returns per-direction values (`eps_inf_rpa_dir`, `n_inf_dir`).

**DFPT / finite field:** GPAW 26.7 has no DFPT dielectric tensor. Its `ConstantElectricField` only works along non-periodic directions; the PW mode is periodic in all three. `gpaw.berryphase` and `gpaw.borncharges` give polarisation and Born charges from ionic displacements, but there is no periodic finite-field electronic response. So a DFPT or finite-field ε∞ is not available in GPAW and was not run. RPA + LFE is the supported route.

**Caveats:**

- The Al2O3 and quartz meshes (4³) are modest.
- MgO and CaF2 come out 6–8 % high, as is typical for RPA@PBE in ionic crystals.

## 5. MPI launcher and YAG runtime

**Problem:** in a Colab kernel GPAW runs as one serial rank.

**New module `gainmat/host/mpi_run.py`:**

- Writes the atoms plus a small per-host driver script.
- Runs `mpiexec -n <ncpu> --allow-run-as-root --oversubscribe python run_<host>.py`.
- Collects the rank-0 JSON.
- **Gotcha found:** GPAW ≥ 25 refuses to run under a launcher unless `GPAW_MPI_BACKEND=cgpaw` is set or it is started as `gpaw python`. The launcher sets this.

**Wiring:**

- `run_host(..., gpaw_nproc=N)` uses the launcher.
- `notebooks/01_host_properties.ipynb`, rebuilt by `gainmat/host/build_notebook.py`, passes `gpaw_nproc=os.cpu_count()`. Its install cell now installs OpenMPI and checks the rank count.
- Tests: `tests/test_host.py::test_caf2_gap_mpi_smoke` and `::test_mpi_command`.

**Measured on this 2-vCPU box** (`bench/yag/`):

- **CaF2 `full` level (PBE + GLLB-SC, ecut 500, path):** 20.4 s serial and 14.1 s with 2 ranks through the launcher, identical gaps (12.98 eV). The Colab run took ~270 s for this. That suggests the in-kernel run was also suffering from BLAS-thread oversubscription (`OPENBLAS_NUM_THREADS = ncpu` in a single rank), not just running serially. The MPI path pins 1 thread per rank.
- **YAG (80 atoms, 292 bands, 2³ Γ-centred mesh = 4 IBZ k-points, ecut 500), PBE:**
  - LCAO initialisation: 36 s.
  - **≈ 27 s per SCF iteration on 2 ranks**, i.e. ≈ 54 core·s per iteration.

**YAG estimate.** Assumptions:

- About 20 SCF iterations each for PBE and GLLB-SC, with GLLB-SC costing ~1.2× PBE.
- Each run's SCF ≈ 20–25 core·min, so ≈ 45 core·min for both.
- The 60-point non-SCF band path costs ~80 core·min per functional.

| Setting | Serial in kernel | 12-vCPU Colab (≈ 6 physical cores, k ∥ 4 × domain) |
|---|---|---|
| Mesh band edges (now the default for YAG: gap direct at Γ) | ~45–60 min | **~10–20 min** |
| 60-point band path (old `full` behaviour) | ~3–5 h | ~40–60 min |

These are extrapolations from measured per-iteration costs, not a completed YAG run. Colab per-core speed is uncertain, so the range is wide. `BAND_EDGES_OVERRIDE = {"Y3Al5O12": "mesh"}` in `run_host.py` applies the cheaper setting.

## 6. Code changes

- **`gainmat/host/gap.py`:**
  - `compute_gap(..., method="gllbsc")`; `method` may be `"pbe" | "tb09" | "hse06"` or a list.
  - Returns `gap_eV` and `gap_method` alongside the per-method keys.
  - `compute_eps_inf` now returns per-direction values, defaults to 6 × n_occ bands, and documents the no-scissor finding.
- **`gainmat/host/run_host.py`:**
  - Every level sets `method="gllbsc"`; `fast` now uses ecut 500.
  - New arguments `gap_method=` and `gpaw_nproc=`.
  - Adds the YAG band-edges override.
  - The UV edge now comes from the primary gap.
- **`gainmat/schema.py`:** new fields `HostProps.gap_eV` and `gap_method`.
- **`gainmat/host/mpi_run.py`:** new (the MPI launcher).
- **`tests/test_host.py`:** new tests for the defaults, the MPI command, a CaF2 MPI smoke run, and a LiF ε∞ smoke run.
- **`SPEC.md`, `README.md`:** updated.
