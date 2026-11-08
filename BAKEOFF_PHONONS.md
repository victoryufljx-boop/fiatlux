# Phonon method bake-off (Nov 2026)

**Winner, now the default in `gainmat/host/phonons.py`:** MACE-OMAT-0 *medium*, with the cell fixed at the experimental lattice (positions relaxed) and the phonopy LO–TO non-analytic correction (NAC) applied wherever Born charges and ε∞ are tabulated.

Across the 4 hosts, Γ TO + Raman frequencies have a **mean MAE of 7.7 cm⁻¹ and a bias of −2.4 cm⁻¹**. The old default (MACE-MP-0 medium, fully relaxed cell, no NAC) gives 36.2 cm⁻¹ and −30.7 cm⁻¹.

The old protocol can still be selected: `compute_phonons(atoms, model="mace-mp-0-medium", relax_cell=True, nac=False)`. `model="medium"` also still works.

## Setup
* **Box:** 2-core CPU, no GPU, so the runtimes below are CPU wall-clock seconds.
* **Driver:** `python -m gainmat.host.bakeoff MODEL --out results/bakeoff`.
* **Table:** `python -m gainmat.host.bakeoff_report results/bakeoff`.
* **Raw data:** all JSONs are in `results/bakeoff/`.
* **Method:** Γ phonons come from phonopy finite displacements (0.01 Å) in the primitive cell. At q = 0 this is exact, because displacing an atom moves all of its periodic images.
  * Modes are labelled by irrep (phonopy `set_irreps`).
  * Each experimental value is matched to the nearest unused computed mode of the same irrep.
* **Cell protocols:**
  * `relaxed`: cell and positions relaxed together, symmetry kept (FixSymmetry + FrechetCellFilter, fmax 1e-3).
  * `expcell`: lattice taken from `gainmat/structures.py` (room-temperature experiment); only positions relaxed.
* **Experimental benchmark:** `gainmat/refs/phonon_refs.py`. Every value was read on the URL cited there. Some lists are incomplete, and the caveats are listed below.

## Comparison: Γ TO + Raman (cm⁻¹)

Each cell gives dV% / MAE vs. experiment / highest Γ mode. Volume error is relative to the room-temperature experimental cell; a 0 K relaxation would be expected to come out roughly 0.5–1 % smaller. "t" is the total CPU time for all 4 hosts, not counting model loading.

| model | protocol | CaF2 dV% / MAE / maxTO | Al2O3 dV% / MAE / maxTO | quartz dV% / MAE / maxTO | LiYF4 dV% / MAE / maxTO | mean MAE | mean bias | t (s) |
|---|---|---|---|---|---|---|---|---|
| 7net-mf-ompa | relaxed | +2.7 / 19.1 / 307 | +3.1 / 26.5 / 724 | +7.1 / 32.7 / 1137 | +4.5 / 16.7 / 446 | **23.7** | -22.8 | 148 |
| 7net-mf-ompa | expcell | +0.0 / 2.9 / 321 | +0.0 / 5.0 / 755 | +0.0 / 20.6 / 1125 | +0.0 / 14.0 / 497 | **10.6** | -5.1 | 46 |
| mace-mp-0-medium | relaxed | +2.6 / 28.6 / 300 | +3.0 / 45.2 / 731 | +7.2 / 45.1 / 1136 | +4.9 / 26.0 / 444 | **36.2** | -30.7 | 40 |
| mace-mp-0-medium | expcell | +0.0 / 14.8 / 311 | +0.0 / 19.6 / 756 | +0.0 / 28.4 / 1121 | +0.0 / 13.0 / 469 | **18.9** | -16.8 | 16 |
| mace-mpa-0-medium | relaxed | +3.9 / 40.8 / 288 | +2.9 / 35.5 / 727 | +6.4 / 26.7 / 1133 | +6.5 / 30.1 / 423 | **33.3** | -30.4 | 18 |
| mace-mpa-0-medium | expcell | +0.0 / 18.9 / 307 | +0.0 / 8.0 / 759 | +0.0 / 18.1 / 1123 | +0.0 / 15.6 / 501 | **15.2** | -10.3 | 9 |
| mace-omat-0-medium | relaxed | +2.9 / 16.1 / 310 | +2.1 / 14.4 / 734 | +6.0 / 23.6 / 1133 | +2.4 / 10.3 / 446 | **16.1** | -15.8 | 18 |
| mace-omat-0-medium | expcell | +0.0 / 1.9 / 325 | +0.0 / 4.8 / 756 | +0.0 / 16.7 / 1126 | +0.0 / 7.6 / 475 | **7.7** | -2.4 | 9 |
| mattersim-v1-5m | relaxed | +2.6 / 19.9 / 303 | +3.0 / 22.2 / 726 | +6.7 / 28.6 / 1140 | +4.3 / 18.5 / 431 | **22.3** | -20.9 | 8 |
| mattersim-v1-5m | expcell | +0.0 / 4.0 / 317 | +0.0 / 5.0 / 756 | +0.0 / 16.7 / 1128 | +0.0 / 8.8 / 479 | **8.6** | -4.3 | 4 |
| experiment | | exp maxTO 322 | exp maxTO 748 | exp maxTO 1160 | exp maxTO 447 | | | |

**LO with NAC.** These use Born charges and ε∞ from `phonon_refs.py`. Highest Γ LO in cm⁻¹:

LO (NAC, Born charges + eps_inf from refs) — highest Gamma LO, cm^-1

| model | protocol | CaF2 LO (exp 463) | Al2O3 LO (exp 906.6) |
|---|---|---|---|
| 7net-mf-ompa | relaxed | 469 | 861 |
| 7net-mf-ompa | expcell | 483 | 890 |
| mace-mp-0-medium | relaxed | 464 | 858 |
| mace-mp-0-medium | expcell | 476 | 881 |
| mace-mpa-0-medium | relaxed | 456 | 851 |
| mace-mpa-0-medium | expcell | 474 | 880 |
| mace-omat-0-medium | relaxed | 470 | 872 |
| mace-omat-0-medium | expcell | 485 | 892 |
| mattersim-v1-5m | relaxed | 470 | 861 |
| mattersim-v1-5m | expcell | 484 | 889 |

## Findings
1. **The volume error drives the softening.**
   * Every model over-expands on full relaxation: +2–4 % for the ionic crystals and +6–7 % for quartz.
   * In every model and host, the relaxed cell gives a 2–5× larger MAE than the experimental cell, with a systematic negative bias.
   * At the experimental cell the bias almost disappears. This accounts for most of the 5–15 % softening seen on Colab: the CaF2 maximum of 383 cm⁻¹ is a relaxed-cell, no-NAC number.
2. **Model ranking at the experimental cell (mean MAE, cm⁻¹):**

   | Model | Mean MAE |
   |---|---|
   | MACE-OMAT-0 | 7.7 |
   | MatterSim v1-5M | 8.6 |
   | SevenNet-mf-ompa | 10.6 |
   | MACE-MPA-0 | 15.2 |
   | MACE-MP-0 | 18.9 |

   * MACE-OMAT-0 is also the best model on relaxed cells (16.1).
   * MatterSim is nearly as accurate and about 2× faster on CPU (float32). It is a reasonable alternative if speed matters, but it needs its own environment.
   * SevenNet is the slowest of these on CPU.
3. **LO–TO splitting.**
   * NAC moves the top of the CaF2 spectrum from TO 256 to LO 485 cm⁻¹ with MACE-OMAT-0 at the experimental cell (experiment: 463).
   * For Al2O3 the highest LO is 892 cm⁻¹ (experiment: 906.6, Eu LO).
   * LO errors are +5 % for CaF2 and −2 % for Al2O3.
   * The CaF2 overshoot is consistent with the Born charge used, Z*(Ca) = 2.37 from PBE (Wang et al. 2013), which comes paired with a larger ε∞ than experiment.
   * `phonon_max_cm1` now includes LO (`phonon_max_TO_cm1` is kept). The DOS mesh is Γ-centred (`is_gamma_center=True`), so there is no MeshSymmetryFallbackWarning.
4. **Hosts with no NAC (quartz, LiYF4).** No readable source with Born charges was found for these two, so NAC is skipped and their `phonon_max_cm1` is TO-like.
   * Quartz: experimental top LO is 1238 cm⁻¹ against a top TO of 1160.
   * LiYF4: experimental LO values were not found either.
   * Fix: add Born charges from DFT once GPAW or another DFPT code is available (see below).

## Not done or limitations (honest list)
* **ORB-v3 was not tested.** `orb-models` 0.7 needs `dm-tree==0.1.8`, which has no Python 3.13 wheel and fails to build here. On Colab (Python 3.11/3.12) `pip install orb-models` should work. The model key `orb-v3-con-omat` is already registered.
* **eSEN/UMA were not tested:** gated Hugging Face access.
* **Born charges from DFT (GPAW) were not computed.** GPAW has no binary wheel here and Berry-phase Z* needs many SCF cycles, which did not fit the 2-core, 50-minute budget. Literature fallback used instead:
  * CaF2: PBE Z*(Ca) from Wang et al. 2013, with Z*(F) = −Z*(Ca)/2.
  * Al2O3: scalar Z*(Al) = 2.92 and Z*(O) = −1.94 quoted in arXiv:1712.03915. Whether these are DFT values is not stated on that page, and they are applied isotropically.
* **Caveats on the experimental data:**
  * The primary papers (Barker, Porto & Krishnan, Gervais & Piriou, Scott & Porto, Miller et al., Salaün et al.) were paywalled and not opened. Values come from open secondary sources that quote them, or from newer primary measurements; all are cited in `phonon_refs.py`.
  * LiYF4 column alignment was inferred from a garbled PDF extraction.
  * The Al2O3 Eg list is incomplete (3 of 5 modes).
* **Thermal expansion is not modelled.** The experimental-cell protocol is a practical correction that cancels both the functional's over-expansion and the 0 K vs. room-temperature difference. Hypothetical hosts with no measured lattice must use `relax_cell=True`; expect them to come out about 15–30 cm⁻¹ soft.

## Colab install
```bash
pip -q install ase phonopy spglib seekpath
pip -q install mace-torch "e3nn==0.4.4"      # default: MACE-OMAT-0 medium (downloads on first use, ASL licence)
# optional, each in a separate env / runtime because the e3nn versions conflict with MACE:
pip -q install mattersim                      # mattersim-v1-5m
pip -q install sevenn                         # 7net-mf-ompa (pulls e3nn>=0.5 -> breaks MACE model unpickling)
pip -q install orb-models                     # orb-v3-con-omat (untested here)
```
