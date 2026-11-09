# Glass forward model: dopant + glass -> fibre gain windows

Code: `gainmat/forward/`. Driver: `python scripts/run_forward.py`. Tests: `tests/test_forward.py`.

## What it computes

For a rare-earth ion in a glass family, the model predicts absorption and emission cross
sections, level lifetimes and quantum efficiencies, and the **gain-window envelope**: the
wavelength ranges where small-signal gain minus background loss exceeds a threshold
(default 0.1 dB/m) for some pump intensity up to ~1 W in a single-mode core.

## Pipeline

1. **Ion** (`ions.py`). Complete Carnall free-ion Hamiltonian (F^k, zeta, alpha, beta, gamma,
   T^i, M^k, P^k; LaF3 parameters) diagonalised once; free-ion eigenvectors below
   32,000 cm^-1 form a truncated basis carrying the crystal-field, U^(lambda) and L+2S operators.
   Ions: Pr, Nd, Dy, Ho, Er, Tm, Yb.
2. **Sites** (`spectra.random_sites` or MD). Each dopant site is a set of ligand positions.
   Stochastic sites draw coordination number, bond lengths and angular distortion from the
   glass family's first-shell statistics; MD-derived sites (`glass.load_md_sites`) can replace
   them.
3. **Crystal field per site** (`spectra.spm_bkq`). Newman superposition model,
   B^k_q = sum_L Bbar_k (R0/R_L)^t_k conj(C^k_q(theta_L, phi_L)).
4. **Stark-resolved intensities.** Every pair of Stark states gets
   S_ED = e^2 sum_lambda Omega_lambda sum_q |<a|U^lambda_q|b>|^2 and S_MD from L+2S. Summed over
   a multiplet pair this reproduces the Judd-Ofelt line strength exactly, so integrated cross
   sections are JO-normalised and independent of the crystal field (a test checks this).
5. **Spectra.** Einstein / Fuechtbauer-Ladenburg per state pair, Boltzmann populations inside
   each multiplet (McCumber-consistent by construction), homogeneous Gaussian of FWHM
   gamma_h, averaged over sites. The inhomogeneous width is not an input: it comes from site
   disorder.
6. **Rates** (`gain.py`). JO radiative rates between all multiplets; multiphonon decay to the
   next lower multiplet from the energy-gap law (family parameters, or a host-agnostic law
   parameterised by the max phonon energy). The gap is the **Stark-edge gap** of each site
   (lowest Stark level of the upper multiplet minus highest of the lower, the definition used
   by Layne et al. 1977), and the rate is averaged over sites, <B exp(-alpha dE_site)>.
   Gap-law parameters vs 26 room-temperature measured rates: rms error 0.40 decades for
   Layne's glasses (`scripts/check_gap_law.py`, `forward_results/gap_law_check.json`).
   Silica, aluminosilicate and germanosilicate fibre glasses have their own laws, calibrated
   on Tm 3H4 lifetimes at <= 200 ppm (Faure et al. 2007): fused silica decays ~1.7x faster than
   Layne's multicomponent silicate, 6 mol% Al2O3 ~2.6x slower than pure silica.
   Pumping acts on every multiplet pair (ground- and excited-state absorption, stimulated
   emission), and a second pump can be added: Er 985 + 1973 nm dual-wavelength pumping
   (Henderson-Sapir et al. 2014) opens the 3.3-3.75 um band in Er:ZBLAN and Er:fluoroindate.
7. **Steady state + gain.** Pump from the ground multiplet into whichever multiplet absorbs at
   the pump wavelength; solve the multi-level rate equations; the small-signal gain
   g = sum [N_i sigma_em(i->j) - N_j sigma_abs(j->i)] includes excited-state absorption.
8. **Loss.** Fibre background floor plus a multiphonon IR edge anchored on silica fibre
   (7.81e11 exp(-48.48/lambda[um]) dB/km) and scaled by the glass's max phonon energy.
9. **Window** = envelope over pump intensities (and pump wavelengths) of gain - loss >= threshold.

## Inputs that are not ab initio (and where they come from)

| Input | Source |
|---|---|
| Free-ion parameters | Carnall 1989 (LaF3), fitted to experiment but host-transferable |
| Judd-Ofelt Omega per (ion, glass) | measured values (`refs/glass_refs.py`); borrowed from the nearest family when missing, flagged in every result |
| SPM intrinsic parameters | literature (`refs/site_refs.py`) |
| Site statistics (CN, R, disorder) | EXAFS/MD literature, or our own MLIP MD (`md_sites/`) |
| Energy-gap-law parameters | literature per family |
| Max phonon energy, n, loss floor | literature per family |

## Known limitations

- No concentration effects: cross-relaxation, upconversion, clustering. This over-predicts
  gain from levels that self-quench at practical doping (notably Pr 1D2, and Er at high doping),
  and misses beneficial cross-relaxation (Tm 3H4 -> 2x 3F4).
- No upconversion or ESA pumping schemes (e.g. Tm S-band pumped at 1060 nm), no Yb->Er or
  Yb->Pr sensitisation.
- One homogeneous width for all transitions; no phonon sidebands.
- Window thresholds are absolute dB/m at a fixed typical doping per ion.
- The energy-gap law is a family-level fit; specific transitions (e.g. Pr 1G4 in ZBLAN) can be
  off by an order of magnitude in W_nr.

## Results (first full run, Nov 2026)

See `forward_results/SUMMARY.md` (validation and grid tables), `forward_results/coverage.png`,
`forward_results/phonon_sweep.png`, and `forward_results/md_vs_stochastic.json`.

- Band positions of known fibre gain media reproduced (Er C+L, Tm 2 um and S band, Pr O band,
  Ho 2.1/2.9 um, Er 2.8 um, Dy 3 um, Yb 1 um); cross sections within ~25% where measured.
- Crystal-field strength: one calibration per ligand (O: Er 4I15/2 splitting 370 cm^-1 in silica;
  F: Carnall Er:LaF3). Independent checks: Ho 5I8 (silica) 332 vs 330 cm^-1, Yb 2F5/2 645 vs 695,
  Yb 2F7/2 765 vs 1210 (miss -> Yb band ends ~1.10 um instead of ~1.2 um).
- All 10 validation gain bands are reproduced, including Er 3.5 um (dual pumping).
- Largest error: lifetimes of the ~2-3 um upper levels. Stark-edge gaps cut the error from
  10-25x to 2-6x (model / measured, ms: Tm 3F4 silica 2.85 / 0.65 (was 6.28); Ho 5I7 silica
  3.21 / 0.7 (6.70); Pr 1G4 ZBLAN 0.21 / 0.11 (0.99); Dy 6H13/2 ZBLAN 3.87 / 0.65 (15.8)).
  With the Faure-calibrated silica laws (Tm 3H4 in silica: 21 us model vs 14 us measured;
  aluminosilicate 37 vs 37 us) the multiphonon rate of Tm 3F4 and Ho 5I7 is still only
  ~100-500 s^-1, while fibres measure 1300-2700 s^-1, and Tm 3F4 lifetimes fall with
  concentration (0.42 ms at 1.1e20 cm^-3, 0.335 ms at 2.9e20).  What is left is not one global factor (Tm 3F4 needs W x25, Ho 5I7 x5), and fused-silica fibre
  rates sit 10-80x above the multicomponent-silicate gap law that Layne measured. Likely
  causes: fused silica is not Layne's silicate, OH and Tm-Tm / Ho-Ho quenching at the
  measured concentrations (WB04: 1.1e20 cm^-3). Not fitted, to avoid tuning to the
  validation set. Gain windows barely move (envelope over pump intensity), but the threshold
  pump power for these levels is optimistic.
- MLIP-MD (MACE-OMAT-0) Er:ZBLAN sites reproduce the stochastic-site spectra (splitting 291 vs
  269 cm^-1, 1.5 um FWHM 40 vs 36 nm), though from a very short quench.


## Level- and device-level validation (`scripts/validate_levels.py`, `scripts/validate_devices.py`)

Primary-paper numbers in `gainmat/refs/` (fiber_spectro_refs, tm_silica_refs, zblan_yb_refs,
midir_ho_refs); results in `forward_results/level_validation.*` and `device_validation.json`.

- Lifetimes, ZBLAN at 0.1 mol% (Wetenkamp 1992) and others: Ho 5I7, 5I6, 5S2; Er 4I13/2, 4I11/2,
  4S3/2; Nd 4F3/2; Tm 3H4 all within 30 %. Misses: Ho 5F5 (model 7 us vs 0.29 ms; the measured
  decay most likely follows its 0.30 ms feeder 5S2), Pr 1G4 (x1.9).
- Tm 3H4 in silica / aluminosilicate / germanosilicate (Faure 2007, the gap-law calibration):
  x1.6 / 1.1 / 1.1.  Radiative lifetimes (11 values) within 25 %.
- Ho 5I7 in Al-silica: model 2.1 (silica law) - 5.3 ms (aluminosilicate law) vs Forster tau0
  1.9 ms (Le Gouet 2020); the 0.7 ms 1/e lifetime is Ho-Ho transfer.
- Peak cross sections (18): mostly within 30 %; Ho 5I7 emission in silica x0.5, Pr 1G4 x2.5.
- Stark levels: Er 4I13/2 lowest level +3 cm^-1; Tm 3H4 / 3H5 edges -64 / -154 cm^-1;
  Yb splittings half the measured (Pask 1995).  Doubling the Yb field fixes the splittings and
  the 1.2 um band edge but broadens the 975 nm zero-phonon line: not applied.
- Demonstrated lasers (Jackson 2012, 13): 8 inside a model gain window.  Misses need physics
  not yet in the model: Tm->Ho sensitisation, Pr deactivation co-doping, pumping in an
  absorption band tail (1064 nm into Tm), and a cryogenic self-terminating Ho line; Ho 3.22 um
  is 80 nm outside the window.
- Devices: Ho:silica laser slope 0.57 vs 0.49, threshold 9.6 vs ~13 W (observed lifetime);
  Er 3.5 um dual-pump laser slope 0.32 vs 0.25 (vs incident P2); Pr:ZBLAN amplifier 10 dB vs
  5.2 dB with the measured lifetime (the remaining factor is the x2.5 cross section);
  Tm cladding-pumped lasers slope 0.17-0.25 vs 0.21-0.46: the model's 3H4 absorption band is
  too narrow (sigma at 805 nm 1.7e-21 vs 4.8e-21), and the YAG-derived cross-relaxation
  coefficients give too little two-for-one; Yb lasers fail at 840 / 1102 nm (band tails and
  Yb Stark splitting).
- Common thread of the remaining misses: band wings.  The Gaussian homogeneous line has no
  phonon sideband and the stochastic sites give too narrow silica bands, so pump tails and
  emission edges are under-predicted.


## Explicit physics added in the last round, and the full scorecard

See `forward_results/SCORECARD.md` (every check in one place) and `forward_results/residuals.jsonl`
(84 model/measured pairs, the target table for an ML correction layer).

New physics: optical pumping on every transition (ground- and excited-state absorption) and
dual-wavelength pumping; fibre amplifier (with forward/backward ASE and noise figure) and
Fabry-Perot laser models; ion-ion energy transfer within one species and between species
(`fiber.CoDopedModel`); temperature-dependent multiphonon decay; per-site nephelauxetic
energy scatter (calibrated on the Er zero-phonon-line width); optional one-phonon sideband.

New checks: Giles & Desurvire 1991 (Er alpha/g* spectra and sigma_e/sigma_a), Golding 2000
(Er,Pr:ZBLAN energy transfer), nLIGHT LIEKKI Er/Yb cross sections (metrics only), Le Verre
Fluore demonstrated lasers (11/15 inside a model window).

What explicit physics has not fixed, and why the next step is a learned correction:
- band wings and widths: model bands are ~1.5x too narrow and absolute cross sections
  ~0.5-0.9x (silica); several demonstrated lines sit 10-150 nm beyond the model's window edge;
- the sigma_e/sigma_a ratio is 0.1-0.2 decades high in Er:silica (McCumber energy too high);
- Yb Stark structure (half the measured splitting), and lineshape tails used by pumps;
- concentration / OH quenching of the ~2 um levels in silica (needs a measured series).
