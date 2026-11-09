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
- Largest error: lifetimes of the ~2-3 um upper levels. Stark-edge gaps cut the error from
  10-25x to 2-6x (model / measured, ms: Tm 3F4 silica 2.85 / 0.65 (was 6.28); Ho 5I7 silica
  3.21 / 0.7 (6.70); Pr 1G4 ZBLAN 0.21 / 0.11 (0.99); Dy 6H13/2 ZBLAN 3.87 / 0.65 (15.8)).
  What is left is not one global factor (Tm 3F4 needs W x25, Ho 5I7 x5), and fused-silica fibre
  rates sit 10-80x above the multicomponent-silicate gap law that Layne measured. Likely
  causes: fused silica is not Layne's silicate, OH and Tm-Tm / Ho-Ho quenching at the
  measured concentrations (WB04: 1.1e20 cm^-3). Not fitted, to avoid tuning to the
  validation set. Gain windows barely move (envelope over pump intensity), but the threshold
  pump power for these levels is optimistic.
- MLIP-MD (MACE-OMAT-0) Er:ZBLAN sites reproduce the stochastic-site spectra (splitting 291 vs
  269 cm^-1, 1.5 um FWHM 40 vs 36 nm), though from a very short quench.
