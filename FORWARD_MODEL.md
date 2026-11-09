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
   parameterised by the max phonon energy).
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
