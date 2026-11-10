# ML layer: results

The plan is in `ML_PLAN.md`. The ML layer does not predict outputs. It learns the few physical
inputs the explicit model cannot compute, and every rule the explicit model enforces still
holds. All evaluations below are out-of-distribution: the held-out unit is an ion, a source, a
publication or a glass family, never a random row. Each one is compared against the explicit
model or a trivial baseline.

## Verdict per layer

| layer | learned quantity | OOD test | result | adopted? |
|---|---|---|---|---|
| L1 | Judd-Ofelt Omega from composition | publication-grouped CV, leave-family-out, cross-ion | no useful gain over "ion median" OOD; the scatter between publications (0.25 / 0.10 / 0.13 dex) is as large as the composition effect | **no**: measured Omega where available, else ion median +- ~0.3 dex |
| L2 | host microstructure theta (cf_scale, disorder, neph_std, gamma_h) | leave-one-ion-out within a host; cross-host | **ZBLAN transfers**: held-out loss 233 -> 133, 4 of 5 ions better, held-out band widths 193 -> 68 (20 bands). Aluminosilicate (233 -> 258) and germanosilicate (138 -> 171) **do not** | ZBLAN only (adoption gate) |
| L3 | extrinsic quenching of ~2 um levels | leave-one-out; Tm-only -> Ho | W_q = 1695 (N/1e20)^0.43 s^-1, silica family, 4800-6300 cm^-1; fluorides none. Fitted on Tm alone it predicts Ho 5I7: 0.66 ms vs 0.70 measured (single-ion 5.35) | yes |
| L4 | per-class error model | leave-one-source-out | oxide integrated cross sections are biased -0.19 dex, and correcting helps out of source (0.205 -> 0.138 dex). Lifetimes: 0.43 dex scatter, no correctable bias. Radiative lifetimes: 0.065 dex | as error bars |
| UQ | gain-window probability | — | Omega spread (measured or predicted) x theta spread from L2 folds | see `ml_results/uq/` |

## L2 in detail: what transfers and why

`ml_results/l2/SUMMARY.json`, `l2_summary.png`. Held-out loss split by what each residual
measures (sum of squared z-scores, explicit -> learned from the other ions):

| host | width | line shape | peak shift | splitting | sigma_e/sigma_a ratio |
|---|---|---|---|---|---|
| ZBLAN | 193 -> 68 (n=20) | 21 -> 18 | 20 -> 48 | — | — |
| aluminosilicate | 24 -> 42 (n=10) | 117 -> 135 | 58 -> 43 | — | 34 -> 39 |
| germanosilicate | 6 -> 4 | 61 -> 67 | 0 -> 1 | 49 -> 67 | 21 -> 31 |

- **ZBLAN.** The explicit bands are too narrow, and the extra broadening is shared by Er, Ho,
  Nd and Tm, so it is a property of the glass. A theta learned from any four ions widens the
  fifth ion's bands correctly. The one fold that gets worse, Dy, gets worse through peak position
  (a single digitised mid-IR curve), not width.
- **Silica family.** The explicit widths are already right. What remains is line shape, peak
  position, the McCumber ratio and the Yb splitting. These depend on each ion's Stark
  structure, so one host-wide theta cannot fix them. Fitting theta to the training ions only
  bends it toward their particular errors, and that hurts the held-out ion. Adoption is gated:
  theta goes into the forward model only if, summed over the folds, the held-out loss falls and
  most folds improve.
- **Identifiability.** Only the total broadening is constrained. Across folds, disorder ranges
  0.2–1.4 and gamma_h ranges 43–108 cm^-1, trading off against each other. In ZBLAN, neph_std
  sits on its upper bound. Read the ZBLAN theta as an effective broadening, not as separate
  measurements of disorder and homogeneous width.
- **Cross-host.** Each family's theta is best on its own data. The ZBLAN theta is the worst
  choice for the oxides (aluminosilicate 156 own vs 328 with ZBLAN's), so theta encodes host
  identity rather than a generic correction.

## The largest remaining errors point at Stark structure in silica, not at ML

- **Tm 3H6 -> 3F4 absorption peak.** Model 1720 nm (aluminosilicate), measured 1633 nm
  (Jackson 1999). The 790 nm peak is right (8.96e-21 vs 8.6e-21 cm^2), but the band is narrow,
  so at the 805 nm pump the model has 1.7e-21 against ~4.8e-21. That alone accounts for most of
  the Tm laser threshold error (model ~10–13 W vs 5–7 W measured).
- **Er sigma_e/sigma_a 25–45 % high (Giles).** Same signature: the thermal weighting of the
  ground-manifold Stark structure.
- **Yb splitting.** It needs about 2x the Er-calibrated crystal field (Pask), and Yb is known to
  be anomalous.

Next step: an ion-resolved crystal-field correction (cf_host x cf_ion), tested by
leave-one-host-out. It is not part of this round.

## Devices (`forward_results/device_validation.json`)

ML-L3 replaces the "observed lifetime" fudge with a predicted quenching rate:

| device | measured slope / threshold | observed-lifetime variant | ML-L3 |
|---|---|---|---|
| Ho:silica (Le Gouet 2020) | 0.49 / 13 W | 0.56 / 9.7 W | 0.56 / 10.0 W |
| Tm fibre 1 (Jackson 2003) | 0.21 / 7.4 W | 0.17 / 11.8 W | 0.17 / 11.9 W |
| Tm fibre 3, +CR | 0.415 / 6.0 W | 0.24 / 9.5 W | 0.24 / 12.2 W |
| Tm fibre 4, +CR | 0.46 / 5.2 W | 0.25 / 8.0 W | 0.25 / 10.9 W |

ML-L3 reproduces the hand-tuned variant without using that device's lifetime. For the heavily
doped Tm fibres 3 and 4 it predicts more quenching than the 0.65 ms assumed there. The Tm
thresholds are dominated by the 805 nm absorption error above.

ZBLAN devices with the learned host theta on (`forward_results_ml/device_validation.json`):

| device | measured | explicit | + ML |
|---|---|---|---|
| Pr:ZBLAN amplifier, net gain at 1310 nm | 5.2 dB | 20.4 dB (11.0 observed-lifetime) | 10.4 dB (8.1) |
| Er:ZBLAN dual-pump 3.5 um laser, single-ion | slope 0.254, 0.26 W at 2.4 W | 0.29, 0.74 W | rolls over at ~0.2 W |
| same, + Golding 2000 ETU1/ETU2 | slope 0.254, 0.26 W | 0.33, 0.81 W | 0.14, 0.38 W (saturating) |

The dual-pump result is mixed. With the broader ZBLAN bands, the 4I11/2 -> 4F7/2 ESA at
985 nm (a real Er transition) rises 26x, from 9.7e-24 to 2.6e-22 cm^2. That drains the 4I11/2
reservoir. Absolute power gets closer to the measurement, but the shape gets worse.

While checking this I fixed the laser solver. Its damped sweeps could stall short of gain
clamping (round-trip log gain 0.026 left over). It now falls back to bracketing the
round-trip-gain zero. Results that had already converged are unchanged.

## Full validation suite: explicit vs explicit + ML

To rerun: `GAINMAT_ML=1` + the validation scripts write into `forward_results_ml/`, then
`scripts/compare_ml.py` writes `forward_results_ml/COMPARISON.md`. Rows are split into three
subsets:

- **independent**: nothing in the ML layer saw these measurements.
- **same source as L2 widths**: absolute cross sections taken from the same ZBLAN spectra whose
  widths trained theta. The absolute value was never used, but it is not independent either.
- **in-sample**: lifetimes that trained L3.

| check | subset | n | rms dex explicit -> +ML | within 30 % | within 2x |
|---|---|---|---|---|---|
| peak cross section | independent | 11 | 0.197 -> 0.179 | 64 -> 64 % | 82 -> 91 % |
| peak cross section | same source as L2 | 7 | 0.140 -> 0.092 | 43 -> 57 % | 100 -> 100 % |
| lifetime | independent | 24 | 0.389 -> 0.387 | 62 -> 67 % | 88 -> 88 % |
| lifetime | in-sample (L3) | 4 | 0.610 -> 0.400 | 25 -> 75 % | 25 -> 75 % |
| integrated cross section, Stark splitting, radiative lifetime | all | 38 | unchanged (as they should be: theta moves widths, not integrals) | | |

**The cleanest transfer test is Pr3+ in ZBLAN.** Pr was not among the five ions used to learn
the ZBLAN theta. With theta switched on:

- the 1G4 lifetime goes from 1.95x to 1.31x measured, because the crystal-field scale moves the
  Stark-edge gap that sets the multiphonon rate;
- the 1310 nm peak cross section goes from 2.47x to 1.91x (Ohishi 1991).

Dy 6H13/2 also improves, from 6.0x to 4.0x. The one row that gets worse is Ho 5F5 in ZBLAN, from
0.03x to 0.02x; the explicit model already gets it wrong by 30x.

## Crystal field across glasses (round 2, `gainmat/ml/cf_ranks.py`, `ml_results/cf/`)

New data: measured total Stark splittings (`gainmat/refs/stark_glass_refs.py`):

- Yb³⁺ ²F₇/₂ at room temperature in seven glasses (Yang 2014; Cao 2018).
- Yb³⁺ ²F₇/₂ at low temperature in phosphate (Zhang 2015).
- Er³⁺ ⁴I₁₅/₂ by FLN at 4.2 K in fluorozirconate, phosphate and silicate (Zemon 1991).

Findings:

1. **No Yb-specific crystal-field anomaly.** Modern Yb ²F₇/₂ splittings in silicates are
   874–965 cm⁻¹, not Pask 1995's approximate 1490. With the Er-calibrated field the model is
   within about 15 % for Yb in oxide glasses. The earlier "Yb needs 2x" was driven by one
   approximate figure.
2. **Ion transfer works in phosphate.** A host field scale set from Er alone (0.84) predicts
   Yb at 594 cm⁻¹ against 589 and 596 measured. In fluorozirconate it is partial: Er wants
   1.12x, Yb 1.45x.
3. **Two independent routes agree for ZBLAN.** The field strength learned by L2 from band
   widths and shapes (cf_scale 1.21) matches the one from Stark splittings (joint 1.24).
4. **Oxygen rank balance does not transfer.** A learned B2:B4:B6 balance (1.22 : 1 : 0.86) fits
   better in-sample, but in the phosphate leave-family-out test it does worse than the handbook
   balance. It is not adopted.
5. **Field strength alone cannot move the Tm ³F₄ absorption peak.** At 1.5x field the peak is
   still at 1705 nm, against 1633 measured.
6. **Adopted (ML on):** per-family field scales from the splittings for silica 1.12,
   aluminosilicate 1.04, phosphosilicate 1.27, germanate 1.29, phosphate 0.88 and tellurite
   0.98 (`ml_results/cf/adopted_family_scale.json`).
    - Only phosphate has a transfer test behind it. The others are single-ion calibrations.
    - ZBLAN keeps its L2 value.

### Round 2, continued: what limits the line shapes in silica, and the measured-absorption mode

- **Stark-resolved intensities (not adopted).** `spectra.aniso_S` computes each Stark
  component's strength from the odd-rank crystal field of the same site. It is renormalised to
  the measured Ω, and the multiplet sums reproduce Judd–Ofelt to better than 1 %. The free
  amplitude ratios fitted on Er/Yb/Ho line shapes do not transfer:
    - held-out Tm: 120 → 171;
    - held-out Er: 109 → 111.

  The fitted ratios also differ completely between folds. Off by default.
- **Free-ion level energies are right.** In ZBLAN, 65 multiplets of Pr, Nd, Dy, Ho, Er and Tm
  match Wetenkamp 1992: median deviation 57 cm⁻¹ (0.36 %), 92 % within 2 %
  (`forward_results/level_energy_validation.json`). Tm sits 0.8 % low. A Tm energy scale fixed
  on ZBLAN makes the silica 3H4/3H5 peaks worse, so it is not the cause of the silica 3F4 error.
- **Neither vibronic sidebands nor field strength fixes the silica peaks.**
    - A sideband with S up to 0.5 leaves every peak in place.
    - A larger Yb field reproduces the measured 1064/1102 nm emission tail, but the emission peak
      moves to 994 nm (measured 975).

  The remaining silica errors need measured Stark levels for Tm (and Yb in germanosilicate) to
  diagnose further.
- **Measured-absorption mode (`glass.apply_measured_absorption`, override
  `measured_absorption=`).**
    - What it does: it replaces a band's absorption with the measured curve and derives the
      emission by McCumber, using ε from the model's Stark structure (or uses a measured emission
      if one is given). It also rescales the band's radiative rate to the measured integral.
    - Test (`gainmat.ml.from_absorption`, 5 glasses with measured abs + em of the same band): the
      emission predicted from the measured absorption beats the model's own emission.
        - shape rms 0.10 vs 0.19;
        - peak position 0.4 vs 7 nm.

      Refitting the host parameters to the absorption helps the model (0.19 → 0.15) but not
      past McCumber. The emission/absorption ratio error (about 0.1 dex) is the same for all
      methods, which locates it in ε (Stark structure).
- **The Yb laser failure is spectroscopy, not the device model.** With Pask's own measured Yb
  spectra the 1102 nm laser gives slope 0.85 and threshold 51 mW (measured 0.90 / 30 mW). The
  composition-only model never reaches threshold, because its emission at 1102 nm is
  essentially zero against 2.2e-21 cm² measured. Ho is unchanged with measured spectra (slope
  0.60); its error is the lifetime, which L3 handles.
- **Independent L3 check (Walsh & Barnes 2004, never used in fitting).** Tm:silica, 1.12e20 cm⁻³:
  ³F₄ lifetime 0.42 ms measured; L3 predicts 0.46 ms (single-ion model 2.9 ms). In Tm:ZBLAN
  (0.36e20) the model gives 8.1 ms against 11.2 ms measured, with no quenching, as expected.
- **MD-derived glass sites do not fix the silica peaks either.** These are the 80 Er sites from
  the MACE MLIP melt-quench on branch forward/md-sites, with the field scaled so Er
  ⁴I₁₅/₂ = 400 cm⁻¹.
    - Tm ³F₄ absorption 1754 nm (random polyhedra 1731; measured 1633).
    - Yb ²F₇/₂ splitting 728 cm⁻¹ (measured 874–965).
    - No Yb emission at 1102 nm.

  The superposition crystal field itself, not the site geometry, is the likely limit for Tm and
  Yb in silica.
- **Usable tool.** `scripts/predict_with_absorption.py --ion --family --band --abs file.csv
  [--em file.csv] [--ml]` writes the calibrated spectra, gain windows and lifetimes, plus a
  comparison plot. Example, Ho:silica with Le Gouët's absorption and the ML layer on: ⁵I₇
  lifetime 0.69 ms (measured 0.70), net-gain window 1807–2214 nm.
- **Judd–Ofelt from measured bands (`gainmat/forward/jo_fit.py`).** Ω are fitted to measured
  band integrals with the model's own band kernels (linear in Ω), so the fit is consistent with
  how the model turns Ω into cross sections.
    - Round trip, Tm:silica: (6.23, 1.91, 1.36) recovered as (6.20, 1.89, 1.37).
    - With peaks only (Jackson 1999 peaks times model shapes) the fit is unreliable: Ω₄ goes to
      0. Full measured curves are needed.
    - The CLI accepts several `--band/--abs` pairs, fits Ω when they constrain all three, and uses
      the fitted Ω for every other transition.
- **The Tm laser error is the pump band, then cross-relaxation (`scripts/validate_tm_cr_transfer.py`).**
    - Pump band: I replaced only the 3H6 → 3H4 absorption with a band built from two published
      numbers (8.6e-21 cm² peak at 790 nm, about 4.8e-21 cm² at 805 nm; 32 nm FWHM). Combined
      with L3 + cross-relaxation, fibre 1 (0.7 wt %) goes from threshold 11.4 W to 7.4 W
      (measured 7.41), slope 0.27 (0.21).
    - Cross-relaxation, learned on fibre 4 only: 8x the Jackson & King coefficient (linear in N),
      giving slope 0.49 / threshold 5.7 W (measured 0.46 / 5.24).
    - Predicting fibre 3 with that coefficient: slope 0.415 against 0.415 measured, threshold
      7.2 W against 6.0 (with the literature coefficient: 0.29 / 9.9 W).
    - It does not transfer to the dilute, low-Al fibre 1 (slope 0.33 vs 0.21). Cross-relaxation
      depends on more than the mean Tm density: the Al:Tm ratio, which sets the Tm–Tm distances,
      differs (5 vs 9–10).
- **The reduced matrix elements match the standard tables.** The model's own U^(λ)² agree with
  Carnall (Er) and Krupke (Nd) to within 1 % for every tabulated transition. So the Judd–Ofelt
  layer is consistent with how literature Ω are derived, except where a paper's own U² differ.
- **Pr³⁺ is the exception.** With Ohishi's Ω, the model's 1G4 radiative rates are uniformly
  about 1.37x the rates Ohishi tabulates (Table 1), and τ_rad is 0.71x theirs.
    - Translating their Ω into the model's basis (NNLS on their six 1G4 rates) gives
      (1.63, 3.26, 3.19) instead of (1.43, 4.22, 4.87).
    - That lowers the 1310 nm peak σe by 32 % (explicit model 9.2 → 6.3e-21; measured 3.5e-21).
    - Not applied by default: it calibrates to one paper's derived rates. It explains about a
      third of the Pr amplifier overestimate.
- **Where the Er σe/σa error comes from (Giles, +0.13 dex).** The zero-line energy is right:
  E(⁴I₁₃/₂ lowest) − E(⁴I₁₅/₂ lowest) = 6548 cm⁻¹, against Desurvire's 6540–6548. So the error
  is in the partition-function ratio, i.e. the spacing of the Stark levels.
    - The model's lowest ⁴I₁₃/₂ levels sit at 0, 19, 41, 67 cm⁻¹; Desurvire's FLN gives 0, 100,
      167, 226 (4 of 7 levels identified).
    - The gap persists with 4x angular disorder and with any tested B2:B4:B6 balance.
    - Either the superposition model misses the low-symmetry ⁴I₁₃/₂ pattern, or the FLN set skips
      closely spaced levels. Crystals such as Er:YAG also have first ⁴I₁₃/₂ spacings near 20 cm⁻¹.
      Unresolved.
- **Yb long-wavelength emission: sidebands don't transfer.** Both real Yb fibres (Pask's Ge-silica
  and a vendor aluminosilicate, the latter a private check) emit 5–9 % of peak at 1100 nm. The
  model emits essentially zero there.
    - A one-phonon sideband fitted on Pask's tail (S = 0.2, tail rms 0.079 → 0.026) makes the
      held-out vendor tail worse (0.047 → 0.067).
    - The whole-curve error (rms 0.22–0.33) is dominated by the zero-phonon region, which no
      tested mechanism fixes.

  Not adopted. For Yb, the measured-spectrum mode is the reliable route.
- **Uncertainty grid, round 2 (`ml_results/uq2/`, ML layer with the per-glass field scales).**
  All 10 known media still reach P = 1.00 inside their measured band.
    - Ho:silica 2.1 µm band coverage at P ≥ 0.5 rises from 0.84 to 0.95.
    - Yb:aluminosilicate stays at 0.64: the model's windows stop near 1095 nm against a measured
      band to 1200 nm, the missing Yb long-wavelength emission above.

### Round 2, late additions: Simpson thesis (Tm in silica) and Malouf 2016 (Er:ZBLAN model parameters)

Numbers are in `gainmat/refs/simpson_malouf_refs.py`; none of them were used in any fit.

- **The L3 concentration law failed an independent test and was corrected.** Simpson measured
  nine Tm aluminosilicate fibres at N = 7e17–1.3e19 cm⁻³; their ³F₄ lifetimes are 384–573 µs
  with no concentration trend.
    - The fitted N^0.43 rule predicted 1.0–2.2 ms (rms log error 1.32).
    - An N-independent rate from the same original points gives 0.18.
    - So the extra ³F₄ decay in silica is a host property, not concentration quenching; Al₂O₃
      lengthens it (384 µs at 1.9 % → 573 µs at 9 %).
    - The exponent is now fixed at 0 (`quench.P_FIXED`). Cost: the Ho ⁵I₇ prediction from Tm
      data moves from 0.66 to 0.50 ms (measured 0.70).
    - Details: `ml_results/l3/simpson_test.json`.
- **Tm absorption bands vs Simpson's tables (aluminosilicate).**
    - ³H₄: peak 789 nm (787.5); FWHM 13 vs 18.7 nm, too narrow by the familiar ~1.4x.
    - ³H₅: peak 1206 nm (1209).
    - ³F₄: FWHM 170 nm (181), peak 1710 nm against 1632–1660.
- **The measured ³F₄ peak moves strongly with composition.** It sits at 1577–1580 nm in Ge-silica
  and 1632 → 1641 → 1660 nm for Al₂O₃ 3.4 → 4.8 → 9 mol %. The model puts it at 1710–1726 nm in
  both. Raising the model's field moves the peak bluer, whereas more Al moves the measured peak
  redder. So the ³F₄ position is not set by field strength in the way the superposition model
  has it. This composition series is the test any future Tm crystal-field model must pass.
- **The ML-layer Er:ZBLAN cross sections agree with Malouf et al.** Malouf et al. modelled the
  same Henderson-Sapir laser (independent fit, H2014 column).
    - The 985 nm 4I11/2 → 4F7/2 ESA: explicit model 9.7e-24 cm² (20x low); with the ML layer
      2.6e-22 against Malouf's 2.0e-22. The 26x increase that made the dual-pump laser saturate
      is correct.
    - 1973 nm ESA: 3.8e-21 (3.0e-21); 3.47 µm emission: 1.0e-21 (1.2e-21); 985 nm GSA: 7.8e-22
      (9.3e-22).
- **The dual-pump laser output stays off with or without the ML layer.** I used Malouf's
  device inputs (1.7 mol %, not the 1 mol % used before, and R_OC 0.87). Both versions lase from
  P2 ≈ 0.3 W, against a measured threshold near 1.4 W:
    - explicit: slope 0.36;
    - ML: slope 0.17;
    - measured: 0.254.

  The remaining gap is in the energy-transfer and loss terms (Malouf's model needs a large
  W5362 cross-relaxation), not in cross sections.
- **Yb Stark levels in a silica fibre (Newell 2007, 77 K).**
    - Measured: ²F₇/₂ at 0, 425, 530, 865 cm⁻¹; ²F₅/₂ first gap 630 cm⁻¹; zero-phonon line
      977 nm (Wang 2014: 975–976).
    - The total ground splitting (865) agrees with Yang's 874 and with the model's silica scale
      (802 with the ML layer), not with Pask's 1490. The highest ground level gives the 1069 nm
      line, so the 1100 nm emission tail is thermal or vibronic, not another Stark level.
    - The model's lowest gaps are too small: ²F₇/₂ second level at 240–300 (measured 425);
      ²F₅/₂ gap 240–350 (630).
    - Rank balance does not fix this. Distorted sites move it toward the data (MD sites: 363/427;
      a 30 % uniaxial distortion: 418/479) but not all the way, while leaving the Er ⁴I₁₅/₂ pattern
      (0, 30, 55, 80 vs Desurvire 0, 30, 55, 129) intact.
    - Real sites are more axially distorted than the jittered polyhedra. This is a site-geometry
      problem (better MD), not a data problem.
- **Lincoln 1991 (Tm silica; old level labels).**
    - FLN at 9 K resolves ground ³H₆ components near 0, 250 and 430 cm⁻¹ and more above
      500 cm⁻¹ in Ge-silica, but not a full set.
    - The 3H4 lifetime (about 20 µs) shows no concentration dependence from 20 to 12 000 ppm,
      consistent with the corrected L3.

## Round 3: MD glass sites and a crystal field learned across ions (interim, `ml_results/md_cf/`)

**Sites.** Pedone-potential melt-quench MD of Er-doped silica glasses, run as four cloud sessions
(branches `md2/*`). Quench: 10 K/ps to 300 K, about 3000 atoms, 8 Er per box. Interim site sets:

| glass | site snapshots | Er–O (Å) | CN |
|---|---|---|---|
| pure silica | 80 | 2.34 ± 0.23 | 4.95 |
| 3.4 % Al₂O₃ | 320 | 2.36 ± 0.24 | 5.0 |
| 10 % GeO₂ | 160 | 2.38 ± 0.28 | 4.9 |

Measured: Er–O ≈ 2.25 Å, CN ≈ 6. The glass network checks pass: Si–O 1.61 Å, O–Si–O 109°,
Si–O–Si 153°.

**The MD sites alone, with the handbook crystal field, don't fix anything**
(`forward_results/md_sites_validation_interim.json`).
- The Yb ²F₅/₂ gap rises from about 245 to about 370 cm⁻¹ (measured 630).
- The Tm ³F₄ absorption peak stays at 1730–1750 nm (measured 1578–1660). The measured
  Ge-silica peak is bluer than the Al-silica one; the model has it redder.

**On the MD sites the rank balance becomes identifiable** (`gainmat/ml/md_cf_fit.py`).
Parameters: overall scale s and multipliers r2, r6 on the handbook B2, B6. Targets:
- Er ⁴I₁₅/₂ levels and totals (Desurvire, Zemon);
- Yb ²F₇/₂ levels and ²F₅/₂ gap (Newell);
- Tm ³F₄ absorption peak and width (Simpson Al 3.4 %) and emission peak (Agger).

Leave-one-ion-out results (χ² of the held-out ion, handbook → fitted on the other two):

| held out | fitted (s, r2, r6) | held-out χ² | held-out quantities |
|---|---|---|---|
| Tm | 0.76, 2.36, 1.86 | 135 → **43** | Al-silica abs 1742 → 1689 nm (1632); em 1789 → 1804 (1825); Ge-silica 1748 → 1660 (1578), now **bluer than Al-silica, as measured** |
| Yb | 1.03, 1.80, 1.67 | 41 → **31** | ²F₅/₂ gap 384 → 497 (630); ²F₇/₂ 0/425/706/955 (0/425/530/865) |
| Er | 1.10, 1.43, 2.86 | 22 → 73 (worse) | ⁴I₁₅/₂ total 646 (420) |

- Two of three held-out ions improve. Er fails because Er is the only ion that pins B6: Yb ²F₇/₂
  has almost no B6 sensitivity, so without Er, r6 drifts. That is a gap in what the data can
  identify, not evidence against the rule.
- All three ions agree on B2 being about 1.4–2.4x the handbook value relative to B4: real oxygen
  sites in silica carry a stronger rank-2 (axial) field than the crystal-derived intrinsic
  parameters assume.
- Never fitted: the Er σe/σa vs Giles error (+0.12 → +0.12–0.14 dex) and the Yb zero-phonon line
  (966–968 vs 976 nm) don't change. Those still need something else.
- The composition trend direction comes out right in 2 of the 4 fits (minus-Tm and minus-Yb), so
  it is fragile.

Final site sets, including 9 % Al₂O₃ as another never-fitted composition, will replace these
interim results.
- **Full line shapes for all four silica-family ions** (curve residuals: mean shape z, mean
  |peak-shift| z; Ho was never fitted):

  | sites / crystal field | Er | Yb | Tm | Ho |
  |---|---|---|---|---|
  | random polyhedra, handbook | 2.95 / 0.41 | 4.77 / 0.27 | 3.54 / 2.22 | 3.02 / 0.67 |
  | MD, handbook | 3.10 / 0.35 | 3.54 / 0.27 | 6.09 / 4.22 | 4.32 / 1.87 |
  | MD, fitted on all three | 3.17 / 0.61 | 3.63 / 0.27 | 3.58 / 1.67 | 3.16 / 2.00 |
  | MD, fitted on Er+Yb | 2.95 / 0.47 | 3.35 / 0.27 | 3.82 / 1.67 | 4.07 / 1.42 |

    - Yb shapes improve with MD sites, and Tm peak positions improve with the learned field.
    - Ho shifts get worse and the Er shapes don't move.
    - Net: the learned field fixes Stark-level patterns and peak positions it was aimed at, but
      not whole line shapes across ions.

### Round 3 final: MD sites on the final site sets, with sampling noise (`ml_results/md_cf_final/`)

**This replaces the interim section above and an intermediate update that claimed every held-out
ion got worse.** The two earlier conclusions contradicted each other because each was a single
draw of 60 sites. The final site sets (`md_sites/`, all from the md2/* branches):

| glass | sites (Er × frames) | CN | Er–O (Å) | quench |
|---|---|---|---|---|
| SiO₂ | 320 | 4.28 | 2.28 ± 0.18 | 10 K/ps |
| 3.4 % Al₂O₃ | 480 | 5.06 | 2.36 ± 0.24 | 10 K/ps |
| 9 % Al₂O₃ | 160 | 5.22 | 2.34 ± 0.19 | 15 K/ps |
| 10 % GeO₂ | 320 | 4.54 | 2.31 ± 0.18 | 10 K/ps |

Each fold was fitted on three independent 60-site subsamples (`--offset 0/1/2`). Held-out χ²,
handbook → fitted:

| held out (fitted on the other two) | subsample 0 | subsample 1 | subsample 2 | verdict |
|---|---|---|---|---|
| Er | 24 → 10 | 25 → 14 | 23 → 9 | **better in 3/3 (about 2×)** |
| Tm | 125 → 115 | 183 → 167 | 126 → 165 | mixed, no real gain |
| Yb | 42 → 38 | 38 → 60 | 37 → 33 | mixed |
| (in-sample, all three ions) | Er 24→10, Yb 42→24, Tm 125→23 | Er 25→13, Yb 38→26, Tm 183→36 | Er 23→10, Yb 37→23, Tm 126→48 | always better |

The all-ion fit settles at s 0.87–1.03, r2 1.6–2.1, r6 1.3–1.7, so the data consistently want a
**stronger rank-2 field (about 2×)** than the handbook balance. That rank-2 balance, learned from
Yb + Tm, predicts the Er Stark pattern. But Tm cannot be predicted from Er + Yb, and Yb cannot be
predicted from Er + Tm. The folds without Tm swing to s ≈ 1.5, r2 ≈ 0.7, an alternative minimum
for Er + Yb alone. **Not adopted:** it passes for one ion of three, and the adoption gate needs
most folds.

**Tm composition trend: not reproduced.** The ³F₄ absorption peak was computed on four disjoint
40-site subsets per glass, with the handbook field and one Er-calibrated scale (1.175):

| glass | model (nm, mean ± s.e.) | shift vs Al 3.4 % | measured shift (Simpson) |
|---|---|---|---|
| Al 3.4 % | 1747 ± 7 | 0 | 0 (1632) |
| Al 9 % | 1732 ± 3 | −15 ± 7 | +28 |
| Ge 10 % | 1734 ± 1.5 | −13 ± 7 | −54 |
| SiO₂ | 1737 ± 3 | −10 ± 8 | — |

The model's shifts are within about 2 s.e. of zero, and the Al shift has the wrong sign. (The
"right direction" in the intermediate update was sampling noise.) With oxygen treated as one
ligand type, the MD geometries do not separate these glasses. The Tm offset of about +100 nm is
the same as on random polyhedra.

**Other checks** (`forward_results/md_sites_validation.json`, final sites):

| quantity | MD sites | random polyhedra | measured |
|---|---|---|---|
| Yb ²F₅/₂ first gap (cm⁻¹) | 356–415 | 251–339 | 630 |
| Yb ²F₇/₂ pattern rms (cm⁻¹) | 74–116 | 67–130 | — |
| Yb ZPL error (nm) | −6 to −7 | −9 | 0 |
| Er ⁴I₁₅/₂ lowest five rms (cm⁻¹) | 53–63 | 43–66 | — |
| Er σe/σa vs Giles | +0.1 dex | +0.1 dex | — |

MD sites close about a third of the Yb ²F₅/₂ gap deficit and slightly improve the ZPL. They change
nothing else.

**What this means.** Better geometry alone does not fix silica. The two robust signals:
(1) the rank-2 field is too weak relative to rank 4 and 6, which is consistent with the Yb gap
deficit;
(2) the composition dependence has to come from ligand chemistry (bridging vs non-bridging O,
Al–O vs Si–O vs Ge–O), which a one-ligand superposition model cannot express.

### Round 3 last test: ligand chemistry on the MD sites (`ml_results/md_nbo/`, `scripts/md_nbo_*.py`)

The MD files tag each first-shell oxygen as bridging or non-bridging (NBO) and record which
formers it bonds to. The chemistry differs a lot between glasses and lines up with Simpson's Tm
series:

| glass | NBO | O bonded to Al | Tm ³F₄ peak (measured) |
|---|---|---|---|
| 10 % GeO₂ | 66 % | 0 | 1578 nm |
| 3.4 % Al₂O₃ | 50 % | 15 % | 1632 nm |
| 9 % Al₂O₃ | 35 % | 41 % | 1660 nm |
| SiO₂ | 72 % | 0 | — |

New code: `load_md_sites_tagged` renames ligands by chemistry, and a `ligand_factors` override
gives them separate superposition strengths. Grids run over bridging strength a and NBO strength
b (× handbook, 60 sites per glass).

- **The NBO factor cannot be learned from Er + Yb.** Fitting the pure-silica Stark levels gives
  loss 34.3 with (a, b), against 35.8 with one oxygen type. The fit prefers a weaker NBO and
  predicts no Tm composition trend (+5, −15 nm against measured −54, +28).
- **Leave one glass out on Tm (a, b fitted on two glasses):** held-out errors are −93, +57, +20 nm.
  One oxygen type gives −13, +58, −8. The extra parameter does not beat the baseline.
- **The real obstacle is field strength.** With one oxygen type:

  | O strength × handbook | 1.0 | 1.4 | 1.8 | 2.2 | 2.6 |
  |---|---|---|---|---|---|
  | Er ⁴I₁₅/₂ total (cm⁻¹; measured about 420) | 338 | 474 | 611 | 749 | 888 |
  | Yb ²F₇/₂ total (cm⁻¹; measured 865) | 635 | 891 | 1147 | 1403 | 1660 |
  | Tm ³F₄ peak, Ge / Al3 / Al9 (nm) | 1740 / 1756 / 1744 | 1722 / 1715 / 1726 | 1629 / 1694 / 1695 | 1597 / 1701 / 1674 | 1565 / 1685 / 1653 |

  The Er and Yb splittings need a strength of about 1.2–1.4. At that strength the MD sites put Tm at
  1715–1740 nm. Tm needs about 2× more field than Er and Yb on the same geometry. With the
  handbook Tm:Er ratio of intrinsic parameters, no single oxygen field fits all three ions. That
  explains why the leave-Tm-out folds fail, and why a geometry change (random polyhedra → MD)
  cannot fix Tm.
- The peak position also hops between Stark components as the field changes (Al3 at 1.8 → 2.2
  goes 1694 → 1701 → 1685), so single-peak targets are a fragile fitting signal.

**Conclusion.** The Tm problem is in the Tm crystal-field parameters (the handbook Tm:Er intrinsic
ratio, or Tm-specific covalency, as with Yb), not in the site geometry. The direct fix is the one
noted in round 2: measured Tm³⁺ Stark levels in a silicate glass. With those, a Tm-specific
strength factor can be set on levels and tested on Simpson's band positions, instead of being
fitted to them.

## Round 4 (Nov 10): toward a frontier model — what helps, measured against the practical alternative

"Frontier" is defined here as beating what a practitioner would do without the model: look up the
same quantity in other papers on the same ion and glass family. Measurement scatter between labs
sets the floor. Every change below was tested on data it never saw.

### 1. Tm crystal field: two more fixes tried, neither adopted (`ml_results/tm_cf/`)

- **A Tm–O field factor.** The handbook has no Tm oxide data and copies Er's values. A factor of
  about 1.9 puts the ³F₄ band right in Al-silica (absorption 1641 nm, emission 1825 nm, matching
  Agger). It breaks everything held out: ³H₅ moves 1209 → 1189 nm, ³H₄ moves blue, the 805 nm
  pump absorption falls, and the ³F₄ widths double.
- **Tm-specific rank multipliers.** Fitted to one glass's ³F₄ band, the best set (r2 = 2,
  r4 = 0.6, r6 = 0.5) keeps that glass's ³H₄ and ³H₅ positions and widths right. The near-equally
  good sets disagree by up to 100 nm in germanate and tellurite, and no measured data can choose
  between them.
- **Lincoln 1991 9 K line-narrowing levels** (³H₆ at ~0, 250, 430 cm⁻¹, Ge-silica) cannot tell
  1× from 2× field: the model has levels near both values.
- **Next step:** measured Tm Stark levels. Turri 2008 (8–300 K) gives none in numbers. Gruber 1989
  (Tm:YAG) would test the "Tm = Er" handbook assumption in an oxide crystal.

### 2. Measured-level mode: adopted (`spectra.ensemble_spectra(level_anchor=...)`, `refs/stark_anchor_refs.py`)

Each site's Stark levels go through one monotone map, so the ensemble mean equals measured levels
from one source. The site-to-site scatter and the crystal-field eigenvectors (intensities) are
kept. It is on by default for Yb and Er in the silica family (`GAINMAT_ANCHORS=0` switches it
off). For Yb in phosphosilicate, phosphate, germanate and tellurite, only the zero-phonon line is
anchored. Tests use other sources (`scripts/validate_level_anchor.py`, full suite rerun):

| held-out test | model | measured-level mode | measured |
|---|---|---|---|
| Pask Yb absorption shape error (rms) | 0.24 | 0.12 | — |
| Pask Yb emission shape error (rms) | 0.24 | 0.14 | — |
| Yb peak (Pask) | 970 nm | 974 nm | 974 nm |
| Yb emission 1035 / 1064 nm relative to peak | 0.17 / 0.02 | 0.25 / 0.14 | 0.24 / 0.14 |
| Pask 5 m free-running Yb laser threshold | no lasing | 54–56 mW | 45 mW |
| Er σe/σa vs Giles, rms (Al- / Ge-silica) | 0.140 / 0.114 dex | 0.116 / 0.091 dex | — |
| Jackson Er emission shape error / FWHM | 0.18 / 20 nm | 0.12 / 40 nm | — / 40 nm |
| Yb zero-phonon line in phosphate (Wang 2017) | 967 nm | 975 nm | 975 nm |

Not fixed:
- The Pask 1102 nm laser still does not lase. Its emission beyond ~1080 nm needs a vibronic
  sideband; S ≈ 0.2 would reproduce it but was not adopted, because that value is tuned on the
  same curve.
- The Er absorption peak is 2 nm shorter than before.
- Tm is not anchored: three partial levels make the spectra worse.

### 3. Bug fix: Yb Ω in non-silica hosts

Yb has a single transition, so its "Ω" is only a line-strength normalisation. For glass families
without Yb data, the fallback averaged other ions' Ω, making Yb in phosphate 2.5× too strong.
Yb now borrows only Yb data from the nearest family.

| check | before | after |
|---|---|---|
| Wang 2017 phosphate σa / σe error | 0.39–0.50 / 0.29–0.40 dex | 0.04–0.15 / under 0.07 dex |
| Wang 2017 phosphate lifetime (measured 1.17 ms) | 0.56 ms | 1.29 ms |
| Pask Ge-silica radiative lifetime | 1.16× | 0.92× |
| Pask Ge-silica peak σe | 0.65× | 0.83× |

### 4. External benchmark on open-access literature (`refs/oa_benchmark_refs.py`, `scripts/oa_benchmark.py`, `forward_results/oa_benchmark.md`)

- **Data:** 99 records from 18 open-access papers. Tm, Er, Yb, Ho and Nd in tellurite, germanate,
  phosphate, ZBLAN and Al-silica; mostly hosts the model had barely been tested on. Only numbers
  stated in text or tables were used. Three Tm values were computed from a paper's stated linear
  fits and are marked as such.
- **Predictors** (each blind to the record):
  - the model with its own Ω;
  - the model given the paper's own Ω;
  - a leave-one-paper-out lookup of the same quantity in other papers on the same ion and family.

| quantity | model (own Ω) | model + paper's Ω | lookup | model beats lookup |
|---|---|---|---|---|
| radiative lifetime | 0.112 dex (n=20) | **0.018 dex** (n=14) | 0.159 dex | **14/20** |
| Ω | 0.087 dex (n=16) | — | 0.092 dex | 7/15 |
| absorption peak | 0.9 nm (n=8) | 0.5 nm | 1.5 nm | 4/4 |
| emission peak | 12.7 nm (n=18) | — | **1.5 nm** | 3/18 |
| emission FWHM | 27 nm (n=7) | 15 nm | **5.3 nm** | 0/2 |
| peak σa | 0.079 dex (n=7) | — | 0.070 dex | 1/2 |
| peak σe | 0.185 dex (n=21) | 0.137 dex (n=10) | **0.124 dex** | 5/17 |

**What this says:**
- **The radiative physics is at measurement precision once Ω is known** (0.018 dex, i.e. 4 %),
  and beats lookup without it.
- **The limit is Ω itself.** Neither the current defaults (0.087 dex) nor LOMS family medians
  (0.149 vs 0.153 dex per Ω, no paper overlap) beat the scatter between publications. This repeats
  the L1 finding.
- **Band positions and widths are where published data beat the model.** For Er, Nd, Ho and Yb
  the errors are small (positions 1–7 nm). Tm ³F₄ emission is 50–90 nm short outside silica.
  Tm emission peaks also scatter by about 70 nm between papers (reabsorption in thick bulk
  samples).
- **Peak σe** follows from the width errors: given the paper's Ω the error is 0.137 dex, close to
  the lookup's 0.124.
- **Hybrid (proof of concept, not adopted):**
  - The model band is shifted and stretched to the peak and FWHM that other papers report, and
    the model's integrated strength is kept.
  - It equals the lookup on positions by construction.
  - It improves σe only where another paper reports a width (2 Nd tellurite rows, 0.137 → 0.077
    dex median over 10 rows), too few to claim.
- **Also tested, not adopted:**
  - A tellurite disorder scale: shrinking it narrows Tm a lot but Nd hardly at all. No shared
    width setting, as in silica (L2).
  - One questionable record: the Ho fluorotellurite paper's Ω give 3.9 ms through the standard
    formula, against its quoted 7.98 ms.

### 5. Usability

- **`gainmat.forward.predict.predict(ion, family | composition=..., omegas=...)`** returns, in one
  call:
  - main bands, lifetimes and gain windows;
  - an empirical uncertainty for each quantity, from the validation residuals;
  - the provenance of each input: measured levels / measured Ω / ML layer / pure model;
  - the known biases found in validation.
- **`gainmat.forward.composition.family_for`** maps a composition (mol%, or a formula such as
  `70TeO2-15K2O-15Nb2O5`) to the model's glass family. It classifies 536 of 928 LOMS glasses;
  borates, bismuthates and oxyfluorides are not covered.
- **ML layer rerun:** with anchors on, it still helps out of sample. Lifetimes within 30 % go from
  62 % to 71 %; peak σ error goes from 0.194 to 0.179 dex.

### Where this leaves the model

**Near or at frontier:**
- radiative lifetimes, branching and quantum efficiency, given Ω or a measured absorption band;
- device simulation with measured spectra (earlier rounds);
- Yb and Er spectra in silica fibres;
- absorption peak positions.

**Not frontier:**
- band positions and widths without published data, and Tm ³F₄ in particular;
- Ω for an unmeasured glass. This is limited by the scatter between papers, so the practical
  route is to measure Ω or absorption.

**Next steps:**
- a larger open-access dataset (about 6 records per paper; extraction works);
- measured Tm Stark levels;
- published widths to make the hybrid testable.

### 6. New capability: energy-transfer microparameters from the model's own spectra (`gainmat/forward/transfer.py`, `scripts/validate_transfer.py`)

**Method.** The Förster–Dexter dipole–dipole C_DA, with phonon assistance (standard Kushida /
Tarelho form, Huang–Rhys S0 = 0.31, the family's maximum phonon energy). It is computed from the
model's predicted donor emission and acceptor absorption cross sections, so cross-relaxation,
migration, upconversion and Tm↔Ho transfer can be predicted for a glass nobody has measured.

**Tests.** Published values for real glasses, in 3 papers, never used. "Overlap" means the
authors' Förster–Dexter overlap of their own measured spectra; "decay fit" means a microparameter
fitted to measured decay curves.

| process | glass | published (cm⁶/s) | kind | model, paper's Ω | model, defaults |
|---|---|---|---|---|---|
| Er ⁴I₁₃/₂ migration | ZBLAN | 4.5e-39 | overlap (Huang 2014) | 1.73× | 0.99× |
| Er ⁴I₁₃/₂ migration | ZBYA | 4.0e-39 | overlap (Huang 2014) | 1.76× | 1.12× |
| Tm ³H₄ migration | Bi-germanate | 3.5e-39 | overlap (Han 2016) | 1.84× | 0.97× |
| Tm cross-relaxation ³H₄,³H₆→³F₄,³F₄ | Bi-germanate | 1.6e-39 | overlap (Han 2016) | 0.82× | 0.64× |
| Tm ³F₄ migration | tellurite | 1.27e-38 | overlap (Gebavi 2012) | 1.03× | 1.11× |
| Tm ³F₄ → Ho ⁵I₇ | tellurite | 5.7e-39 | overlap (Gebavi 2012) | 0.94× | 0.95× |
| Ho ⁵I₇ → Tm ³F₄ | tellurite | 1.5e-40 | overlap (Gebavi 2012) | 1.13× | 1.12× |
| Tm cross-relaxation | tellurite | 4.9e-39 | **decay fit** (Gebavi 2012) | 0.22× | 0.24× |

**Results.**
- **Against the published overlap calculations: all 7 within 2×.** Fully predictive, the median
  error is 0.05 dex (rms 0.08). Given the paper's Ω, it is 0.09 dex (rms 0.17). The phonon-assisted
  share also matches: Tm cross-relaxation is mostly one-phonon, as Han reports (84 %).
- **Against the one decay-fitted value: 4.5× low.** A fitted decay folds in other quenching
  channels and diffusion, so this is the honest number for device use. Use it as a lower bound on
  the cross-relaxation rate until more decay-fitted data are collected.
- **Now part of `predict()`:** a standard process set per ion, with this uncertainty attached.

**Cross-relaxation in silica fibres: the physics sets a floor, and the excess points to clustering.**
The model's Tm C_DA (cross-relaxation) and C_DD (³H₄ migration) go into the Burshtein hopping
model, k = π(2π/3)^{5/2} √(C_DA C_DD) N, at N = 4.5e20 cm⁻³:

| quantity | cm³/s |
|---|---|
| model k, Al-silica / silica | 1.3e-17 / 2.2e-17 |
| Jackson & King literature value | 1.8e-16 (8–14× more) |
| value the Jackson & Mossman laser fits needed (round 2) | 1.44e-15 (65–110× more) |

A uniform random distribution of ions underestimates real silica fibres by one to two orders of
magnitude. That is the signature of Tm clustering in silica, and it agrees with three earlier
observations:
- the learned cross-relaxation depends on the Al:Tm ratio (Al disperses rare-earth ions);
- the decay-fitted tellurite value is 4.5× above the overlap value;
- the L3 quenching rule, which is concentration-independent, points to the same host effect.

For devices, the physics value is a lower bound. A cluster fraction (or the Al:Tm-dependent factor)
is the missing input, and it is a glass-processing property, not a spectroscopic one.

**Update: 5 more published values (Cao 2016 Er→Ho in Ba-fluorosilicate; Balaji 2016 Yb→Ho in
tellurite; Sontakke Nd in Ba-Al metaphosphate).** Full table in `forward_results/transfer_validation.json`.

| process | glass | published (cm⁶/s) | kind | model (defaults) |
|---|---|---|---|---|
| Er ⁴I₁₃/₂ → Ho ⁵I₇ | Ba-fluorosilicate | 5.5e-40 | overlap | 0.78× |
| Er ⁴I₁₁/₂ → Ho ⁵I₆ | Ba-fluorosilicate | 6.6e-41 | overlap | **0.29×** |
| Yb ²F₅/₂ → Ho ⁵I₆ | tellurite | 1.4e-40 | overlap | **2.24×** |
| Nd ⁴F₃/₂ migration | Ba-Al metaphosphate | 3.2–3.8e-39 | overlap | 1.26× |
| Nd self-quenching ⁴F₃/₂,⁴I₉/₂ → ⁴I₁₅/₂,⁴I₁₅/₂ | Ba-Al metaphosphate | 3–8e-42 | decay fit | 4.5× |

**Totals.**
- **Published overlap values: 9 of 11 within 2×** (median 0.05 dex, rms 0.21 dex), fully
  predictive.
- **The two misses** are the strongly non-resonant, multi-phonon-assisted processes (Er ⁴I₁₁/₂ →
  Ho ⁵I₆; Yb → Ho ⁵I₆). There the overlap sits in band tails that a Gaussian-line model does not
  get right.
- **Against decay-fitted values: both about 4.5× off, in opposite directions** (Tm 0.22×, Nd 4.5×).
  Decay-fit microparameters depend on the fitting model (Sontakke's range from different models
  is 3–8e-42). Treat predictions as ±0.3 dex for resonant and near-resonant processes and as
  roughly ±0.5 dex otherwise.

### 7. Vibronic sideband: tested on the external benchmark, not adopted

`GAINMAT_SIDEBAND=S` (default 0) switches the sideband on.

| peak σe error, external benchmark | S = 0 | S = 0.1 | S = 0.25 |
|---|---|---|---|
| model, own Ω | 0.185 dex | 0.150 dex | 0.155 dex |
| model, paper's Ω (the clean test of line shape) | 0.137 dex | 0.182 dex | 0.259 dex |

The gain with the model's own Ω only offsets Ω errors. Given measured Ω, the sideband makes peak
cross sections worse. So a uniform sideband is not the missing piece for widths. The far tails
(Yb > 1080 nm; multi-phonon transfer) need an ion-specific vibronic treatment.

## Round 5 (Nov 10, 03:00–08:00): papers supplied by the user (Gruber 1989, Yang 2014, Turri 2008)

### 5.1 Tm:YAG settles the "Tm = Er" question (`scripts/validate_tm_yag.py`, `ml_results/tm_yag/`, `refs/tm_yag_refs.py`)

- **Machinery check.** Gruber's fitted B_kq go through our Tm Hamiltonian and reproduce his
  calculated Stark patterns to a few cm⁻¹; the rms against his measured levels is 16 cm⁻¹ (his own
  fit: 11). Our free-ion centroids of ³F₄ and ³H₄ sit about 70 cm⁻¹ above YAG's.
- **The handbook superposition parameters (copied from Er:YAG, never fitted to Tm) on the real YAG
  site** (8 O at 2.303 / 2.432 Å) reproduce the measured Tm:YAG multiplet splittings:

  | multiplet | predicted splitting (cm⁻¹) | measured (cm⁻¹) |
  |---|---|---|
  | ³H₆ | 782 | 730 |
  | ³F₄ | 679 | 677 |
  | ³H₅ | 546 | 543 |
  | ³H₄ | 523 | 552 |

  Level rms is 55 cm⁻¹. Fitting the Tm intrinsics instead gives ratios to Er of 0.89 / 0.93 / 1.22
  (ranks 2 / 4 / 6).
- **So "Tm = Er" holds in an oxide.** The Tm ³F₄ error in silica is not a Tm crystal-field
  strength problem; this withdraws the round-3 reading that "Tm needs about 2× the field".
- **Line intensities.** For ground-state absorption into the ³F₄ components (Gruber, 15 K), neither
  intensity rule matches the measured pattern (log correlation: isotropic per-component Judd–Ofelt
  0.34, site-derived odd-rank 0.28). The site-derived rule does zero the symmetry-forbidden lines.
- **An equal-strength rule** (`intensity="uniform"`, every component of a multiplet pair equally
  strong) leaves the silica Tm ³F₄ absorption peak at 1711 nm.
- **Conclusion:** the band position comes from where the Stark levels sit within the ground ³H₆
  multiplet, and that is set by the glass site, not by the ion's parameters. Practical route:
  measured-absorption anchoring for Tm.

### 5.2 Yb level patterns in eight glasses (Yang 2014; `scripts/validate_yb_yang.py`, `refs/yang2014_yb_refs.py`)

- **Shape of the pattern** (³F₇/₂ levels relative to the top level d), random-polyhedron crystal
  field, before any per-glass scale:

  | glass (family) | b/d, model vs measured | c/d, model vs measured | rms after scaling to d |
  |---|---|---|---|
  | tellurite | 0.41 vs 0.40 | 0.71 vs 0.70 | 6 cm⁻¹ |
  | fluorozirconate | 0.39 vs 0.41 | 0.71 vs 0.73 | 9 cm⁻¹ |
  | bismuthate (germanate proxy) | 0.38 vs 0.39 | 0.70 vs 0.67 | 10 cm⁻¹ |
  | germanate | 0.38 vs 0.35 | 0.70 vs 0.56 | 61 cm⁻¹ |
  | lead aluminosilicate | 0.36 vs 0.30 | 0.66 vs 0.56 | 52 cm⁻¹ |
  | phosphate | 0.36 vs 0.35 | 0.67 vs 0.52 | 47 cm⁻¹ |

  The crystal-field shape is right in 3 of 6 glasses; in the other 3 the c level sits too high.
  What varies most between glasses is the overall strength: fluorozirconate needs about 1.5×
  (414 vs 624 cm⁻¹).
- **As anchors** (tested on Wang 2017 phosphate, 1005/975 nm cross-section ratios, measured
  emission 0.59 and absorption 0.14):

  | anchor | emission ratio | absorption ratio |
  |---|---|---|
  | none | 1.48 | 0.28 |
  | zero-phonon line only (default) | 0.38 | 0.10 |
  | Yang phosphate levels | 0.29 | 0.07 |

  Yang's full phosphate pattern does worse. Yb levels vary with composition within a family (Wang's
  is a Ba phosphosilicate, Yang's a borophosphate), so the per-family Yang anchors are kept opt-in
  (`auto=False`).

### 5.3 Temperature dependence (Turri 2008; `scripts/validate_temperature.py`)

Turri gives linear fits over 100–300 K of the Tm ³F₄→³H₆ peak emission cross section and peak
wavelength in germanate, Al-silica and borophosphate. The model was run at 100, 200 and 300 K with
each glass's measured Ω. Compared: slope of the peak cross section (%/K).

| homogeneous width vs T | germanate (meas −0.153) | silica (meas −0.064) | phosphate (meas −0.083) | mean slope error |
|---|---|---|---|---|
| constant (default) | −0.033 | −0.034 | −0.032 | 0.45 dex |
| ∝ T | −0.073 | −0.104 | −0.102 | 0.21 dex |
| ∝ T² (Raman theory) | −0.230 | −0.318 | −0.310 | 0.48 dex |

- **Room-temperature values:**
  - silica: σ 3.77e-21 vs 3.61e-21, peak 1796 vs 1808 nm;
  - germanate: σ 5.1 vs 4.0e-21, peak 1779 vs 1870 nm (the known Tm position error);
  - phosphate: σ 3.75 vs 3.12e-21, peak 1782 vs 1836 nm.
- **Peak-wavelength slopes:** silica −0.29 vs −0.21 nm/K; phosphate and germanate too steep.
- **Reading:**
  - with a constant width, the cross section changes 2–5× too little with temperature;
  - scaling the whole room-temperature width as T² overshoots, because that width is an effective
    one that also absorbs static broadening;
  - T¹ fits best, but that was chosen after seeing the data.
- **Not adopted.** The default stays a constant width (`GAINMAT_GAMMA_T_EXP` switches the law). The
  physically right fix is a two-part width: a static part plus a dynamic part ∝ T².

### 5.4 Device-level energy-transfer rates predicted from spectra (`scripts/validate_transfer_rates.py`)

**Method.** The macroscopic rate-equation coefficients (cm³/s) of Er:ZBLAN measured by Golding,
Jackson, King & Pollnau, PRB 62, 856 (2000), Table I (2–14e20 cm⁻³), compared with fully
predictive values:
- C_DA and C_DD from the model's own spectra; phonon absorption is now included as well as
  emission, as the standard Miyakawa–Dexter form requires;
- k = k_hop + k_dir, where k_hop = π(2π/3)^{5/2} √(C_DA C_DD) N (Burshtein migration-assisted) and
  k_dir = (4/3) π^{3/2} √(C_DA / τ_D) (direct Förster, as a rate), with τ_D the model's own lifetime;
- the cross-relaxation uses Golding's lumped thermalised ⁴S₃/₂ + ²H₁₁/₂ pair: channel a through the
  ²H₁₁/₂ Boltzmann share, plus phonon-absorbing channel b.

| process | model / measured at 0.4 / 2 / 8 / 14e20 cm⁻³ |
|---|---|
| upconversion ⁴I₁₃/₂+⁴I₁₃/₂ (the key Er device process) | — / 0.60 / 0.93 / 0.66 |
| upconversion ⁴I₁₁/₂+⁴I₁₁/₂ | — / 0.75 / 0.33 / 0.26 |
| ⁴S₃/₂/²H₁₁/₂ cross-relaxation | 0.46 / 0.52 / 0.19 / 0.17 |

**Results.**
- All 10 values are between 0.17× and 0.93×; median error 0.31 dex (a factor of 2), with nothing
  fitted.
- The direct term matters: hopping alone gives 0.01–0.07× for the cross-relaxation, whose measured
  coefficient is concentration-independent at low N.
- The remaining shortfall is systematic (always low, worse at high N). That points to processes
  not modelled: dipole–quadrupole terms, clustering, or diffusion-limited transfer. No factor was
  tuned.
- Also updated: `validate_transfer.py` after adding phonon absorption. 9/11 published overlap values
  remain within 2× (median 0.10 dex, rms 0.22).

### 5.5 Data from open-access papers, collected by extraction agents (`gainmat/refs/data/`, `refs/oa_agent_data.py`)

Five background agents read open-access full texts (arXiv, MDPI, Sci. Rep., RSC Adv., university
repositories) and wrote one JSON record per measured value, each with citation, URL read and notes.
Committed to the repo as numbers + citations only.  `oa_agent_data.py` validates every record (glass
family, quantity, multiplet labels against the ion model, value range), drops papers already in the
hand-collected set and abstract-only entries.

| file | records | content |
|---|---|---|
| `oa_records.jsonl` | 355 (351 used) | peaks, widths, cross sections, Ω, τrad, τ, n for Er/Yb/Nd/Tm/Ho/Pr/Dy in all ten families, incl. 69 fibre-host records (Melkumov Yb, Simpson thesis Tm, Le Gouët Ho …) |
| `et_records.jsonl` | 31 | published Förster–Dexter microparameters (overlap and decay-fit) |
| `et_macro_records.jsonl` | 46 | lifetime-vs-concentration series and cross-relaxation coefficients (Taher Tm tellurite, Pugliese Er phosphate …) |
| `stark_records.jsonl` | 23 | measured Stark levels / zero-line energies in glass (Faure Tm, Melkumov Yb, Peretti/Staudt Er) |

Access limits: Optica pages (bot check), PMC captchas and a shared web-search budget capped the
yield; Optica-journal papers are under-represented.

### 5.6 Energy-transfer microparameters on 44 published values (`scripts/validate_transfer.py`)

Hand-collected 13 + agent 31 values, all predicted from the model's own spectra (model defaults):
- overlap integrals of measured spectra, n = 36: median error 0.21 dex, rms 0.36, 69% within 2×
  (unflagged subset similar);
- microparameters fitted to measured decays, n = 8: median 0.56 dex, 38% within 2×; real glasses
  quench faster than the dipole–dipole floor (clustering, multipolar terms).
`predict()` now quotes both numbers with the transfer parameters.

### 5.7 Stark structure vs measured levels in glass (`scripts/validate_stark_agent.py`)

Pure prediction (anchors off) vs the agents' glass records: zero-line energies E0 of the lowest
excited level, n = 8: median |error| 65 cm⁻¹ (Er ⁴I₁₃/₂ in Al/Ge-silica +15/+19 cm⁻¹; Tm ³H₄ in
Al/P/Ge-silica −73/−42/−128 cm⁻¹; Yb ²F₅/₂ +62/+92 cm⁻¹; Tm ³F₄ germanate +69 cm⁻¹).  Total
splittings, n = 7: median model/measured 0.91, 71% within 25% (Yb ²F₇/₂ Al-silica 0.69 against a
1210 cm⁻¹ Gaussian-decomposition value; Er chalcogenide 0.5 against a three-band deconvolution).
The Tm ³H₄ level sits 40–130 cm⁻¹ too low in silicates — the same direction as the Tm band-position
bias already flagged in `predict()`.

### 5.8 New capability: lifetime vs concentration (`gainmat/forward/quenching.py`, `scripts/validate_quenching.py`)

Levels with a cross-relaxation channel (Tm ³H₄, Nd ⁴F₃/₂, Pr ¹D₂, Ho ⁵S₂) quench with dopant
concentration N.  Prediction with nothing fitted: decay I(t) = exp(−t/τ0 − γ√t − Wt), static
Inokuti–Hirayama γ = (4/3)π^{3/2} N √C_CR plus Burshtein hopping W = π(2π/3)^{5/2} √(C_CR C_mig) N²,
with C_CR and C_mig from the model's spectra (§5.6) and τ0 the model's own lifetime; reported τ is
∫I dt.

Test: Tm ³H₄ in 75TeO₂–20ZnO–5Na₂O (Taher, PhD thesis, Swansea 2011; six glasses 0.82–9.06e20 cm⁻³,
τ falling 317 → 12 µs).

| N (1e20 cm⁻³) | 0.82 | 1.65 | 2.28 | 4.87 | 6.84 | 9.06 |
|---|---|---|---|---|---|---|
| measured τ (µs) | 317 | 220 | 130 | 34 | 19.9 | 12 |
| model, fully predictive | 146 | 100 | 75 | 27 | 15 | 9 |
| model, τ0 = most dilute measured | 211 | 130 | 91 | 29 | 15.7 | 9.2 |

- Median error 0.18 dex fully predictive, 0.13 dex given τ0.  At ≥ 4.9e20 cm⁻³ the predicted quench
  rate is 1.16–1.30× the measured one — the cross-relaxation physics is right to ~25%.
- At low N the model quenches ~2× too fast (part of this is using the slightly quenched 0.82e20
  sample as τ0); the fully predictive τ0 (199 µs vs ~330 µs) is the larger error.
- Outside the model, shown for honesty: Er ⁴I₁₃/₂ in phosphate (Pugliese 2015, 7.0 → 2.1 ms over
  0.02–13.7e20 cm⁻³) and Tm ³F₄ in the same tellurite (2.5 → 0.07 ms) quench by migration to OH⁻ and
  defects; the model, which has no impurity concentration, predicts them flat.
- `predict(ion, family, N_cm3=...)` now returns this self-quenched lifetime.

Also: the ion models are now cached on disk (`ions.get_ion`): building the 4f⁹ (Dy) model peaked at
2.8 GB and killed the benchmark twice; the cached model peaks at 0.45 GB.

### 5.9 External benchmark on 440 values, and the data library as a model input (`forward_results/oa_benchmark.md`)

With the agent records the open-access benchmark grows from 246 to 440 values (99 hand + 341 agent,
~90 papers, all seven ions, ten families).  Median absolute errors; "lookup" = mean of the same
quantity for the same ion/transition/family in OTHER papers; "library Ω" = the model given the
geometric-mean Ω of the same ion and family from other papers (still leave-one-paper-out).

| quantity | model (own Ω) | model + library Ω | model + paper's Ω | lookup |
|---|---|---|---|---|
| absorption peak (nm) | 2.9 | 3.7* | 4.6 (n=14) | 2.0 |
| emission peak (nm) | 12.0 | 12.2* | 12.4 (n=31) | 13.5 |
| emission FWHM (nm) | 29 | 28* | 19 (n=10) | 24 |
| σ_abs peak (dex) | 0.112 | 0.115* | — | 0.173 |
| σ_em peak (dex) | 0.188 | **0.144*** | 0.150 | 0.072 |
| τ_rad (dex) | 0.135 | **0.085*** | **0.047** (n=51) | 0.122 |
| τ measured (dex) | 0.318 | 0.266* | 0.20 | 0.167 |
| refractive index | 0.070 | — | — | 0.131 |
| Ω (dex) | 0.119 | (= lookup) | — | 0.074 |

\* on the rows where a lookup exists.

- The larger set is a harder test than the first 246 values: it adds many Tm fibre absorption peaks
  (the known Tm ³F₄ band-position bias; Simpson, Chen, Cajzl: model 1710–1724 nm vs measured
  1577–1650 nm) and Pr chalcogenide bands (+75–170 nm), which move the absorption-peak median
  from 1.8 to 2.9 nm.
- The model alone beats lookup on emission peak, σ_abs and n; lookup wins on σ_em and Ω.
- **Adopted:** physics + data.  Feeding the model the library Ω makes it beat lookup on τ_rad
  (0.085 vs 0.122 dex) and cuts σ_em and τ errors by 25% and 16%.  `predict()` now uses the
  library Ω by default when the library has the ion in that family (provenance says how many papers),
  the user's own Ω when given (τ_rad then 0.047 dex), and the old defaults otherwise.
- Measured lifetimes stay the weakest quantity (0.27–0.32 dex vs 0.17 for lookup): they carry each
  sample's concentration and OH quenching.  §5.8 covers the cross-relaxation part when N is known;
  OH/defect quenching is not modelled.
