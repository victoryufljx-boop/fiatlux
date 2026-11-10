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
