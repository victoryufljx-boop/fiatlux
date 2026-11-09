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
