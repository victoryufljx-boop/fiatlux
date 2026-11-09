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
