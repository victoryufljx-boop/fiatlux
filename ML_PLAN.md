# ML layer: plan

**Status (2026-11-09):** layers L1-L4 and UQ run; results and adoption decisions in `ml_results/RESULTS.md`, comparison against the explicit model in `forward_results_ml/COMPARISON.md`.

## The problem, stated carefully

We have an explicit forward model, composition + dopant -> (Stark levels, cross sections,
lifetimes, gain windows, device performance), and about 120 checks against measurements. It
gets the physics right in kind: band positions, lifetimes in ZBLAN, energy transfer, device
slopes. The errors that remain are systematic: band widths ~1.5x too narrow, absolute cross
sections 0.5-0.9x, Er sigma_e/sigma_a too high, Yb Stark structure, quenching of the ~2 um
levels in silica.

There is very little data per output: tens of lifetimes, tens of cross sections, a handful of
devices. A black-box model trained on those outputs would memorise them and say nothing
about a new glass. So the ML does not predict outputs. **It predicts the few physical inputs
the explicit model cannot compute from first principles, as functions of host descriptors, and
the explicit model turns them into outputs.** Every physical rule the explicit model enforces
holds by construction: McCumber reciprocity, the Judd-Ofelt sum rule, Boltzmann populations,
energy-gap-law monotonicity, rate-equation conservation.

## What is learned (latent physical quantities)

| layer | learned quantity | from | data |
|---|---|---|---|
| L1 | Judd-Ofelt Omega_2,4,6 (ion, composition) | composition descriptors with physical meaning (optical basicity, cation field strength, anion type, network former) | LOMS (~1500 records, 9 ions, hundreds of glasses) |
| L2 | host microstructure: crystal-field strength scale, site disorder (radial / angular), nephelauxetic scatter, homogeneous width | host descriptors | band widths, Stark splittings, spectra of several ions in the same host (Wetenkamp linewidths, digitised spectra, Desurvire, Pask, vendor data) |
| L3 | non-radiative channels beyond the gap law (quenching ~ concentration, codopant) | concentration, gap, host | sparse: report, constrain, do not over-fit |
| L4 | residual correction with calibrated uncertainty | small GP on (ion, transition energy, host descriptors) | residuals after L1-L2 |

L2 is the core "rule" test: one set of host parameters must explain *every* ion in that host
(Er, Tm, Ho, Nd, Yb...) at once. If parameters fitted on some ions predict the widths and
splittings of a held-out ion, the model has learned something about the glass, not about
the data points.

## Evaluation (what "grasps the rules" means here)

1. **Out-of-distribution splits**, never random splits:
   - leave-one-glass-family-out (e.g. train without tellurites, predict tellurites),
   - leave-one-ion-out (host parameters from other ions),
   - leave-one-publication-out (LOMS),
   - extrapolation in descriptor space (train on low optical basicity, test on high).
2. **Baselines** on the same splits: explicit model alone; family-mean / nearest-neighbour;
   a black-box regressor on the outputs directly. The structured model must beat them OOD.
3. **Physical consistency checks** on predictions: McCumber, JO sum rule, monotonic
   gap-law, sign of known composition trends (Omega_2 grows with covalency/asymmetry;
   fluorides have small Omega_2; Omega_6 falls with optical basicity in oxides). Learned
   trends are checked against the literature, not imposed blindly.
4. **Uncertainty**: every learned quantity carries an ensemble / GP spread, and OOD
   predictions must be wider where the data are thin (coverage of held-out points).

## Deliverables

- `gainmat/ml/` package: descriptors, observation table, L1-L4 models, evaluation harness.
- `ml_results/`: OOD scorecards per layer, baselines, consistency checks, plots.
- Forward model switch: `use_ml=True` takes Omega and host microstructure from the learned
  layers; the full validation suite is rerun with and without.
