# gainmat

First-principles optical properties of laser and fiber-amplifier gain media, built to
generate training data for an ML layer later. See `SPEC.md` for the full input/output
definition.

## Quick start: predict an ion in a glass

```python
from gainmat.forward.predict import predict
p = predict("Er3+", "aluminosilicate")                     # by glass family
p = predict("Tm3+", composition="70TeO2-15K2O-15Nb2O5")    # or by composition (mol%)
p = predict("Tm3+", "tellurite", omegas=(4.51, 0.76, 1.13))  # measured Judd-Ofelt Omegas, if you have them
```
`p` holds the main bands (peak wavelengths, peak cross sections, emission FWHM), level lifetimes and
quantum efficiencies, gain windows, energy-transfer microparameters (migration, cross-relaxation,
upconversion), an empirical uncertainty for each quantity, the provenance of every input, and the
systematic biases found in validation. Command line: `python -m gainmat.forward.predict Er3+ aluminosilicate`.

How far to trust it (held-out literature, `forward_results/oa_benchmark.md`, `transfer_validation.json`):

| quantity | typical error | best input to improve it |
|---|---|---|
| radiative lifetimes | 0.11 dex; 0.02 dex with measured Omegas | `omegas=` |
| band positions | 1-7 nm (Er, Nd, Ho, Yb); Tm 3F4 up to ~90 nm | `measured_absorption=` override in `glass.run` |
| peak cross sections | 0.19 dex; 0.14 dex with measured Omegas | `omegas=`, measured absorption |
| energy-transfer microparameters | within 2x for 9 of 11 published values | none needed; multi-phonon processes ~3x |

Switches: `GAINMAT_ML=1` turns on the learned layer (ZBLAN host broadening, silica 2 um quenching);
`GAINMAT_ANCHORS=0` turns off the measured-Stark-level mode (Yb, Er in silica fibres). Results log:
`ml_results/RESULTS.md`; validation scorecard: `forward_results/SCORECARD.md`.

## Running on Colab (A100)

1. Upload `gainmat.zip` to the Colab session (Files pane), or put it in Drive.
2. Open any notebook in `notebooks/` and run top to bottom. Each one unzips the
   package to `/content/gainmat` if needed, installs its own dependencies, and writes
   results to `/content/results/`.

| Notebook | What it computes | Hardware | Rough time |
|---|---|---|---|
| `01_host_properties.ipynb` | Band gaps, phonon cutoffs (+ glass DOS) for the validation hosts | GPU for MACE; CPU for GPAW | see notebook |
| `02_lanthanide_4f.ipynb` | 4f^N levels, Judd–Ofelt rates, multiphonon quenching | CPU, seconds | < 1 min |
| `03_embedded_clusters.ipynb` | Ce³⁺ 4f→5d and Ti³⁺ d–d transitions | CPU (many cores) | see notebook |
| `00_validation_summary.ipynb` | Collects all results into one computed-vs-experiment table | CPU | seconds |

Every experimental number lives in `gainmat/refs/*.py` with its citation.

## Layout
```
gainmat/
  schema.py        inputs/outputs shared by all tiers
  structures.py    validation crystal structures
  host/            GPAW gap, MACE phonons, glass MD
  ln4f/            4f^N CI, Judd-Ofelt, energy-gap law
  cluster/         embedded-cluster CASSCF/NEVPT2
  refs/            experimental references with citations
tests/             pytest smoke tests (CPU)
notebooks/         Colab notebooks
```

## GPAW under MPI
In a notebook kernel GPAW runs on one core. `run_host(key, gpaw_nproc=os.cpu_count())`
(or `gainmat.host.mpi_run.run_gpaw_mpi`) runs the gap / ε∞ steps through
`mpiexec -n N --allow-run-as-root python <driver>` with `GPAW_MPI_BACKEND=cgpaw` and
collects the JSON. Command line: `python -m gainmat.host.mpi_run CaF2 full 12`.
Method choice and benchmarks: `BAKEOFF_GAPS.md`.
