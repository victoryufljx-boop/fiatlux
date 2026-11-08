# gainmat

First-principles optical properties of laser and fiber-amplifier gain media, built to
generate training data for an ML layer later. See `SPEC.md` for the full input/output
definition.

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
