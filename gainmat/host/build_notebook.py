"""Build notebooks/01_host_properties.ipynb programmatically (nbformat) and validate it.
Run:  python gainmat/host/build_notebook.py
"""
import os
import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

cells.append(md("""# gainmat 01 — host properties (band gap, n, phonon cutoff)

Runs the host tier for the validation set on a Colab **A100** and compares with experiment.

* **Band gap**: GPAW plane-wave GLLB-SC (with derivative discontinuity) is the default; PBE is reported alongside. TB09 (mBJ) and HSE06 are selectable with `GAP_METHOD`. In the bake-off (`BAKEOFF_GAPS.md`) GLLB-SC had MAE ≈ 0.4–0.6 eV against the *fundamental* gaps of LiF, MgO, CaF2, Al2O3 and quartz, tending to overshoot. TB09 and HSE06@PBE came in 1–3 eV low and cost more.
* **ε∞ / n**: GPAW RPA (optional, off by default: expensive, and it runs on CPU only).
* **Max phonon**: MACE-MP-0 *medium* + phonopy (GPU). No LO–TO correction, so the max is TO-like. Expect ±20–60 cm⁻¹ against experiment.
* **Glasses**: MACE melt-quench MD (GPU) → VDOS cutoff.

**Setup**: upload the `gainmat` project folder (the one containing `gainmat/`, `notebooks/`, `tests/`), or a zip of it, to `/content`. The project then lives at `/content/gainmat`.

**GPAW runs under MPI.** The notebook kernel is one serial rank (MPI is compiled in but unused, OpenMP is off). So `run_host(..., gpaw_nproc=GPAW_NPROC)` writes a small per-host driver script and runs it with `mpiexec -n <ncpu> --allow-run-as-root python ...`. It sets `GPAW_MPI_BACKEND=cgpaw`, which GPAW ≥ 25 requires under a launcher, and collects the JSON. The code is in `gainmat/host/mpi_run.py`. MACE uses the GPU in the kernel as before."""))

cells.append(md("## 1. Install (≈5–8 min; GPAW compiles from source)"))
cells.append(code("""%%bash
set -e
apt-get -qq update
apt-get -qq install -y libxc-dev libopenblas-dev libfftw3-dev openmpi-bin libopenmpi-dev > /dev/null
mkdir -p ~/.gpaw
cat > ~/.gpaw/siteconfig.py <<'EOF'
libraries = ['xc', 'openblas', 'fftw3']
fftw = True
EOF
pip -q install ase phonopy spglib seekpath nbformat
pip -q install gpaw          # builds _gpaw against libxc/openblas/fftw; pulls gpaw-data (PAW setups)
pip -q install mace-torch    # uses Colab's preinstalled CUDA torch
gpaw info | grep -E 'MPI|FFTW|libxc|scalapack'
# MPI check: must print the number of ranks
GPAW_MPI_BACKEND=cgpaw mpiexec -n $(nproc) --allow-run-as-root --oversubscribe \
    python -c "from gpaw.mpi import world; world.rank == 0 and print('GPAW MPI ranks:', world.size)"""))

cells.append(code("""# If you uploaded a zip instead of the folder:
import os, zipfile, glob
if not os.path.isdir('/content/gainmat/gainmat'):
    zips = glob.glob('/content/gainmat*.zip')
    if zips:
        zipfile.ZipFile(zips[0]).extractall('/content')
assert os.path.isdir('/content/gainmat/gainmat'), 'upload the gainmat folder to /content first'
import sys; sys.path.insert(0, '/content/gainmat')
import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available(),
                    torch.cuda.get_device_name(0) if torch.cuda.is_available() else '')"""))

cells.append(md("## 2. Structure checks + smoke tests"))
cells.append(code("""import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', str(os.cpu_count()))
from gainmat.structures import STRUCTURES, verify_structure
for k in STRUCTURES: print(verify_structure(k))
!cd /content/gainmat && python -m pytest -q tests/test_host.py -x 2>&1 | tail -5"""))

cells.append(md("""## 3. Settings
`LEVEL="full"` uses path band edges, ecut 500 eV, and 10 Å phonon supercells. GPAW runs under `mpiexec -n GPAW_NPROC`, all vCPUs by default. The 80-atom YAG primitive cell (PBE + GLLB-SC, Γ-centred 2×2×2, ecut 500) is the long pole: about **YAG_EST** on a 12-vCPU Colab runtime, compared with several hours serially in the kernel (see `BAKEOFF_GAPS.md`). Set `SKIP_GAP` to leave it out."""))
cells.append(code("""LEVEL = 'full'
CRYSTALS = ['CaF2', 'LiYF4', 'Al2O3', 'SiO2-alpha-quartz', 'Y3Al5O12']
GLASSES = ['SiO2-glass', 'ZBLAN']
DO_EPS = False            # RPA eps_inf (+ scissor to the GLLB-SC gap); CPU, runs under MPI too
GAP_METHOD = 'gllbsc'     # 'gllbsc' (default, bake-off winner) | 'pbe' | 'tb09' | 'hse06'
GPAW_NPROC = os.cpu_count()   # MPI ranks for GPAW (1 = serial in-kernel, slow)
SKIP_GAP = []             # e.g. ['Y3Al5O12'] if short on time
MACE_MODEL = 'medium'     # MACE-MP-0 medium; alternatives: 'medium-mpa-0', 'medium-omat-0'
RESULTS = '/content/results/host_results.json'
os.makedirs('/content/results', exist_ok=True)"""))

cells.append(md("## 4. Run crystals"))
cells.append(code("""import json, time
from dataclasses import asdict
from gainmat.host.run_host import run_host, compare_to_experiment

results, rows = {}, []
def save():
    with open(RESULTS, 'w') as f:
        json.dump(dict(results=results, comparison=rows), f, indent=1, default=float)

for key in CRYSTALS:
    t0 = time.time()
    hp = run_host(key, level=LEVEL, do_gap=key not in SKIP_GAP, do_eps=DO_EPS,
                  mace_model=MACE_MODEL, gap_method=GAP_METHOD, gpaw_nproc=GPAW_NPROC)
    d = asdict(hp); d['method'].pop('tracebacks', None)
    results[key] = d
    rows += compare_to_experiment(key, hp)
    print(f'{key}: {time.time()-t0:.0f}s  PBE={hp.gap_pbe_eV}  gap[{hp.gap_method}]={hp.gap_eV}  n={hp.n_inf}  '
          f'wmax={hp.phonon_max_cm1}  weff={hp.phonon_eff_cm1}  errors={hp.method.get(\"errors\")}')
    save()"""))

cells.append(md("## 5. Run glasses (MACE melt-quench MD on the GPU)"))
cells.append(code("""for key in GLASSES:
    t0 = time.time()
    hp = run_host(key, level=LEVEL, mace_model=MACE_MODEL)
    d = asdict(hp); d['method'].pop('tracebacks', None)
    results[key] = d
    rows += compare_to_experiment(key, hp)
    print(f'{key}: {time.time()-t0:.0f}s  w_eff={hp.phonon_eff_cm1}  top peak={hp.phonon_max_cm1}  '
          f'errors={hp.method.get(\"errors\")}')
    save()"""))

cells.append(md("## 6. Computed vs experiment"))
cells.append(code("""import pandas as pd
df = pd.DataFrame(rows)
pd.set_option('display.max_colwidth', 60)
show = df[['host','quantity','computed','experiment','err_pct','exp_kind','ref_verified']].copy()
show['computed'] = show['computed'].map(lambda x: None if x is None else round(x, 3))
show['err_pct'] = show['err_pct'].map(lambda x: None if x is None else round(x, 1))
display(show)
save(); print('saved', RESULTS)"""))

cells.append(code("""# Phonon DOS plots
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(7, 4))
for key, d in results.items():
    dos = d['method'].get('phonon_dos') or d['method'].get('vdos')
    if dos:
        import numpy as np
        g = np.array(dos['dos']); ax.plot(dos['freq_cm1'], g / g.max(), label=key)
ax.set_xlabel('frequency (cm$^{-1}$)'); ax.set_ylabel('DOS (norm.)'); ax.set_xlim(0, 1400); ax.legend()
plt.savefig('/content/results/host_phonon_dos.png', dpi=150)"""))

cells.append(md("""## Notes on accuracy
* GLLB-SC is compared against the *fundamental* gap. Wide-gap fluorides have exciton peaks about 1 eV below it, and the `vs_exciton` row shows that gap.
* MACE-MP-0 tends to soften phonons (PES softening). Expect top optical modes up to ~5–10 % low; `medium-mpa-0` or `medium-omat-0` usually do better.
* `phonon_max_cm1` has no LO–TO correction, so compare it with the TO / Raman top mode rather than the LO value. The comparison row uses the LO reference and so understates the result.
* Glass VDOS from ps-scale quenches is semi-quantitative."""))

nb = nbf.v4.new_notebook(cells=cells, metadata={
    "kernelspec": {"name": "python3", "display_name": "Python 3"},
    "accelerator": "GPU", "colab": {"gpuType": "A100", "provenance": []}})
nbf.validate(nb)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "notebooks", "01_host_properties.ipynb")
nbf.write(nb, out)
print("wrote", out)
