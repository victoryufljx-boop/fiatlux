"""Build notebooks/01_host_properties.ipynb programmatically (nbformat) and validate it.
Run:  python gainmat/host/build_notebook.py
"""
import os
import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

cells.append(md("""# gainmat 01 — host properties (band gap, n, phonon cutoff)

Runs the host tier for the validation set on a Colab **A100** and compares with experiment.

* **Band gap**: GPAW plane-wave PBE and GLLB-SC (with derivative discontinuity). Expect PBE −30…45 %, GLLB-SC within ~10–15 % of the *fundamental* gap.
* **ε∞ / n**: GPAW RPA (optional, off by default: expensive, and it runs on CPU only).
* **Max phonon**: MACE-OMAT-0 *medium* + phonopy (GPU) at the experimental lattice (positions relaxed), with the LO–TO non-analytic correction where Born charges + ε∞ are tabulated (CaF2, Al2O3), so the max includes LO. Bake-off: Γ TO/Raman MAE ≈ 8 cm⁻¹ (see BAKEOFF_PHONONS.md). Legacy MACE-MP-0 protocol: `MACE_MODEL='medium'` + `relax_cell=True, nac=False`.
* **Glasses**: MACE melt-quench MD (GPU) → VDOS cutoff.

**Setup**: upload the `gainmat` project folder (the one containing `gainmat/`, `notebooks/`, `tests/`), or a zip of it, to `/content`. The project then lives at `/content/gainmat`.

Note: GPAW runs on the CPU in serial mode inside the notebook kernel (BLAS threads via `OPENBLAS_NUM_THREADS`). MACE uses the GPU. For faster GPAW, run `gainmat/host/run_host.py` logic under `mpiexec -n 12 python ...` from a `%%bash` cell."""))

cells.append(md("## 1. Install (≈5–8 min; GPAW compiles from source)"))
cells.append(code("""%%bash
set -e
apt-get -qq update
apt-get -qq install -y libxc-dev libopenblas-dev libfftw3-dev > /dev/null
mkdir -p ~/.gpaw
cat > ~/.gpaw/siteconfig.py <<'EOF'
libraries = ['xc', 'openblas', 'fftw3']
fftw = True
EOF
pip -q install ase phonopy spglib seekpath nbformat
pip -q install gpaw          # builds _gpaw against libxc/openblas/fftw; pulls gpaw-data (PAW setups)
pip -q install mace-torch 'e3nn==0.4.4'   # uses Colab's preinstalled CUDA torch; MACE-OMAT-0 downloads on first use
# optional bake-off models (separate envs recommended: sevenn needs e3nn>=0.5, breaks MACE pickles):
# pip -q install mattersim     |   pip -q install sevenn
gpaw info | head -25"""))

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
`LEVEL="full"` uses path band edges, ecut 500 eV, and 10 Å phonon supercells. On Colab, GPAW for the 80-atom YAG primitive cell with GLLB-SC takes the longest (tens of minutes). Set `SKIP_GAP` to leave it out."""))
cells.append(code("""LEVEL = 'full'
CRYSTALS = ['CaF2', 'LiYF4', 'Al2O3', 'SiO2-alpha-quartz', 'Y3Al5O12']
GLASSES = ['SiO2-glass', 'ZBLAN']
DO_EPS = False            # RPA eps_inf: slow; set True to try (CPU)
SKIP_GAP = []             # e.g. ['Y3Al5O12'] if short on time
MACE_MODEL = 'medium-omat-0'   # bake-off winner; legacy: 'medium' (MACE-MP-0); also 'medium-mpa-0'
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
                  mace_model=MACE_MODEL)
    d = asdict(hp); d['method'].pop('tracebacks', None)
    results[key] = d
    rows += compare_to_experiment(key, hp)
    print(f'{key}: {time.time()-t0:.0f}s  PBE={hp.gap_pbe_eV}  GLLBSC={hp.gap_gllbsc_eV}  '
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
