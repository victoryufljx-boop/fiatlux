"""Builds 03_embedded_clusters.ipynb (run: python build_03.py inside notebooks/)."""
import json, nbformat as nbf
nb = nbf.v4.new_notebook()
md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
nb.cells = [
md("""# 03 - Embedded-cluster CASSCF/NEVPT2: Ce3+:LiYF4 (4f-5d) and Ti3+:Al2O3 (3d-3d)

Tier `gainmat.cluster`. Dopant + first anion shell in QM; next cation shell as ECP-only caps
(Li+, Al3+) or, where no full-core ECP exists (Y3+), as **frozen-ion AIMP-like model potentials**
(nuclear + ECP + exact Coulomb/exchange of the frozen Y3+ density + core-projection shift, no basis
functions on Y; Seijo & Barandiaran style); ~20 A sphere of Ewald-fitted point charges.

* **Ti3+:Al2O3**: ROHF -> SA-CASSCF(1,5) -> SC-NEVPT2, vertical 2T2g -> 2Eg at the relaxed ground-state
  geometry; 2Eg Jahn-Teller minimum on the (a1g breathing, e_g theta) grid from MOM Delta-SCF, NEVPT2
  re-run there -> emission and Stokes shift.
* **Ce3+:LiYF4**: ROHF -> SA-CASSCF(1,12) (4f+5d) -> SC-NEVPT2 -> 12x12 ab initio ligand-field matrix +
  spin-orbit. First F shell relaxed in the ground state (UHF) and lowest 5d (MOM Delta-SCF) along the S4
  breathing coordinates -> relaxed emission and Stokes shift. Optional check of Ce 5s5p semicore
  correlation in NEVPT2.

`LEVEL='full'` uses Ce def2-TZVP (ECP28MWB, bundled), F def2-TZVPD, Ti/O def2-TZVP(D), 20 A PC sphere.
Flags below switch each improvement off for comparison. Expected runtimes on a 12-core Colab
high-RAM CPU runtime: Ti `full` ~1-2 h (9-point JT grid + 2 NEVPT2), Ce `full` ~4-8 h (9-point
2-coordinate scan + SA-CASSCF/NEVPT2 + semicore NEVPT2); `fast+` ~15-30 min total. gpu4pyscf is not used."""),
code("""# Install. Note: pyscf already includes mcscf, mrpt (NEVPT2) and qmmm
!pip -q install pyscf ase nbformat matplotlib"""),
code("""# Upload gainmat.zip (the repo folder) to /content, or mount Drive and copy it there
import os, sys, zipfile, glob
if not os.path.isdir('/content/gainmat') and os.path.exists('/content/gainmat.zip'):
    zipfile.ZipFile('/content/gainmat.zip').extractall('/content')
root = '/content/gainmat' if os.path.isdir('/content/gainmat/gainmat') else '/content'
sys.path.insert(0, root)
import gainmat; print('gainmat from', gainmat.__file__)"""),
code("""# Threads + level + improvement flags
import os, multiprocessing
os.environ['OMP_NUM_THREADS'] = str(multiprocessing.cpu_count())
from pyscf import lib
lib.num_threads(multiprocessing.cpu_count())
print('threads', lib.num_threads())
LEVEL = 'full'       # 'fast' = original baseline, 'fast+' = improvements with small basis, 'full'
# Flags (None = level default). Set e.g. CAP_MODE='pc', RELAX=None to reproduce the old scheme.
CAP_MODE = None      # 'aimp' | 'partial' | 'pc'        (Y3+ next shell, Ce only)
R_CAP    = None      # 5.3 A includes all 8 nearest Y3+; 4.0 only 4
CE_RELAX = None      # 'breathing2' | 'breathing' | None
TI_RELAX = None      # 'jt' | None
SEMICORE = None      # True: extra NEVPT2 with Ce 5s5p frozen
# os.environ['CE_FULL_BASIS'] = 'bundled:ce_def2-svp.nw'   # cheaper Ce basis
# os.environ['F_FULL_BASIS'] = 'def2-svpd'
CE_FLAGS = dict(cap_mode=CAP_MODE, r_cap=R_CAP, relax=CE_RELAX, semicore=SEMICORE)
TI_FLAGS = dict(relax=TI_RELAX)"""),
code("""# Sanity: embedding geometry + Madelung potential vs Ewald
from gainmat.cluster import embedding as emb
rpc = 20.0 if LEVEL == 'full' else 12.0
for args, kw in [(('Al2O3','Al','Ti',3), dict(r_qm=2.2, r_cap=3.6)),
                 (('LiYF4','Y','Ce',3), dict(r_qm=2.6, r_cap=R_CAP or 5.3, cap_mode=CAP_MODE or 'aimp'))]:
    cl = emb.make_cluster(*args, r_pc=rpc, **kw)
    s = cl.summary(); print({k: s[k] for k in ['name','n_qm','n_cap','n_aimp','n_pc','dopant_ligand_A','total_charge','madelung_ewald_V','madelung_cluster_V','fit_rms_V','notes']})"""),
code("""import time, json
from dataclasses import asdict
from gainmat.cluster.run_cluster import run_center, print_table
results, raw = {}, {}
for name, fl in [('Ti:Al2O3', TI_FLAGS), ('Ce:LiYF4', CE_FLAGS)]:
    t = time.time(); r = {}
    dp = run_center(name, LEVEL, raw_out=r, **fl)
    results[name] = asdict(dp); raw[name] = r
    print_table(name, dp); print(f'  wall {time.time()-t:.0f} s')"""),
code("""# Computed-vs-experiment table (pandas): absorption, emission, Stokes shift
import pandas as pd
from gainmat.cluster.run_cluster import compare_table
from gainmat.schema import DopantProps
rows = []
for name, d in results.items():
    dp = DopantProps(**d)
    for q, c, e, u in compare_table(name, dp):
        st = q == 'Stokes shift'
        rows.append(dict(center=name, quantity=q, calc_cm1=c, exp_cm1=e,
                         calc_nm=(1e7/c if c and not st else None), exp_nm=(None if st else 1e7/e),
                         diff_cm1=(c-e if c else None)))
df = pd.DataFrame(rows); df.round(1)"""),
code("""# Relaxation and semicore details
for name in results:
    m = results[name]['method']
    print(name, 'relax:', {k: v for k, v in (m.get('relax') or {}).items() if k != 'scan'})
    if m.get('semicore'): print('  semicore 5s5p:', m['semicore'])"""),
code("""# Stick spectra: computed (up) vs experiment (down)
import matplotlib.pyplot as plt
from gainmat.refs import cluster_refs as R
fig, axs = plt.subplots(1, 2, figsize=(11, 3.6))
ce = results['Ce:LiYF4']; ti = results['Ti:Al2O3']
ax = axs[0]
ax.vlines([1e7/x for x in R.CE_LIYF4['abs_5d_cm1']], 0, -1, color='k', label='exp. 4f-5d')
ax.vlines(ce['abs_peaks_nm'] or [], 0, 1, color='C0', label='calc. NEVPT2+SO')
ax.vlines(R.CE_LIYF4['em_nm'], 0, -0.5, color='C3', label='exp. emission')
ax.vlines([ce['em_peak_nm']] if ce['em_peak_nm'] else [], 0, 0.5, color='C1', label='calc. emission')
ax.axhline(0, color='grey', lw=.5); ax.set_xlabel('wavelength (nm)'); ax.set_title('Ce3+:LiYF4'); ax.legend(fontsize=7); ax.set_yticks([])
ax = axs[1]
ax.vlines(R.TI_AL2O3['abs_nm'], 0, -1, color='k', label='exp. abs (main, JT shoulder)')
ax.vlines([R.TI_AL2O3['em_nm']], 0, -0.5, color='C3', label='exp. emission')
ax.vlines(ti['abs_peaks_nm'] or [], 0, 1, color='C0', label='calc. 2T2g->2Eg')
ax.vlines([ti['em_peak_nm']] if ti['em_peak_nm'] else [], 0, 0.5, color='C1', label='calc. emission (JT min)')
ax.axhline(0, color='grey', lw=.5); ax.set_xlabel('wavelength (nm)'); ax.set_title('Ti3+:Al2O3'); ax.legend(fontsize=7); ax.set_yticks([])
plt.tight_layout(); os.makedirs('/content/results', exist_ok=True); plt.savefig('/content/results/cluster_sticks.png', dpi=150); plt.show()"""),
code("""os.makedirs('/content/results', exist_ok=True)
with open('/content/results/cluster_results.json', 'w') as f:
    json.dump(dict(level=LEVEL, results=results, raw=raw), f, indent=1,
              default=lambda x: x.tolist() if hasattr(x, 'tolist') else str(x))
print('saved /content/results/cluster_results.json')"""),
]
nbf.validate(nb)
nbf.write(nb, '03_embedded_clusters.ipynb')
print('ok')
