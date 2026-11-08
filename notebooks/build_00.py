import nbformat as nbf
nb = nbf.v4.new_notebook()
md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
nb.cells = [
md("""# gainmat 00 — validation summary
Collects `/content/results/{host,ln,cluster}_results.json` from notebooks 01–03 into one
computed-vs-experiment view. Run 01–03 first (any subset works; missing tiers are skipped).

**Suggested order on Colab:** 02 (seconds, CPU) → 01 (A100 for MACE, CPU for GPAW) → 03 (CPU-heavy;
pick the high-RAM runtime and expect hours at `LEVEL='full'`; try `'fast'` first)."""),
code("""import os, json, pandas as pd
R = '/content/results' if os.path.isdir('/content') else '../results'
def load(name):
    p = os.path.join(R, name)
    if not os.path.exists(p): print('missing', p); return None
    return json.load(open(p))
host, ln, cl = load('host_results.json'), load('ln_results.json'), load('cluster_results.json')"""),
md("## Host tier: gap, phonon cutoff"),
code("""if host:
    df = pd.DataFrame(host['comparison'])
    display(df)"""),
md("## Lanthanide tier: levels, Judd–Ofelt rates, multiphonon quenching"),
code("""if ln:
    v = ln['validation']
    print(json.dumps(v, indent=1, default=str)[:6000])"""),
md("## Embedded clusters: Ce³⁺ 4f→5d, Ti³⁺ d–d"),
code("""if cl:
    print('level:', cl.get('level'))
    for k, r in cl['results'].items():
        print(k, 'abs_peaks_nm =', r.get('abs_peaks_nm'), ' em_peak_nm =', r.get('em_peak_nm'),
              ' error =', r.get('method', {}).get('error'))"""),
md("""## Coupling the tiers
The host tier's phonon cutoff feeds the lanthanide tier's quantum efficiency, and its band gap
bounds where Ce³⁺ 5d levels can sit without photoionizing."""),
code("""import sys
for cand in ['/content/gainmat', '..']:
    if os.path.isdir(os.path.join(cand, 'gainmat')): sys.path.insert(0, cand); break
from gainmat.ln4f.run_ln import run_dopant
if host:
    hw = {k: v.get('phonon_max_cm1') for k, v in host['results'].items()}
    print('computed phonon cutoffs (cm-1):', hw)
    # (ion, ln-tier host key, host-tier results key, n)
    for ion, key, hkey, n in [('Er3+', 'LiYF4', 'LiYF4', 1.45), ('Nd3+', 'YAG', 'Y3Al5O12', 1.82)]:
        w = hw.get(hkey)
        if w:
            dp = run_dopant(ion, key, n_index=n, phonon_max_cm1=w, nr_model='phonon')
            print(ion, key, 'eta =', {k: round(x, 3) for k, x in dp.quantum_eff.items()})"""),
]
nbf.validate(nb)
nbf.write(nb, '00_validation_summary.ipynb')
print('ok')
