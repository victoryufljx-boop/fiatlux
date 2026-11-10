"""Grid of Tm 3F4 absorption peak (and Er/Yb Stark patterns on pure silica) against separate
superposition strengths for bridging (a) and non-bridging (b) oxygen, on MD sites.

    python scripts/md_nbo_grid.py GLASS FAMILY OUT.json [--stark]
"""
import sys, json, time, numpy as np
sys.path.insert(0, __import__("os").path.abspath(__import__("os").path.join(__import__("os").path.dirname(__file__), "..")))
from gainmat.forward.glass import load_refs, apply_site_refs, calibrate_spm, CARNALL_ER_LAF3_BKQ, run, load_md_sites_tagged, spm_for
from gainmat.forward.spectra import spm_bkq, fwhm_nm
from gainmat.forward.ions import get_ion
load_refs(); apply_site_refs(); calibrate_spm(370.0, CARNALL_ER_LAF3_BKQ, verbose=False)
glass, fam, out = sys.argv[1:4]; stark = "--stark" in sys.argv
A = [round(0.6 + 0.2 * i, 1) for i in range(11)]; B = [round(1.0 + 0.4 * i, 1) for i in range(10)]   # a 0.6-2.6, b 1.0-4.6
lab = lambda k, f: "Onb" if k == "NBO" else "O"
st = load_md_sites_tagged(f"md_sites/{glass}.json", lab)
sub = [st[i] for i in np.random.default_rng(1).permutation(len(st))[:60]]


def pattern(ion, label, a, b):
    m = get_ion(ion); p = dict(spm_for(ion, "O")); res = []
    pa = dict(p, **{f"B{k}": p[f"B{k}"] * a for k in (2, 4, 6)}); pb = dict(p, **{f"B{k}": p[f"B{k}"] * b for k in (2, 4, 6)})
    for s in sub:
        E, W, par = m.stark(spm_bkq(s, {"O": pa, "Onb": pb}))
        e = np.sort(E[par == m.mult_index(label)]); res.append((e - e[0])[::2])
    return np.mean(res, axis=0).round(1).tolist()


rows = []; t0 = time.time()
for a in A:
    for b in B:
        r = dict(a=a, b=b)
        if stark:
            r.update(er15=pattern("Er3+", "4I15/2", a, b), er13=pattern("Er3+", "4I13/2", a, b),
                     yb7=pattern("Yb3+", "2F7/2", a, b), yb5=pattern("Yb3+", "2F5/2", a, b))
        else:
            sp = run("Tm3+", fam, n_sites=60, pumps=[], overrides=dict(sites=sub, cf_scale=a, ligand_factors={"Onb": b / a}))["_spec"]
            lam = 1e7 / sp["nu"]; y = sp["sigma_abs"][("3F4", "3H6")]; m = (lam > 1450) & (lam < 1950)
            r.update(tm_abs=float(lam[m][np.argmax(y[m])]), tm_fwhm=float(fwhm_nm(sp["nu"], np.where(m, y, 0.0))))
        rows.append(r); print(glass, r if not stark else (a, b, r["er15"][-1], r["yb7"][-1]), round(time.time() - t0), flush=True)
json.dump(dict(glass=glass, family=fam, A=A, B=B, rows=rows), open(out, "w"), indent=1)
