import sys, json, numpy as np
sys.path.insert(0, "/home/claude/mainrepo")
from gainmat.forward.glass import load_refs, apply_site_refs, calibrate_spm, CARNALL_ER_LAF3_BKQ, run, load_md_sites
from gainmat.forward.spectra import fwhm_nm
load_refs(); apply_site_refs(); calibrate_spm(370.0, CARNALL_ER_LAF3_BKQ, verbose=False)
S = sys.argv[1]; s = 1.175; out = {}
for g, fam in (("Er-Al3-SiO2", "aluminosilicate"), ("Er-Al9-SiO2", "aluminosilicate"), ("Er-Ge10-SiO2", "germanosilicate"), ("Er-SiO2", "silica")):
    st = load_md_sites(f"{S}/{g}.json")
    # split by replica-independent strided blocks: 4 disjoint subsets of 40
    rng = np.random.default_rng(1); idx = rng.permutation(len(st))
    pk = []
    for b in range(4):
        sub = [st[i] for i in idx[b*40:(b+1)*40]]
        r = run("Tm3+", fam, n_sites=40, pumps=[], overrides=dict(sites=sub, cf_scale=s))
        sp = r["_spec"]; lam = 1e7/sp["nu"]; y = sp["sigma_abs"][("3F4","3H6")]; m = (lam>1450)&(lam<1950)
        pk.append(float(lam[m][np.argmax(y[m])]))
    out[g] = dict(peaks=pk, mean=float(np.mean(pk)), sem=float(np.std(pk, ddof=1)/2))
    print(g, out[g], flush=True)
json.dump(out, open(f"{S}/tm_noise.json", "w"), indent=1)
