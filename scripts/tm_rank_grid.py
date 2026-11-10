import sys, numpy as np, json, itertools
sys.path.insert(0, "/home/claude/mainrepo")
from gainmat.forward import glass as g
from gainmat.forward.spectra import fwhm_nm
g.load_refs(); g.apply_site_refs(); g.calibrate_spm(370.0, g.CARNALL_ER_LAF3_BKQ, verbose=False)
fam = sys.argv[1]; r2s = [float(x) for x in sys.argv[2].split(",")]; out = sys.argv[3]
R4 = [0.6, 1.0, 1.5, 2.2, 3.0]; R6 = [0.5, 1.0, 2.0, 3.0]
def band(sp, pair, which, lo, hi):
    lam = 1e7/sp["nu"]; y = sp["sigma_abs" if which=="abs" else "sigma_em"][pair]; m=(lam>lo)&(lam<hi)
    return float(lam[m][np.argmax(y[m])]), float(fwhm_nm(sp["nu"], np.where(m,y,0))), y
rows = []
for r2, r4, r6 in itertools.product(r2s, R4, R6):
    g.ION_LIG_RANK["Tm:O"] = (r2, r4, r6)
    sp = g.run("Tm3+", fam, pumps=[], n_sites=40)["_spec"]
    a = band(sp,("3F4","3H6"),"abs",1450,1950); e = band(sp,("3F4","3H6"),"em",1600,2150)
    h4 = band(sp,("3H4","3H6"),"abs",740,850); h5 = band(sp,("3H5","3H6"),"abs",1100,1300)
    r805 = float(np.interp(1e7/805, sp["nu"], h4[2]) / h4[2].max())
    r = dict(r2=r2, r4=r4, r6=r6, abs=a[0], abs_fwhm=a[1], em=e[0], em_fwhm=e[1], h4=h4[0], h4_fwhm=h4[1], r805=r805, h5=h5[0])
    rows.append(r); print(json.dumps({k: round(v, 2) for k, v in r.items()}), flush=True)
json.dump(rows, open(out, "w"))
