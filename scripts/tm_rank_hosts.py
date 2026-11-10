import sys, numpy as np, json
sys.path.insert(0, "/home/claude/mainrepo")
from gainmat.forward import glass as g
from gainmat.forward.spectra import fwhm_nm
g.load_refs(); g.apply_site_refs(); g.calibrate_spm(370.0, g.CARNALL_ER_LAF3_BKQ, verbose=False)
def band(sp, pair, which, lo, hi):
    lam = 1e7/sp["nu"]; y = sp["sigma_abs" if which=="abs" else "sigma_em"][pair]; m=(lam>lo)&(lam<hi)
    return round(float(lam[m][np.argmax(y[m])])), round(float(fwhm_nm(sp["nu"], np.where(m,y,0))))
for rk in [(1,1,1), (2,0.6,0.5), (2,1.5,1), (1,0.6,3)]:
    g.ION_LIG_RANK["Tm:O"] = rk
    for fam in ["silica", "germanosilicate", "germanate", "tellurite", "phosphate"]:
        sp = g.run("Tm3+", fam, pumps=[], n_sites=60)["_spec"]
        print(rk, fam, "3F4 abs", band(sp,("3F4","3H6"),"abs",1450,1950), "em", band(sp,("3F4","3H6"),"em",1600,2150),
              "3H4", band(sp,("3H4","3H6"),"abs",740,850), "3H5", band(sp,("3H5","3H6"),"abs",1100,1300)[0], flush=True)
