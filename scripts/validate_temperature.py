"""Temperature dependence of the Tm3+ 3F4->3H6 emission (Turri et al., JAP 103, 093104 (2008), Sec. III.C):
linear fits over 100-300 K of the peak stimulated-emission cross section and of the peak wavelength,
germanate / silica (10.3 wt% Al2O3) / borophosphate.  The model is run at 100, 200, 300 K with each glass's
measured Omegas (Turri Table I) and refractive index of the family; compared: relative slope of the peak
cross section (%/K) and peak-wavelength slope (nm/K).

    python scripts/validate_temperature.py -> forward_results/temperature_validation.json
"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from gainmat.forward import glass as g
g.load_refs(); g.apply_site_refs(); g.calibrate_spm(370.0, g.CARNALL_ER_LAF3_BKQ, verbose=False)
TURRI = {  # (family, Omegas, sigma(T) = a + b T [cm2], lam(T) = c + d T [nm])
    "germanate": ("germanate", (4.4, 1.2, 0.9), (5.9e-21, -6.2e-24), (1870.0, 0.0)),
    "silica": ("aluminosilicate", (3.7, 2.3, 0.6), (4.3e-21, -2.3e-24), (1871.0, -0.21)),
    "phosphate": ("phosphate", (3.6, 0.7, 1.0), (3.9e-21, -2.6e-24), (1857.0, -0.07)),
}
Ts = (100.0, 200.0, 300.0); rows = []
for name, (fam, om, (a, b), (c, d)) in TURRI.items():
    sig, lam = [], []
    for T in Ts:
        sp = g.run("Tm3+", fam, pumps=[], n_sites=60, overrides=dict(omegas=om, T=T))["_spec"]
        y = sp["sigma_em"][("3F4", "3H6")]; nu = sp["nu"]; m = (1e7 / nu > 1600) & (1e7 / nu < 2150)
        i = np.argmax(np.where(m, y, 0)); sig.append(float(y[i])); lam.append(float(1e7 / nu[i]))
    ps = np.polyfit(Ts, sig, 1); pl = np.polyfit(Ts, lam, 1)
    r = dict(glass=name, family=fam, model_sigma=sig, model_lam=lam,
             model_rel_slope_pct_per_K=float(100 * ps[0] / np.polyval(ps, 300)), meas_rel_slope_pct_per_K=float(100 * b / (a + b * 300)),
             model_sigma_300K=float(np.polyval(ps, 300)), meas_sigma_300K=a + b * 300,
             model_lam_slope_nm_per_K=float(pl[0]), meas_lam_slope_nm_per_K=d, model_lam_300K=float(np.polyval(pl, 300)), meas_lam_300K=c + d * 300)
    rows.append(r)
    print(f"{name:10s} sigma slope model {r['model_rel_slope_pct_per_K']:+.3f} %/K meas {r['meas_rel_slope_pct_per_K']:+.3f} | sigma300 model {r['model_sigma_300K']:.2e} meas {r['meas_sigma_300K']:.2e}"
          f" | lam slope model {r['model_lam_slope_nm_per_K']:+.3f} meas {d:+.3f} nm/K | lam300 model {r['model_lam_300K']:.0f} meas {c+d*300:.0f}", flush=True)
json.dump(dict(src="Turri et al., J. Appl. Phys. 103, 093104 (2008), Sec. III.C fits and Table I", rows=rows),
          open("forward_results/temperature_validation.json", "w"), indent=1, default=float)
