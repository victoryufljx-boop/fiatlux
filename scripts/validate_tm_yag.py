"""Tm3+:YAG (Gruber et al., PRB 40, 9464 (1989)) as a test of the crystal-field machinery and of the
handbook assumption 'Tm intrinsic superposition parameters = Er'.

  1. Gruber's fitted B_kq through our Tm Hamiltonian -> his calculated Stark pattern (machinery check).
  2. YAG dodecahedral site (Y 24c, Ia-3d, a = 12.000 A, O 96h (-0.0306, 0.0512, 0.1500)) + handbook O
     superposition parameters (taken from Er:YAG; NOT fitted to Tm) -> Tm:YAG levels vs measured.
  3. Tm-O intrinsic B2, B4, B6 fitted to the measured levels -> ratio to the Er values.
  4. Ground-state absorption line intensities into 3F4 (15 K, Gruber Table III) vs isotropic per-component
     Judd-Ofelt and site-derived (odd-rank superposition) intensities.

    python scripts/validate_tm_yag.py -> ml_results/tm_yag/tm_yag.json
"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from ase.spacegroup import crystal
from ase.geometry import get_distances
from scipy.optimize import minimize
from gainmat.forward.ions import get_ion
from gainmat.forward.spectra import spm_bkq, intensity_B, site_omegas, aniso_S
from gainmat.refs.spm_handbook import SPM_FOR_MODEL as H
from gainmat.refs.tm_yag_refs import GRUBER1989 as G

yag = crystal(["Y", "Al", "Al", "O"], basis=[(0.125, 0, 0.25), (0, 0, 0), (0.375, 0, 0.25), (-0.0306, 0.0512, 0.1500)],
              spacegroup=230, cellpar=[12.0, 12.0, 12.0, 90, 90, 90])
pos = yag.get_positions(); sym = yag.get_chemical_symbols(); iy = sym.index("Y")
D, d = get_distances([pos[iy]], pos, cell=yag.cell, pbc=True)
site = [("O", *map(float, D[0, j])) for j in np.argsort(d[0]) if sym[j] == "O" and 0.1 < d[0, j] < 2.6]
m = get_ion("Tm3+"); MULT = ("3H6", "3F4", "3H5", "3H4"); out = {}


def levels(bkq):
    E, W, par = m.stark(bkq)
    return {lab: np.sort(E[par == m.mult_index(lab)]) for lab in MULT}, (E, W, par)


def resid(lv):
    r = []
    for lab in MULT:
        e = lv[lab] - lv[lab][0]; ex = G["levels"][lab]
        r += [e[i] - (x - ex[0]) for i, x in enumerate(ex) if x is not None and i > 0]
    return np.array(r)


lvG, _ = levels({k: complex(v) for k, v in G["bkq"].items()})
out["machinery"] = {lab: dict(ours=(lvG[lab] - lvG[lab][0]).round(0).tolist(),
                              gruber_calc=(np.array(G["calc_levels"][lab]) - G["calc_levels"][lab][0]).tolist()) for lab in MULT}
out["machinery_rms_vs_expt"] = float(np.sqrt(np.mean(resid(lvG) ** 2)))
hb = dict(H["O"]["Tm"]); lvH, _ = levels(spm_bkq(site, {"O": hb}))
spl = lambda lv: {lab: float(lv[lab][-1] - lv[lab][0]) for lab in MULT}
out["handbook_on_yag"] = dict(rms=float(np.sqrt(np.mean(resid(lvH) ** 2))), splitting=spl(lvH),
                              measured_splitting={lab: max(x for x in G["levels"][lab] if x is not None) - G["levels"][lab][0] for lab in MULT})
f = lambda x: float(np.sqrt(np.mean(resid(levels(spm_bkq(site, {"O": dict(hb, B2=x[0], B4=x[1], B6=x[2])}))[0]) ** 2)))
b = minimize(f, [hb["B2"], hb["B4"], hb["B6"]], method="Nelder-Mead", options=dict(xatol=1, fatol=0.05, maxfev=400))
out["fitted_tm_intrinsics"] = dict(B2=b.x[0], B4=b.x[1], B6=b.x[2], rms=b.fun,
                                   ratio_to_er=(b.x / np.array([hb["B2"], hb["B4"], hb["B6"]])).tolist())
# line intensities from the ground level into 3F4
_, (E, W, par) = levels(spm_bkq(site, {"O": hb}))
i0 = int(np.argmin(E)); iF = np.where(par == m.mult_index("3F4"))[0]; iF = iF[np.argsort(E[iF])]
Uw = {lam: [W.conj().T @ Uq @ W for Uq in m.Uq[lam]] for lam in (2, 4, 6)}
meas = np.array([0.05, 0.57, np.nan, 0.94, 0.06, np.nan, 0.93, 0.13, 0.01])
OM = np.array([0.7, 1.3, 1.0])
iso = np.array([sum(OM[i] * sum(abs(X[i0, j]) ** 2 for X in Uw[lam]) for i, lam in enumerate((2, 4, 6))) for j in iF])
B = intensity_B(site, None, None, 2.3); om = site_omegas(B)
S = aniso_S(Uw, B, {lam: OM[i] / om[lam] for i, lam in enumerate((2, 4, 6))}); an = np.array([S[i0, j] for j in iF])
ok = ~np.isnan(meas)
corr = lambda p: float(np.corrcoef(np.log10(np.maximum(p[ok] / p.max(), 1e-3)), np.log10(meas[ok]))[0, 1])
out["intensities_3F4"] = dict(measured=meas.tolist(), isotropic=(iso / iso.max()).round(3).tolist(), site_derived=(an / an.max()).round(3).tolist(),
                              log_corr_isotropic=corr(iso), log_corr_site=corr(an))
os.makedirs("ml_results/tm_yag", exist_ok=True)
json.dump(out, open("ml_results/tm_yag/tm_yag.json", "w"), indent=1, default=float)
print(json.dumps({k: out[k] for k in ("machinery_rms_vs_expt", "handbook_on_yag", "fitted_tm_intrinsics")}, indent=1, default=float))
print("intensity log-correlations: isotropic", round(out["intensities_3F4"]["log_corr_isotropic"], 2), " site-derived", round(out["intensities_3F4"]["log_corr_site"], 2))
