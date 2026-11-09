"""Forward model calibrated by a measured absorption spectrum ("one spectrum" mode).

    python scripts/predict_with_absorption.py --ion Er3+ --family aluminosilicate \
        --band 4I13/2,4I15/2 --abs my_abs.csv [--em my_em.csv] [--ml] [--N 1e20] [--out out_dir]

CSV files: two columns, wavelength (nm) and cross section (cm^2); '#' comment lines allowed.
What changes with the measured absorption: that band's absorption is the measured curve, its
emission is the McCumber transform (eps from the model's Stark structure) unless --em is given,
and its radiative rate is rescaled to the measured band strength (which feeds the lifetimes).
Everything else (other bands, ESA, non-radiative rates, gain windows) comes from the model.
Writes <out>/prediction.json and <out>/spectra.png comparing model-only and calibrated runs.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np


def read_csv(p):
    rows = [l.replace(",", " ").split() for l in open(p) if l.strip() and not l.lstrip().startswith("#")]
    a = np.array([[float(x[0]), float(x[1])] for x in rows if len(x) >= 2])
    return a[:, 0].tolist(), a[:, 1].tolist()


def summary(res, band):
    spec = res["_spec"]; nu = spec["nu"]; lam = 1e7 / nu
    sa, se = spec["sigma_abs"][band], spec["sigma_em"][band]
    out = dict(abs_peak_nm=float(lam[np.argmax(sa)]), sigma_abs_peak=float(sa.max()),
               em_peak_nm=float(lam[np.argmax(se)]), sigma_em_peak=float(se.max()),
               windows_nm=res.get("windows"), measured_bands=spec.get("measured_bands"))
    taus = {}
    for k, v in res["pumps"].items():
        if "tau_ms" in v:
            taus = v["tau_ms"]; break
    out["lifetimes_ms"] = {k: round(float(v), 4) for k, v in taus.items()}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ion", required=True)
    ap.add_argument("--family", required=True)
    ap.add_argument("--band", required=True, help="upper,lower multiplet labels, e.g. 4I13/2,4I15/2")
    ap.add_argument("--abs", required=True)
    ap.add_argument("--em", default="")
    ap.add_argument("--N", type=float, default=None, help="dopant density, cm^-3")
    ap.add_argument("--ml", action="store_true", help="also switch on the ML layer")
    ap.add_argument("--sites", type=int, default=30)
    ap.add_argument("--out", default="forward_results/with_absorption")
    a = ap.parse_args()
    if a.ml:
        os.environ["GAINMAT_ML"] = "1"
    from gainmat.forward import glass as g
    g.USE_ML = a.ml
    g.load_refs(); g.apply_site_refs(); g.calibrate_spm(370.0, g.CARNALL_ER_LAF3_BKQ, verbose=False)
    band = tuple(x.strip() for x in a.band.split(","))
    meas = read_csv(a.abs)
    val = (meas[0], meas[1]) + ((read_csv(a.em),) if a.em else ())
    base = dict(N_cm3=a.N) if a.N else {}
    r0 = g.run(a.ion, a.family, n_sites=a.sites, overrides=dict(base))
    r1 = g.run(a.ion, a.family, n_sites=a.sites, overrides=dict(base, measured_absorption={band: val}))
    os.makedirs(a.out, exist_ok=True)
    out = dict(ion=a.ion, family=a.family, band=list(band), model_only=summary(r0, band),
               calibrated=summary(r1, band), abs_file=a.abs, em_file=a.em or None, ml=a.ml)
    json.dump(out, open(os.path.join(a.out, "prediction.json"), "w"), indent=1, default=float)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    for r, ls, tag in ((r0, "--", "model only"), (r1, "-", "with measured absorption")):
        s = r["_spec"]; lam = 1e7 / s["nu"]
        ax.plot(lam, s["sigma_abs"][band], ls, color="#2b6cb0", label=f"absorption, {tag}")
        ax.plot(lam, s["sigma_em"][band], ls, color="#c05621", label=f"emission, {tag}")
    ax.plot(meas[0], meas[1], "k.", ms=3, label="measured absorption")
    lo, hi = min(meas[0]), max(meas[0]); pad = 0.15 * (hi - lo)
    ax.set_xlim(lo - pad, hi + pad); ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("cross section (cm$^2$)")
    ax.set_title(f"{a.ion} {a.family}: {band[0]} - {band[1]}", fontsize=9); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "spectra.png"), dpi=140)
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
