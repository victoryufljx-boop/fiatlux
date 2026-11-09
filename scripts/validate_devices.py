"""Device-level validation: model fibre amplifiers and lasers built to the published
geometry, doping and pumping, compared with the measured gain / slope / threshold.

    python scripts/validate_devices.py [--sites 30]  ->  forward_results/device_validation.json

Every number describing a device comes from the paper (see gainmat/refs/*_refs.py); nothing
is fitted to the device result.  Variants:
  single-ion : the forward model as is (intrinsic lifetimes)
  observed   : plus the extra first-order decay needed to reach the lifetime measured in a
               comparable fibre (concentration / OH quenching that the single-ion model lacks)
  +CR        : plus Tm 3H4,3H6 -> 3F4,3F4 cross-relaxation (and its inverse) with the
               Jackson & King (1999) coefficients scaled linearly with Tm concentration
"""
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np

from gainmat.forward.glass import load_refs, apply_site_refs, calibrate_spm, CARNALL_ER_LAF3_BKQ
from gainmat.forward.fiber import local_model_for, amplifier, laser, slope_and_threshold

ap = argparse.ArgumentParser()
ap.add_argument("--sites", type=int, default=30)
ap.add_argument("--out", default="forward_results")
ap.add_argument("--only", default="")
a = ap.parse_args()
load_refs(); apply_site_refs(); calibrate_spm(370.0, CARNALL_ER_LAF3_BKQ, verbose=False)
only = set(a.only.split(",")) - {""}
rows = []


def done(name, **kw):
    kw["case"] = name
    rows.append(kw)
    print(json.dumps(kw, default=float))
    json.dump(rows, open(f"{a.out}/device_validation.json", "w"), indent=1, default=float)


def tau_quench(lm, level, tau_obs_s):
    i = lm.labels.index(level)
    k0 = lm.A[i].sum() + lm.W[i]
    return max(1.0 / tau_obs_s - k0, 0.0)


# ------------------------------------------------------------------ 1. Ohishi 1991 Pr:ZBLAN amp
if not only or "pr_amp" in only:
    t0 = time.time()
    N = 500e-6 * 4.33 / 140.9 * 6.022e23            # 500 wppm Pr in ZBLAN (rho 4.33 g/cm3)
    NA = 1.5 * np.sqrt(2 * 0.006)                     # Delta n 0.6 %
    for var in ("single-ion", "observed"):
        lm = local_model_for("Pr3+", "ZBLAN", n_sites=a.sites)
        if var == "observed":
            lm = local_model_for("Pr3+", "ZBLAN", n_sites=a.sites,
                                 quench={"1G4": tau_quench(lm, "1G4", 110e-6)})
        g = [amplifier(lm, 8, 4.0, NA, N, [(1017, P)], 1310, 5e-7, nz=80, backward_pump=True)["gain_dB"]
             for P in (0.06, 0.12, 0.18)]
        done("Pr:ZBLAN amplifier (Ohishi et al. 1991)", variant=var, quantity="net gain at 1310 nm, 180 mW",
             model=g[-1], measured=5.2, model_slope_dB_per_mW=(g[2] - g[0]) / 120, measured_slope_dB_per_mW=0.03,
             note="8 m, 4 um core, Delta n 0.6 %, 500 wppm, 1017 nm backward pump; measured internal gain 7.1 dB",
             seconds=time.time() - t0)

# ------------------------------------------------------- 2. Henderson-Sapir 2014 Er:ZBLAN 3.5 um
if not only or "er_dwp" in only:
    t0 = time.time()
    N = 0.01 * 4.33 / 142.4 * 6.022e23               # 1 mol% ErF3 in ZBLAN (~1.8e20 cm^-3)
    lm = local_model_for("Er3+", "ZBLAN", n_sites=a.sites)
    out = []
    for P2 in (0.3, 0.8, 1.4, 2.0, 2.4):
        r = laser(lm, 0.18, 10.0, 0.15, N, [(985, 0.175), (1973, P2)], 3500, R1=0.99, R2=0.90,
                  nz=40, pump_overlap=0.8)
        out.append((P2, r["P_out_W"]))
    x = np.array([p for p, o in out if o > 1e-4]); y = np.array([o for p, o in out if o > 1e-4])
    s = float(np.polyfit(x, y, 1)[0]) if len(x) > 1 else None
    done("Er:ZBLAN 3.5 um dual-pump laser (Henderson-Sapir et al. 2014)", variant="single-ion",
         quantity="slope vs launched 1973 nm power; output at 2.4 W", model_slope=s, measured_slope=0.254,
         model_P_out_W=out[-1][1], measured_P_out_W=0.26, curve=out,
         note="18 cm, 10 um core, NA 0.15, 1 mol% Er, 985 nm 175 mW launched, R 0.99 / 0.90; "
              "measured slope is vs INCIDENT P2 (52 % vs absorbed); pump overlap 0.8 assumed (multimode)",
         seconds=time.time() - t0)

# --------------------------------------------------------------- 3. Jackson & Mossman 2003 Tm lasers
if not only or "tm_laser" in only:
    from gainmat.refs.tm_silica_refs import JACKSON2003 as J
    fib = {f["fibre"]: f for f in J["fibres"]}
    fig5 = J["fig5_output_vs_launched_pump_digitized"]
    cases = [(4, "fibre4_1.6wt", 0.46, 5.24), (3, "fibre3_1.3wt", 0.415, 5.99), (1, "fibre1_0.7wt", 0.213, 7.41)]
    for k, key, s_meas, th_meas in cases:
        f = fib[k]; L = fig5[key]["L_m_est"]
        N = f["N_Tm_cm3_from_alpha_core"]; core = min(max(f["core_diam_um_est_from_area_ratio"], 20.0), 25.0)
        Gp = (core / 300.0) ** 2 * 4 / np.pi / (1 - 0.1)  # D-shaped 300 um clad (~10 % cut), core/clad area
        k3101 = 1.8e-16 * N / 4.5e20                       # Jackson & King 1999, linear in N
        for var in ("single-ion", "observed", "observed+CR", "ML-L3", "ML-L3+CR"):
            t0 = time.time()
            et = [("3H4", "3F4", "3H6", "3F4", k3101), ("3F4", "3H4", "3F4", "3H6", 0.084 * k3101)] \
                if "CR" in var else None
            lm0 = local_model_for("Tm3+", "aluminosilicate", n_sites=a.sites)
            if var.startswith("ML"):
                from gainmat.ml.quench import extra_rates
                q = extra_rates("aluminosilicate", lm0.spec["bary"], N)
            else:
                q = {"3F4": tau_quench(lm0, "3F4", 0.65e-3)} if var != "single-ion" else None
            lm = local_model_for("Tm3+", "aluminosilicate", n_sites=a.sites, quench=q, energy_transfer=et)
            s, th, curve = slope_and_threshold(lm, L, core, 0.19, N, 805, [6, 9, 12, 15], 1990,
                                               R1=0.99, R2=0.035, alpha_bg_dB_m=0.1, nz=40,
                                               pump_overlap=Gp)
            done(f"Tm:silica cladding-pumped laser, Jackson & Mossman 2003 fibre {k}", variant=var,
                 quantity="slope (launched) / threshold W", model_slope=s, measured_slope=s_meas,
                 model_threshold_W=th, measured_threshold_W=th_meas, curve=curve,
                 note=f"L {L} m (derived), core {core:.1f} um, N {N:.3g}, Gamma_p {Gp:.4f}, 805 nm, R 0.99/0.035, "
                      f"bg 0.1 dB/m; k3101 {k3101:.2g} cm3/s", seconds=time.time() - t0)

# ------------------------------------------------------------------ 4. Le Gouet 2020 Ho:silica laser
if not only or "ho_laser" in only:
    from gainmat.refs.midir_ho_refs import LEGOUET2020 as LG
    N = 5.3e19; core = 20.0
    Gp = np.pi * (core / 2) ** 2 / (0.8284 * 250.0 ** 2)   # octagon 250 um flat-to-flat
    for var in ("single-ion", "observed", "ML-L3"):
        t0 = time.time()
        lm0 = local_model_for("Ho3+", "aluminosilicate", n_sites=a.sites)
        if var == "ML-L3":
            from gainmat.ml.quench import extra_rates
            q = extra_rates("aluminosilicate", lm0.spec["bary"], N)
        else:
            q = {"5I7": tau_quench(lm0, "5I7", 0.7e-3)} if var == "observed" else None
        lm = local_model_for("Ho3+", "aluminosilicate", n_sites=a.sites, quench=q)
        s, th, curve = slope_and_threshold(lm, 16, core, 0.08, N, 1940, [50, 100, 150, 200], 2120,
                                           R1=0.99, R2=0.03, alpha_bg_dB_m=0.04, nz=40, pump_overlap=Gp)
        done("Ho:silica cladding-pumped laser (Le Gouet et al. 2020)", variant=var,
             quantity="slope (coupled) / threshold W", model_slope=s, measured_slope=0.49,
             model_threshold_W=th, measured_threshold_W=13, curve=curve,
             note=f"16 m, 20 um core NA 0.08, 250 um octagonal clad (Gamma_p {Gp:.4f}), 5.3e19 cm^-3, 1940 -> 2120 nm, "
                  "FBG 0.99 / Fresnel 0.03, bg 0.04 dB/m; 15-23 % of ions in pairs not modelled",
             seconds=time.time() - t0)

# ----------------------------------------------------------------------- 5. Pask 1995 Yb lasers
if not only or "yb_laser" in only:
    N = 550e-6 * 2.2e22                                 # 550 ppm (molar, per SiO2 unit)
    lm = local_model_for("Yb3+", "germanosilicate", n_sites=a.sites)
    t0 = time.time()
    s, th, curve = slope_and_threshold(lm, 90, 3.75, 0.17, N, 1047, [0.05, 0.1, 0.2, 0.3], 1102,
                                       R1=0.98, R2=0.04, alpha_bg_dB_m=0.005, nz=60)
    done("Yb:germanosilicate 1102 nm laser (Pask et al. 1995)", variant="single-ion",
         quantity="slope (launched) / threshold W", model_slope=s, measured_slope=0.90,
         model_threshold_W=th, measured_threshold_W=0.030, curve=curve,
         note="90 m, 3.75 um core, NA 0.17, 550 ppm (taken as 1.2e19 cm^-3), 1047 nm pump, R 0.98 / Fresnel",
         seconds=time.time() - t0)
    t0 = time.time()
    s, th, curve = slope_and_threshold(lm, 5, 3.75, 0.17, N, 840, [0.08, 0.15, 0.25, 0.35], 1040,
                                       R1=0.04, R2=0.04, alpha_bg_dB_m=0.005, nz=40)
    done("Yb:germanosilicate free-running laser, 5 m (Pask et al. 1995 Fig. 4)", variant="single-ion",
         quantity="threshold W (launched)", model_slope=s, model_threshold_W=th,
         measured_threshold_W=0.045, curve=curve,
         note="5 m, Fresnel both ends, 840 nm pump; measured minimum threshold ~45 mW (digitised)",
         seconds=time.time() - t0)
print("wrote", f"{a.out}/device_validation.json")
