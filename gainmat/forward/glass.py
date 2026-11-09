"""Glass host descriptions for the forward model, and the one-call driver.

A host is described by what the physics needs, not by a name:
  n          refractive index (local-field factors)
  hw_max     max phonon energy (multiphonon decay)          gap_family: energy-gap-law set
  ligand / CN / R / R_std / ang_std : the dopant's first-shell geometry ensemble
  loss_dB_m  background loss floor of real fibre in that family
  omegas     Judd-Ofelt (Omega2, Omega4, Omega6) per ion (measured; the one non-ab-initio input)

Defaults below are overwritten by gainmat.refs.glass_refs / site_refs when those modules
provide sourced values (see load_refs()).  Every value records where it came from.
"""
from __future__ import annotations

import copy
import os

import numpy as np

from .spectra import random_sites, ensemble_spectra, nu_grid, peak, fwhm_nm
from .gain import gain_spectra, envelope, windows

# Superposition-model intrinsic parameters (Wybourne, cm^-1, R0 in Angstrom).
# Placeholder until site_refs provides sourced values; scaled per ion in spm_for().
SPM_DEFAULT = {
    "O": dict(R0=2.30, B2=400.0, B4=350.0, B6=150.0, t2=5.0, t4=6.0, t6=6.0, src="placeholder"),
    "F": dict(R0=2.30, B2=300.0, B4=300.0, B6=130.0, t2=5.0, t4=6.0, t6=6.0, src="placeholder"),
    "S": dict(R0=2.80, B2=350.0, B4=250.0, B6=100.0, t2=5.0, t4=6.0, t6=6.0, src="placeholder"),
}

FAMILIES = {
    # name: n, hw_max, gap_family, ligand, CN, R, R_std, ang_std, loss floor (dB/m)
    "silica":          dict(n=1.444, hw_max=1100, gap_family="fused_silica", ligand="O", CN=6.0, R=2.25, R_std=0.07, ang_std=12, loss_dB_m=0.005),
    "aluminosilicate": dict(n=1.46, hw_max=1100, gap_family="aluminosilicate", ligand="O", CN=6.5, R=2.27, R_std=0.09, ang_std=15, loss_dB_m=0.01),
    "germanosilicate": dict(n=1.46, hw_max=1000, gap_family="germanosilicate", ligand="O", CN=6.0, R=2.26, R_std=0.08, ang_std=13, loss_dB_m=0.01),
    "phosphosilicate": dict(n=1.46, hw_max=1300, gap_family="phosphate", ligand="O", CN=6.5, R=2.27, R_std=0.08, ang_std=13, loss_dB_m=0.01),
    "phosphate":       dict(n=1.53, hw_max=1200, gap_family="phosphate", ligand="O", CN=6.5, R=2.28, R_std=0.08, ang_std=12, loss_dB_m=0.05),
    "germanate":       dict(n=1.70, hw_max=850, gap_family="germanate", ligand="O", CN=7.0, R=2.30, R_std=0.10, ang_std=15, loss_dB_m=0.5),
    "tellurite":       dict(n=2.00, hw_max=750, gap_family="tellurite", ligand="O", CN=7.0, R=2.32, R_std=0.12, ang_std=18, loss_dB_m=0.2),
    "ZBLAN":           dict(n=1.50, hw_max=600, gap_family="fluoride", ligand="F", CN=8.0, R=2.27, R_std=0.08, ang_std=12, loss_dB_m=0.01),
    "fluoroindate":    dict(n=1.50, hw_max=510, gap_family=None, ligand="F", CN=8.0, R=2.27, R_std=0.08, ang_std=12, loss_dB_m=0.05),
    "chalcogenide":    dict(n=2.40, hw_max=350, gap_family=None, ligand="S", CN=6.5, R=2.80, R_std=0.10, ang_std=15, loss_dB_m=0.5),
}
for _f in FAMILIES.values():
    _f.setdefault("src", "placeholder (see load_refs)")

# Judd-Ofelt parameters, 1e-20 cm^2.  Seeded with ln_refs values; glass_refs overrides.
OMEGA = {
    ("Er3+", "silica"): (4.23, 1.04, 0.61),
    ("Pr3+", "ZBLAN"): (1.431, 4.221, 4.872),     # Ohishi et al., Opt. Lett. 16, 1747 (1991)
    # Wetenkamp, West & Toebben, J. Non-Cryst. Solids 140, 35 (1992), Table 1
    ("Ho3+", "ZBLAN"): (2.30, 2.30, 1.71),
    ("Tm3+", "ZBLAN"): (2.70, 1.59, 1.09),
    ("Nd3+", "ZBLAN"): (2.10, 3.71, 4.62),
    ("Tm3+", "silica"): (6.23, 1.91, 1.36),       # Walsh & Barnes, Appl. Phys. B 78, 325 (2004)
    ("Er3+", "ZBLAN"): (2.54, 1.39, 0.97),
}

# typical dopant density used for dB/m (ions/cm^3) and pump wavelengths per ion
N_DEFAULT = {"Er3+": 1e19, "Yb3+": 5e19, "Tm3+": 5e19, "Ho3+": 5e19, "Pr3+": 2e19,
             "Nd3+": 3e19, "Dy3+": 5e19}
PUMPS = {"Er3+": [980.0, 1480.0], "Yb3+": [915.0, 976.0], "Tm3+": [790.0, 1210.0, 1600.0],
         "Ho3+": [1150.0, 1950.0], "Pr3+": [1017.0, 444.0], "Nd3+": [808.0],
         "Dy3+": [1100.0, 1700.0]}
EMAX = {"Er3+": 21000.0, "Yb3+": 12000.0, "Tm3+": 22000.0, "Ho3+": 19000.0,
        "Pr3+": 23000.0, "Nd3+": 19000.0, "Dy3+": 14000.0}


def yb_omega_from_lifetime(tau_s: float, n: float):
    """Yb3+ has a single 2F5/2 -> 2F7/2 transition, so Judd-Ofelt Omegas are not separately
    defined.  Use equal Omega_2 = Omega_4 = Omega_6 = x chosen so that the radiative rate
    (ED + MD) equals 1/tau of the measured 2F5/2 lifetime (Yb has no nonradiative channel in
    these glasses, so tau ~ tau_rad).  This is an input normalisation, documented as such."""
    from .ions import get_ion
    from .gain import multiplet_rates
    m = get_ion("Yb3+")
    E = {mm["label"]: mm["energy"] for mm in m.multiplets}
    a_md = multiplet_rates("Yb3+", (0, 0, 0), n, E)[1, 0]
    a_1 = multiplet_rates("Yb3+", (1, 1, 1), n, E)[1, 0] - a_md
    x = max((1.0 / tau_s - a_md) / a_1, 0.0)
    return (x, x, x)


def load_refs():
    """Merge sourced values from refs modules if present (silently keeps defaults otherwise)."""
    loaded = []
    try:
        from ..refs import glass_refs as G
        for fam, d in getattr(G, "GLASSES", {}).items():
            if fam not in FAMILIES:
                continue
            for key, tgt in (("n_D", "n"), ("phonon_max_cm1", "hw_max"),
                             ("fiber_loss_dB_per_m_min", "loss_dB_m")):
                v = d.get(key)
                v = v.get("value") if isinstance(v, dict) else v
                if isinstance(v, (int, float)) and v > 0:
                    FAMILIES[fam][tgt] = float(v)
            FAMILIES[fam]["src"] = "glass_refs"
        for key, d in getattr(G, "OMEGA", {}).items():
            ion_el, fam = key.split(":")
            ion = ion_el if ion_el.endswith("3+") else ion_el + "3+"
            try:
                o = tuple(float(d[k]["value"] if isinstance(d[k], dict) else d[k]) for k in ("O2", "O4", "O6"))
                OMEGA[(ion, fam)] = o
            except Exception:
                tau = d.get("tau_2F52_s")
                tau = tau.get("value") if isinstance(tau, dict) else tau
                if ion == "Yb3+" and tau and fam in FAMILIES:
                    OMEGA[(ion, fam)] = yb_omega_from_lifetime(float(tau), FAMILIES[fam]["n"])
        OMEGA.update(OMEGA_PRIMARY)          # values read from the primary papers win
        loaded.append("glass_refs")
    except ImportError:
        pass
    try:
        from ..refs import site_refs as S
        for lig, d in getattr(S, "SPM_FOR_MODEL", {}).items():
            SPM_DEFAULT[lig] = dict(d)
        for fam, d in getattr(S, "SITES_FOR_MODEL", {}).items():
            if fam in FAMILIES:
                FAMILIES[fam].update({k: d[k] for k in ("CN", "R", "R_std", "ang_std") if k in d})
        loaded.append("site_refs")
    except ImportError:
        pass
    return loaded


# Multiphonon IR absorption edge of the glass:  alpha(lambda) = exp(a (1/l1 - 1/lambda)) dB/m,
# with l1 the wavelength where it reaches 1 dB/m.  Anchored on silica fibre,
# alpha = 7.81e11 exp(-48.48/lambda[um]) dB/km (the standard multiphonon fit used for silica
# fibre loss), giving a = 48.48 um and l1 = 2.37 um.  For other glasses the edge is scaled by
# the max phonon energy (edge wavenumber ~ multiple of hw_max): l1 = 2.37 um * 1100/hw_max,
# a = 48.48 um * 1100/hw_max.  This places ZBLAN at ~4.3 um, fluoroindate ~5.1 um,
# tellurite ~3.5 um, germanate ~3.1 um and sulfide glass ~7.4 um.
def ir_loss_dB_m(lam_nm, hw_max):
    lam_um = np.asarray(lam_nm) / 1000.0
    s = 1100.0 / hw_max
    a, l1 = 48.48 * s, 2.37 * s
    return np.exp(np.clip(a * (1.0 / l1 - 1.0 / lam_um), -50, 50))


# Extra inhomogeneous broadening: per-site fractional scatter of multiplet energies
# (nephelauxetic variation).  0.27 % makes the Er 4I13/2-4I15/2 zero-phonon line 49 cm^-1 wide
# (FWHM) in aluminosilicate, as measured at 4.2 K by Desurvire & Simpson, OL 15, 547 (1990);
# the crystal-field disorder alone gives 27 cm^-1.  Applied to every family (fluorides assumed).
NEPH_STD = 0.0027
# One-phonon vibronic sideband strength (fraction of each line).  0 = off unless calibrated.
SIDEBAND_S = 0.0
# Learned host microstructure (ML layer L2, gainmat/ml/host_calib.py): when USE_ML is set,
# run() applies ml_results/l2/learned_theta.json (cf_scale, disorder, neph_std, gamma_h)
# to the families that have been calibrated.  Off by default.
USE_ML = False
_LEARNED_THETA = None


def learned_theta(family):
    global _LEARNED_THETA
    if _LEARNED_THETA is None:
        import json as _json
        p = os.path.join(os.path.dirname(__file__), "..", "..", "ml_results", "l2", "learned_theta.json")
        _LEARNED_THETA = _json.load(open(p)) if os.path.exists(p) else {}
    return _LEARNED_THETA.get(family)


# dual-wavelength pumping: (first pump from the ground state, second pump between excited states)
# Er 985 + 1973 nm: 4I15/2 -> 4I11/2, then 4I11/2 -> 4F9/2 for 3.5 um gain on 4F9/2 -> 4I9/2
# (Henderson-Sapir, Munch & Ottaway, Opt. Lett. 39, 493 (2014)).
DUAL_PUMPS = {"Er3+": [(985.0, 1973.0)]}
# pump intensities: up to ~1 W in a ~50 um^2 single-mode core (2e6 W/cm^2)
I_PUMP = np.logspace(2, np.log10(2e6), 10)


# measured sets the data agent flagged as unreliable (garbled exponent / implausibly low)
OMEGA_SUSPECT = {("Er3+", "chalcogenide")}
# Omega read directly from primary papers (kept over compilation values by load_refs)
OMEGA_PRIMARY = {k: OMEGA[k] for k in [("Pr3+", "ZBLAN"), ("Ho3+", "ZBLAN"), ("Tm3+", "ZBLAN"),
                                       ("Nd3+", "ZBLAN"), ("Tm3+", "silica")]}


def omegas_for(ion, family):
    """Measured Omega for (ion, family) when available and not flagged; otherwise the
    LOMS-trained composition model's prediction for the family's representative composition;
    otherwise fall back to the same ion in the nearest family (flagged)."""
    if (ion, family) in OMEGA and (ion, family) not in OMEGA_SUSPECT:
        return OMEGA[(ion, family)], "measured"
    try:
        from ..refs.omega_family_pred import OMEGA_PRED
        v = OMEGA_PRED.get(f"{ion}:{family}")
        if v:
            return tuple(v), "predicted (LOMS-trained model)"
    except ImportError:
        pass
    near = {"aluminosilicate": "silica", "germanosilicate": "silica", "phosphosilicate": "silica",
            "fluoroindate": "ZBLAN", "germanate": "tellurite", "tellurite": "germanate",
            "phosphate": "silica", "chalcogenide": "tellurite", "silica": "aluminosilicate",
            "ZBLAN": "fluoroindate"}
    f2 = near.get(family)
    if f2 and (ion, f2) in OMEGA:
        return OMEGA[(ion, f2)], f"borrowed from {ion}:{f2}"
    same = [v for (i, f), v in OMEGA.items() if f == family]
    if same:
        return tuple(np.mean(same, axis=0)), f"family mean over other ions ({family})"
    allv = [v for (i, f), v in OMEGA.items() if i == ion]
    if allv:
        return tuple(np.mean(allv, axis=0)), f"mean of {ion} over families"
    return (3.0, 1.5, 1.0), "generic default"


def run(ion, family, n_sites=40, pumps=None, overrides=None, seed=0, threshold_dB_m=0.1,
        lam=(400.0, 5000.0), dnu=2.0, intensities=I_PUMP):
    """Forward model for one (ion, glass): spectra, lifetimes, efficiencies, gain windows."""
    fam = copy.deepcopy(FAMILIES[family])
    if USE_ML and learned_theta(family):
        th = learned_theta(family)
        fam.update(cf_scale=th["cf_scale"], R_std=fam["R_std"] * th["disorder"],
                   ang_std=fam["ang_std"] * th["disorder"], neph_std=th["neph_std"], gamma_h=th["gamma_h"])
    fam.update(overrides or {})
    om, om_src = omegas_for(ion, family) if "omegas" not in fam else (fam["omegas"], "override")
    lig = fam["ligand"]
    sites = fam.get("sites") or random_sites(n_sites, lig, fam["CN"], fam["R"], fam["R_std"],
                                             fam["ang_std"], seed=seed)
    spm = {lig: spm_for(ion, lig)}
    if fam.get("cf_scale", 1.0) != 1.0:                 # host crystal-field strength multiplier
        spm[lig] = dict(spm[lig], **{f"B{k}": spm[lig][f"B{k}"] * fam["cf_scale"] for k in (2, 4, 6)})
    emax = fam.get("emax", EMAX[ion])
    spec = ensemble_spectra(ion, sites, spm, om, fam["n"], nu=nu_grid(lam[0], lam[1], dnu),
                            emax=emax, T=fam.get("T", 295.0), gamma_h=fam.get("gamma_h", 20.0),
                            neph_std=fam.get("neph_std", NEPH_STD), sideband_S=fam.get("sideband_S", SIDEBAND_S),
                            hw_max=fam["hw_max"], seed=seed)
    N = fam.get("N_cm3", N_DEFAULT[ion])
    res = dict(ion=ion, family=family, omegas=om, omega_src=om_src, n=fam["n"], hw_max=fam["hw_max"],
               loss_dB_m=fam["loss_dB_m"], N_cm3=N, pumps={}, split=spec["split"])
    env_all = None
    loss = fam["loss_dB_m"] + ir_loss_dB_m(1e7 / spec["nu"], fam["hw_max"])
    xW = None
    if USE_ML:                                          # L3 extrinsic quenching
        from ..ml.quench import extra_rates
        xW = extra_rates(family, spec["bary"], N) or None
    for p in (PUMPS[ion] + DUAL_PUMPS.get(ion, []) if pumps is None else pumps):
        p1, p2 = (p, None) if not isinstance(p, (tuple, list)) else (float(p[0]), float(p[1]))
        key = p if p2 is None else f"{p1:g}+{p2:g}"
        try:
            gs = gain_spectra(ion, spec, om, fam["n"], fam["hw_max"], fam.get("gap_family"),
                              N_tot_cm3=N, pump_nm=p1, intensities_W_cm2=intensities, emax=emax,
                              pump2_nm=p2, T=fam.get("T", 295.0), extra_W=xW)
        except ValueError as e:
            res["pumps"][key] = dict(error=str(e)); continue
        env = envelope(gs)
        env_all = env if env_all is None else np.maximum(env_all, env)
        res["pumps"][key] = dict(pump_level=gs["pump_level"],
                               windows=windows(gs["nu"], env, loss, threshold_dB_m),
                               tau_ms={k: v * 1e3 for k, v in gs["tau"].items()},
                               eta=gs["eta"])
    res["windows"] = windows(spec["nu"], env_all, loss, threshold_dB_m) if env_all is not None else []
    res["_spec"], res["_env"], res["_loss"] = spec, env_all, loss
    return res


def key_transitions(res, pairs):
    """Peak cross sections / FWHM for named (upper, lower) pairs."""
    spec = res["_spec"]; out = {}
    for p in pairs:
        if p not in spec["sigma_abs"]:
            continue
        la, sa = peak(spec["nu"], spec["sigma_abs"][p]); le, se = peak(spec["nu"], spec["sigma_em"][p])
        out[f"{p[0]}->{p[1]}"] = dict(abs_peak_nm=round(la, 1), sigma_abs=sa, em_peak_nm=round(le, 1),
                                      sigma_em=se, fwhm_em_nm=round(fwhm_nm(spec["nu"], spec["sigma_em"][p]), 1))
    return out


def load_md_sites(path, max_sites=None):
    """Sites from a doped-glass MD JSON (md_sites/<glass>.json on branch forward/md-sites):
    er_sites[i]['ligands'] = [[element, x, y, z], ...] relative to the dopant (Angstrom).
    Returns the list-of-ligand-tuples format used by spectra.ensemble_spectra."""
    import json as _json
    d = _json.load(open(path))
    sites = [[(l[0], float(l[1]), float(l[2]), float(l[3])) for l in s["ligands"]]
             for s in d.get("er_sites", []) if s.get("ligands")]
    return sites[:max_sites] if max_sites else sites


# ------------------------------------------------------------------ crystal-field parameters
# Ion-specific superposition-model parameters from Newman & Ng (eds.), Crystal Field Handbook
# (CUP 2000), ch. 5 (gainmat.refs.spm_handbook), multiplied by one strength scale per ligand
# type set by calibrate_spm().  Garnet/fluorite crystal values over-predict glass splittings
# about 2x, hence the scale; the per-ion ratios and the distance laws t_k are kept.
SPM_SCALE = {"O": 1.0, "F": 1.0, "S": 1.0}
USE_HANDBOOK = True


# Ion-specific strength factor on top of the per-ligand scale.  Yb3+: the Er-calibrated scale
# gives half the measured Stark splittings in silica glass (Pask et al., IEEE JSTQE 1, 2 (1995),
# Fig. 1(a): 2F7/2 1490, 2F5/2 1370 cm^-1, germanosilicate).  Yb3+ has the lowest
# ligand->metal charge-transfer energy of the trivalent lanthanides after Eu3+, and the
# charge-transfer (covalent) contribution to the crystal field grows as that energy falls,
# so a larger effective field is physical.  Set from the 2F7/2 total splitting only;
# 2F5/2 is then a check.  calibrate_yb_spm() recomputes it.
ION_SPM_SCALE = {}


def calibrate_yb_spm(target_cm1=1490.0, family="germanosilicate", n_sites=30, verbose=True):
    """Scale Yb3+ SPM so the site-averaged 2F7/2 splitting matches target (Pask 1995)."""
    ION_SPM_SCALE["Yb"] = 1.0
    out = {}
    for _ in range(4):
        r = run("Yb3+", family, n_sites=n_sites, pumps=[])
        sp = r["split"]["2F7/2"]
        ION_SPM_SCALE["Yb"] *= target_cm1 / sp
        out = dict(split_before=sp, scale=ION_SPM_SCALE["Yb"])
        if abs(sp - target_cm1) < 10:
            break
    r = run("Yb3+", family, n_sites=n_sites, pumps=[])
    out.update(split_2F72=r["split"]["2F7/2"], split_2F52=r["split"]["2F5/2"],
               target_2F72=target_cm1, measured_2F52=1370.0)
    if verbose:
        print("Yb SPM scale %.3f: 2F7/2 %.0f (target %.0f), 2F5/2 %.0f (Pask 1370)" % (
            ION_SPM_SCALE["Yb"], out["split_2F72"], target_cm1, out["split_2F52"]))
    return out


def spm_for(ion: str, ligand: str) -> dict:
    el = ion.replace("3+", "")
    if USE_HANDBOOK:
        try:
            from ..refs.spm_handbook import SPM_FOR_MODEL as H
            if ligand in H and el in H[ligand]:
                p = dict(H[ligand][el])
                f = SPM_SCALE.get(ligand, 1.0) * ION_SPM_SCALE.get(el, 1.0)
                for k in (2, 4, 6):
                    p[f"B{k}"] = p[f"B{k}"] * f
                p["src"] = f"Newman & Ng handbook x{f:.3f}"
                return p
        except ImportError:
            pass
    return SPM_DEFAULT[ligand]


# ------------------------------------------------------------------ crystal-field calibration
def cf_strength(bkq):
    """Scalar crystal-field strength S = sqrt(sum_k 1/(2k+1) sum_q |B^k_q|^2) (cm^-1)."""
    s = 0.0
    for (k, q), B in bkq.items():
        s += (1 if q == 0 else 2) * abs(B) ** 2 / (2 * k + 1)
    return float(np.sqrt(s))


def _er_ground_split(fam, scale, n_sites=30, seed=11):
    from .ions import get_ion
    from .spectra import spm_bkq
    f = FAMILIES[fam]; lig = f["ligand"]
    p = dict(spm_for("Er3+", lig)); [p.__setitem__(f"B{k}", p[f"B{k}"] * scale) for k in (2, 4, 6)]
    m = get_ion("Er3+"); g = m.mult_index("4I15/2")
    sp = []
    for s in random_sites(n_sites, lig, f["CN"], f["R"], f["R_std"], f["ang_std"], seed=seed):
        E, W, par = m.stark(spm_bkq(s, {lig: p}))
        e = E[par == g]; sp.append(e.max() - e.min())
    return float(np.mean(sp))


def calibrate_spm(target_er_split_O=370.0, carnall_er_laf3=None, verbose=True):
    """One strength scale per ligand type; ion-specific ratios and t_k come from the handbook
    (or the generic SPM_DEFAULT shape when USE_HANDBOOK is False):
    O: Er3+ 4I15/2 total Stark splitting in silica sites = 370 cm^-1
       (glasses: 335-400 cm^-1, Zemon et al., J. Appl. Phys. 69, 6799 (1991));
    F: mean scalar crystal-field strength of Er3+ ZBLAN sites = that of Carnall's fitted
       Er:LaF3 crystal field (transfer from a fluoride-crystal fit);
    S: no glass data; the O scale is applied (unconstrained, flagged)."""
    from .spectra import spm_bkq
    out = {}
    for _ in range(2):                                  # splitting is ~linear in the scale
        a = _er_ground_split("silica", 1.0)
        r = target_er_split_O / a
        if USE_HANDBOOK:
            SPM_SCALE["O"] *= r
        else:
            for k in (2, 4, 6):
                SPM_DEFAULT["O"][f"B{k}"] *= r
    out["O"] = dict(scale=SPM_SCALE["O"], er_split_check=_er_ground_split("silica", 1.0))
    if carnall_er_laf3:
        S_ref = cf_strength(carnall_er_laf3)
        f = FAMILIES["ZBLAN"]
        S_now = np.mean([cf_strength(spm_bkq(s, {"F": spm_for("Er3+", "F")}))
                         for s in random_sites(30, "F", f["CN"], f["R"], f["R_std"], f["ang_std"], seed=5)])
        r = S_ref / S_now
        if USE_HANDBOOK:
            SPM_SCALE["F"] *= r
        else:
            for k in (2, 4, 6):
                SPM_DEFAULT["F"][f"B{k}"] *= r
        out["F"] = dict(scale=SPM_SCALE["F"], S_ref=S_ref)
    SPM_SCALE["S"] = SPM_SCALE["O"]
    out["S"] = dict(scale=SPM_SCALE["S"], note="O scale applied (no data)")
    out["source"] = "Newman & Ng handbook per-ion parameters" if USE_HANDBOOK else "generic shape"
    if verbose:
        print("SPM calibration:", out, flush=True)
    return out


# Er3+:LaF3 crystal field (Wybourne, cm^-1), Carnall et al. 1989 via qlanth data/carnall.py
CARNALL_ER_LAF3_BKQ = {(2, 0): -238.0, (4, 0): 453.0, (6, 0): 373.0, (2, 2): -91.0, (4, 2): 308.0,
                       (4, 4): 417.0, (6, 2): -489.0, (6, 4): -240.0, (6, 6): -536.0}


def apply_site_refs():
    """Measured first-shell statistics (site_refs.SITES) and RT homogeneous width."""
    used = {}
    try:
        from ..refs import site_refs as S
    except ImportError:
        return used
    v = lambda x: x.get("value") if isinstance(x, dict) else x
    # silica: CN 6.0, Er-O 2.28 A, spread 0.09 A (EXAFS, Marcus & Polman 1991)
    # aluminosilicate: CN 6.3, Er-O 2.21-2.25 A (Peters & Houde-Walter 1998)
    # phosphate: CN 6.4-7, Er-O 2.23-2.29 A (Karabulut 2002; Mountjoy 2007)
    meas = {"silica": dict(CN=6.0, R=2.28, R_std=0.09), "aluminosilicate": dict(CN=6.3, R=2.23),
            "phosphate": dict(CN=6.7, R=2.26)}
    for fam, d in meas.items():
        FAMILIES[fam].update(d); FAMILIES[fam]["site_src"] = "site_refs (EXAFS)"
        used[fam] = d
    # RT homogeneous width: Er 1.53 um in Al-silica EDFA, 11.5 nm ~ 49 cm^-1 (IEEE PTL 2(4) 1990)
    for fam in FAMILIES:
        FAMILIES[fam].setdefault("gamma_h", 49.0)
    used["gamma_h_cm1"] = 49.0
    return used
