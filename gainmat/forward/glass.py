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
    "silica":          dict(n=1.444, hw_max=1100, gap_family="silicate", ligand="O", CN=6.0, R=2.25, R_std=0.07, ang_std=12, loss_dB_m=0.005),
    "aluminosilicate": dict(n=1.46, hw_max=1100, gap_family="silicate", ligand="O", CN=6.5, R=2.27, R_std=0.09, ang_std=15, loss_dB_m=0.01),
    "germanosilicate": dict(n=1.46, hw_max=1000, gap_family="silicate", ligand="O", CN=6.0, R=2.26, R_std=0.08, ang_std=13, loss_dB_m=0.01),
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


# pump intensities: up to ~1 W in a ~50 um^2 single-mode core (2e6 W/cm^2)
I_PUMP = np.logspace(2, np.log10(2e6), 10)


def omegas_for(ion, family):
    """Measured Omega for (ion, family); fall back to the same ion in the nearest family,
    then to the same family for another ion (flagged)."""
    if (ion, family) in OMEGA:
        return OMEGA[(ion, family)], "measured"
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
    fam = copy.deepcopy(FAMILIES[family]); fam.update(overrides or {})
    om, om_src = omegas_for(ion, family) if "omegas" not in fam else (fam["omegas"], "override")
    lig = fam["ligand"]
    sites = fam.get("sites") or random_sites(n_sites, lig, fam["CN"], fam["R"], fam["R_std"],
                                             fam["ang_std"], seed=seed)
    spm = {lig: SPM_DEFAULT[lig]}
    emax = fam.get("emax", EMAX[ion])
    spec = ensemble_spectra(ion, sites, spm, om, fam["n"], nu=nu_grid(lam[0], lam[1], dnu),
                            emax=emax, T=fam.get("T", 295.0), gamma_h=fam.get("gamma_h", 20.0))
    N = fam.get("N_cm3", N_DEFAULT[ion])
    res = dict(ion=ion, family=family, omegas=om, omega_src=om_src, n=fam["n"], hw_max=fam["hw_max"],
               loss_dB_m=fam["loss_dB_m"], N_cm3=N, pumps={}, split=spec["split"])
    env_all = None
    loss = fam["loss_dB_m"] + ir_loss_dB_m(1e7 / spec["nu"], fam["hw_max"])
    for p in (pumps or PUMPS[ion]):
        try:
            gs = gain_spectra(ion, spec, om, fam["n"], fam["hw_max"], fam.get("gap_family"),
                              N_tot_cm3=N, pump_nm=p, intensities_W_cm2=intensities, emax=emax)
        except ValueError as e:
            res["pumps"][p] = dict(error=str(e)); continue
        env = envelope(gs)
        env_all = env if env_all is None else np.maximum(env_all, env)
        res["pumps"][p] = dict(pump_level=gs["pump_level"],
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
