"""Concentration self-quenching of a level by cross-relaxation, predicted from the model's own spectra:
lifetime vs dopant concentration with no fitted parameter.

Donor decay (Inokuti-Hirayama static dipole-dipole transfer + Burshtein migration-assisted transfer):

    I(t) = exp(-t/tau0 - gamma sqrt(t) - W t)
    gamma = (4/3) pi^(3/2) N_A sqrt(C_DA)              (static, s^-1/2)
    W     = pi (2 pi/3)^(5/2) sqrt(C_DA C_DD) N_A N_D   (hopping, s^-1)

For self-quenching the acceptor is the same ion in its ground state, N_A = N_D = N.  C_DA (cross-
relaxation) and C_DD (donor-donor migration) come from gainmat.forward.transfer.c_da on the model's
cross sections; tau0 is the model's own low-concentration lifetime of the level (radiative +
multiphonon).  Reported lifetimes are the mean lifetime int I dt and the 1/e time.

Not included (so the prediction is a floor on quenching): clustering / ion pairs, OH- and other
impurity quenching, upconversion (needs excited acceptors), radiation trapping.

    from gainmat.forward.quenching import lifetime_vs_N
    lifetime_vs_N("Tm3+", "silica", [1e19, 1e20, 5e20])
"""
from __future__ import annotations

import math

import numpy as np

from . import glass as g
from .transfer import c_da

# level: (cross-relaxation channels [(donor emission (u,l), acceptor absorption (u,l))], migration transition)
SELF_QUENCH = {
    "Tm3+": {"3H4": ([(("3H4", "3F4"), ("3F4", "3H6"))], ("3H4", "3H6"))},
    "Nd3+": {"4F3/2": ([(("4F3/2", "4I15/2"), ("4I15/2", "4I9/2"))], ("4F3/2", "4I9/2"))},
    "Pr3+": {"1D2": ([(("1D2", "1G4"), ("3F4", "3H4")), (("1D2", "1G4"), ("3F3", "3H4"))], ("1D2", "3H4"))},
    "Ho3+": {"5S2": ([(("5S2", "5I4"), ("5I7", "5I8")), (("5S2", "5I7"), ("5I4", "5I8"))], ("5S2", "5I8"))},
}
K_HOP = math.pi * (2 * math.pi / 3) ** 2.5
K_STAT = (4.0 / 3.0) * math.pi ** 1.5


def microparameters(ion, family, level, spec=None, overrides=None, n_sites=60):
    chans, mig = SELF_QUENCH[ion][level]
    if spec is None:
        spec = g.run(ion, family, n_sites=n_sites, pumps=[], overrides=overrides)["_spec"]
    cda = 0.0
    for d, a in chans:
        try:
            cda += c_da(ion, family, d, a, spec=spec, overrides=overrides)["total"]
        except KeyError:
            continue
    cdd = c_da(ion, family, mig, mig, spec=spec, overrides=overrides)["total"]
    return cda, cdd


def decay_lifetimes(tau0, gamma, W, n=4000):
    """mean lifetime int I dt and 1/e time (s) of I(t) = exp(-t/tau0 - gamma sqrt t - W t)."""
    t = np.linspace(0.0, 8.0 * tau0, n)
    I = np.exp(-t / tau0 - gamma * np.sqrt(t) - W * t)
    mean = float(np.trapezoid(I, t)) if hasattr(np, "trapezoid") else float(np.trapz(I, t))
    t_e = float(np.interp(math.exp(-1.0), I[::-1], t[::-1]))     # I decreases monotonically
    return mean, t_e


def lifetime_vs_N(ion, family, N_list, level=None, overrides=None, spec=None, lm=None, n_sites=60, tau0_s=None):
    """Predicted lifetime (ms) of a self-quenching level at each concentration N (ions/cm^3).
    tau0_s: low-concentration lifetime to use instead of the model's (e.g. a measured one)."""
    level = level or next(iter(SELF_QUENCH[ion]))
    if spec is None or (lm is None and tau0_s is None):
        r = g.run(ion, family, n_sites=n_sites, pumps=[], overrides=overrides)
        spec = r["_spec"]
        if lm is None and tau0_s is None:
            from .fiber import LocalModel
            F = g.FAMILIES[family]
            n = (overrides or {}).get("n", F["n"])
            lm = LocalModel(ion, spec, r["omegas"], n, (overrides or {}).get("hw_max", F["hw_max"]),
                            F.get("gap_family"), emax=g.EMAX[ion])
    if tau0_s is None:
        i = lm.labels.index(level)
        tau0_s = 1.0 / (float(lm.A[i].sum()) + float(lm.W[i]))
    cda, cdd = microparameters(ion, family, level, spec=spec, overrides=overrides)
    out = []
    for N in N_list:
        gam = K_STAT * N * math.sqrt(cda)
        W = K_HOP * math.sqrt(cda * cdd) * N * N
        mean, te = decay_lifetimes(tau0_s, gam, W)
        out.append(dict(N_cm3=float(N), tau_mean_ms=1e3 * mean, tau_1e_ms=1e3 * te, gamma_s_half=gam, W_hop_s=W))
    return dict(ion=ion, family=family, level=level, tau0_ms=1e3 * tau0_s, C_CR_cm6_s=cda, C_mig_cm6_s=cdd, points=out)
