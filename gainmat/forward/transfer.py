"""Energy-transfer microparameters from the model's own cross sections (Forster-Dexter, dipole-dipole,
with phonon assistance), so that cross-relaxation, migration and upconversion strengths can be
predicted for a glass where none were measured.

    C_DA = 3c / (8 pi^4 n^2) * exp(-(2 nbar + 1) S0) * sum_m S0^m (nbar + 1)^m / m!
                 * integral sigma_em^D(lambda) sigma_abs^A(lambda_m^+) d lambda,   1/lambda_m^+ = 1/lambda - m h w

(Kushida 1973; phonon-assisted form as used by Tarelho, Gomes et al. and by most glass papers that
report C_DA, e.g. Han et al., Sci. Rep. 6, 31207 (2016); Huang et al., Sci. Rep. 4, 5053 (2014)).
sigma in cm^2, lambda in cm -> C_DA in cm^6 s^-1.  m = 0 is the resonant term.

    from gainmat.forward.transfer import c_da
    c_da("Tm3+", "germanate", donor=("3H4", "3F4"), acceptor=("3F4", "3H6"))   # 3H4,3H6 -> 3F4,3F4
"""
from __future__ import annotations

import math

import numpy as np

from . import glass as g

C_CGS = 2.99792458e10
KB_CM = 0.6950348


def c_da(ion, family, donor, acceptor, S0=0.31, hw=None, m_max=3, T=295.0, n_sites=60, overrides=None,
         ion_acceptor=None, spec=None, spec_acceptor=None):
    """donor = (upper, lower) emitting transition of the donor; acceptor = (upper, lower) transition
    whose ABSORPTION (lower -> upper) the acceptor makes.  Returns dict(total, per_m, fractions)."""
    if spec is None:
        spec = g.run(ion, family, n_sites=n_sites, pumps=[], overrides=overrides)["_spec"]
    ia = ion_acceptor or ion
    if spec_acceptor is None:
        spec_acceptor = spec if ia == ion else g.run(ia, family, n_sites=n_sites, pumps=[], overrides=overrides)["_spec"]
    F = g.FAMILIES[family]
    n = (overrides or {}).get("n", F["n"])
    hw = hw if hw is not None else (overrides or {}).get("hw_max", F["hw_max"])
    nu = spec["nu"]; se = spec["sigma_em"][donor]
    nuA = spec_acceptor["nu"]; sa = spec_acceptor["sigma_abs"][acceptor]
    nbar = 1.0 / (math.exp(hw / (KB_CM * T)) - 1.0)
    pref = 3 * C_CGS / (8 * math.pi ** 4 * n ** 2)
    per_m = []
    for m in range(m_max + 1):
        w = math.exp(-(2 * nbar + 1) * S0) * S0 ** m * (nbar + 1) ** m / math.factorial(m)
        # the donor's surplus energy goes into m emitted phonons: acceptor absorption at nu_A = nu - m hw
        sa_m = np.interp(nu - m * hw, nuA, sa, left=0.0, right=0.0)
        # integral over lambda (cm): d lambda = d nu / nu^2 (nu in cm^-1, lambda = 1/nu cm)
        integ = float(np.trapezoid(se * sa_m / nu ** 2, nu)) if hasattr(np, "trapezoid") else float(np.trapz(se * sa_m / nu ** 2, nu))
        per_m.append(pref * w * abs(integ))
    tot = float(sum(per_m))
    return dict(total=tot, per_m=per_m, fractions=[p / tot if tot else 0.0 for p in per_m], n=n, hw=hw, S0=S0)


# Standard single-ion-species processes: (label, donor emission (u, l), acceptor absorption (u, l)).
PROCESSES = {
    "Er3+": [("4I13/2 migration", ("4I13/2", "4I15/2"), ("4I13/2", "4I15/2")),
             ("4I11/2 migration", ("4I11/2", "4I15/2"), ("4I11/2", "4I15/2")),
             ("upconversion 4I13/2+4I13/2 -> 4I9/2+4I15/2", ("4I13/2", "4I15/2"), ("4I9/2", "4I13/2"))],
    "Tm3+": [("3F4 migration", ("3F4", "3H6"), ("3F4", "3H6")),
             ("3H4 migration", ("3H4", "3H6"), ("3H4", "3H6")),
             ("cross-relaxation 3H4+3H6 -> 3F4+3F4", ("3H4", "3F4"), ("3F4", "3H6"))],
    "Ho3+": [("5I7 migration", ("5I7", "5I8"), ("5I7", "5I8"))],
    "Yb3+": [("2F5/2 migration", ("2F5/2", "2F7/2"), ("2F5/2", "2F7/2"))],
    "Nd3+": [("4F3/2 migration", ("4F3/2", "4I9/2"), ("4F3/2", "4I9/2")),
             ("self-quenching 4F3/2+4I9/2 -> 4I15/2+4I15/2", ("4F3/2", "4I15/2"), ("4I15/2", "4I9/2"))],
}


def standard_parameters(ion, family, spec=None, overrides=None):
    """C_DA (cm^6/s) for the standard processes of one ion in one glass, from the model's own spectra."""
    out = []
    for lab, d, a in PROCESSES.get(ion, []):
        try:
            r = c_da(ion, family, d, a, spec=spec, overrides=overrides)
        except KeyError:
            continue
        out.append(dict(process=lab, C_DA_cm6_s=r["total"], phonon_fractions=[round(f, 3) for f in r["fractions"]]))
    return out
