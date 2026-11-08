"""Judd-Ofelt intensities from intermediate-coupled 4f^N eigenvectors.

Squared reduced matrix elements |<aJ||U^(l)||bJ'>|^2 are obtained basis-independently:
by Wigner-Eckart, sum_{M,M',q} |<JM|T_q|J'M'>|^2 = |<J||T||J'>|^2 (sum of 3j^2 = 1),
so we take the Frobenius norm of V_J^+ U_q V_J' summed over q, where V_J spans the
(2J+1)-dimensional multiplet. Same for the magnetic dipole operator L+2S.
Conventions follow Carnall et al. tables (U^(l) = sum_i u^(l)(i), <f||u^l||f>=1).
"""
from __future__ import annotations
import numpy as np

# CGS constants
E_ESU = 4.80320471e-10
H_CGS = 6.62607015e-27
M_E = 9.1093837e-28
C_CGS = 2.99792458e10
MU_B_ESU = E_ESU * H_CGS / (4 * np.pi * M_E * C_CGS)   # e*hbar/(2mc)


def _rme2(ops, Va, Vb):
    return float(sum(np.sum(np.abs(Va.conj().T @ O @ Vb) ** 2) for O in ops))


def U2(basis, Va, Vb, lam):
    return _rme2([basis.U(lam, q) for q in range(-lam, lam + 1)], Va, Vb)


def md_ops(basis):
    a = basis.angular()
    Mz = a["Lz"] + 2 * a["Sz"]
    Mp = a["Lp"] + 2 * a["Sp"]
    return [Mz, -Mp / np.sqrt(2), Mp.T / np.sqrt(2)]


def LS2(basis, Va, Vb):
    """|<aJ||L+2S||bJ'>|^2 (units of hbar^2)."""
    return _rme2(md_ops(basis), Va, Vb)


def reduced_elements(res, emax=None):
    """Dict {(i_upper, i_lower): dict(U2,U4,U6,MD)} over multiplets of an FnResult."""
    m = res.multiplets
    idx = [i for i in range(len(m)) if emax is None or m[i]["energy"] <= emax]
    out = {}
    b = res.basis
    for i in idx:
        for j in idx:
            if m[j]["energy"] >= m[i]["energy"]:
                continue
            Vi, Vj = m[i]["vecs"], m[j]["vecs"]
            d = {f"U{l}": U2(b, Vi, Vj, l) for l in (2, 4, 6)}
            d["MD"] = LS2(b, Vi, Vj) if abs(m[i]["J"] - m[j]["J"]) <= 1 else 0.0
            out[(i, j)] = d
    return out


def chi_ed(n):
    return n * (n * n + 2) ** 2 / 9.0


def chi_md(n):
    return n ** 3


def rates(res, omegas, n_index, emax=30000.0, red=None):
    """Radiative rates between all multiplets below emax.

    omegas: (O2,O4,O6) in 1e-20 cm^2. Returns list of dicts and per-level summaries.
    A_ED = 64 pi^4 e^2 sigma^3/(3h(2J+1)) chi_ED sum Omega U^2
    A_MD = 64 pi^4 sigma^3/(3h(2J+1)) chi_MD mu_B^2 |<L+2S>|^2      (sigma in cm^-1)
    """
    m = res.multiplets
    red = red if red is not None else reduced_elements(res, emax)
    O = np.array(omegas, float) * 1e-20
    trans = []
    for (i, j), d in red.items():
        sig = m[i]["energy"] - m[j]["energy"]
        if sig <= 0:
            continue
        g = 2 * m[i]["J"] + 1
        pref = 64 * np.pi ** 4 * sig ** 3 / (3 * H_CGS * g)
        S_ed = E_ESU ** 2 * (O[0] * d["U2"] + O[1] * d["U4"] + O[2] * d["U6"])
        S_md = MU_B_ESU ** 2 * d["MD"]
        A_ed = pref * chi_ed(n_index) * S_ed
        A_md = pref * chi_md(n_index) * S_md
        trans.append(dict(i=i, j=j, upper=m[i]["label"], lower=m[j]["label"],
                          sigma_cm1=sig, wavelength_nm=1e7 / sig, A_ed=A_ed, A_md=A_md,
                          A=A_ed + A_md, S_ed_cm2=S_ed / E_ESU ** 2, **d))
    tot = {}
    for t in trans:
        tot[t["i"]] = tot.get(t["i"], 0.0) + t["A"]
    for t in trans:
        t["branching"] = t["A"] / tot[t["i"]] if tot[t["i"]] > 0 else 0.0
    tau = {i: (1.0 / a if a > 0 else np.inf) for i, a in tot.items()}
    return trans, tot, tau


def find(res, label, nth=0):
    """Index of the nth multiplet with a given label (ordered by energy)."""
    hits = [i for i, m in enumerate(res.multiplets) if m["label"] == label]
    if len(hits) <= nth:
        raise KeyError(label)
    return hits[nth]
