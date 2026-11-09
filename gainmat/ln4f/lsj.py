"""Free-ion 4f^N solver in the coupled |f^N tau S L J> basis (Carnall/Crosswhite style).

Complete Carnall effective-operator Hamiltonian:
  F^2,F^4,F^6 ; zeta ; alpha L(L+1), beta G(G2), gamma G(R7) ;
  T^2,T^3,T^4,T^6,T^7,T^8 (three-body) ; M^0,M^2,M^4 (spin-spin + spin-other-orbit) ;
  P^2,P^4,P^6 (electrostatically-correlated spin-orbit).

Matrix elements are NOT computed here: they are read from the machine-readable LS(J)
tables of the qlanth project (zia-lab/qlanth, commit 14563fd, data/*.m; arXiv:2503.17377),
converted to numeric JSON (data/qlanth_f1to4.json.gz) by a mechanical parse of the
Mathematica expressions (script reproduced in tests/docs; no hand typing). The qlanth
tables are built from Nielson-Koster CFPs and reproduce the Crosswhite fNcross tables.
Labels |tau S L> use the Nielson-Koster (W)(U) classification; Casimir eigenvalues
G(G2)=U-label, G(R7)=W-label (Wybourne tables 2-6/2-7, as in qlanth).

Hole configurations N>7 use N'=14-N with zeta and all T^i sign-flipped and the
"t2Switch=0" form of T^2 (qlanth HoleElectronConjugation convention; M^k, P^k unchanged).

U^(k) reduced matrix elements between intermediate-coupled J multiplets are computed
directly from the tabulated doubly-reduced <tau S L||U^(k)||tau' S L'> (Nielson-Koster
normalisation <f||u^k||f>=1), and the magnetic-dipole <J||L+2S||J'> analytically.
"""
from __future__ import annotations
import gzip
import json
import os
import re
from functools import lru_cache
import numpy as np
from sympy.physics.wigner import wigner_6j
from sympy import Rational

_DATA = os.path.join(os.path.dirname(__file__), "data", "qlanth_f1to4.json.gz")
LETTERS = "SPDFGHIKLMNOQ"
GG2U = {"00": 0, "10": 6 / 12, "11": 12 / 12, "20": 14 / 12, "21": 21 / 12, "22": 30 / 12,
        "30": 24 / 12, "31": 32 / 12, "40": 36 / 12}
GSO7W = {"000": 0, "100": 3 / 5, "110": 1.0, "111": 6 / 5, "200": 7 / 5, "210": 9 / 5,
         "211": 2.0, "220": 12 / 5, "221": 13 / 5, "222": 3.0}
T_KEYS = ("T2", "T3", "T4", "T6", "T7", "T8")


@lru_cache(maxsize=1)
def _tables():
    with gzip.open(_DATA, "rt") as f:
        return json.load(f)


def _hf(x):
    return Rational(int(round(2 * x)), 2)


@lru_cache(maxsize=None)
def sixj(a, b, c, d, e, f):
    return float(wigner_6j(_hf(a), _hf(b), _hf(c), _hf(d), _hf(e), _hf(f)))


def term_SL(t):
    m = re.match(r"(\d)([A-Z])", t)
    return (int(m.group(1)) - 1) / 2.0, LETTERS.index(m.group(2))


def base_term(t):
    """'2H1' -> '2H' (Carnall-style label without the tau index)."""
    return t[:2]


@lru_cache(maxsize=None)
def terms(n):
    d = _tables()["electrostatic"]
    out = []
    for k in d:
        nn, a, b = k.split("|")
        if int(nn) == n and a == b:
            out.append(a)
    # any term with all-zero electrostatic (only f1) -> fall back to spin-orbit keys
    if not out:
        out = sorted({k.split("|")[1] for k in _tables()["spin_orbit"] if int(k.split("|")[0]) == n})
    return sorted(out, key=lambda t: (term_SL(t)[0], term_SL(t)[1], t), reverse=False)


@lru_cache(maxsize=None)
def nk_labels(t):
    """(W, U) Nielson-Koster labels for term t, first match as in qlanth FindNKLSTerm."""
    for lab, v, W, U in _tables()["terms_NK"]:
        if t in lab:
            return W, U
    raise KeyError(t)


def _get(tab, n, a, b, J=None):
    k = f"{n}|{a}|{b}" + ("" if J is None else f"|{float(J)}")
    v = tab.get(k)
    if v is None:
        k = f"{n}|{b}|{a}" + ("" if J is None else f"|{float(J)}")
        v = tab.get(k)
    return v


def _dot(v, p):
    if v is None:
        return 0.0
    if isinstance(v, (int, float)):
        return float(v)
    return sum(c * p.get(s, 0.0) for s, c in v.items())


def carnall_params(p):
    """Complete the Carnall parameter set: M2=0.56 M0, M4=0.31 M0, P4=0.5 P2, P6=0.1 P2
    unless given explicitly (Carnall 1989 Table 4 footnotes)."""
    q = dict(p)
    M0, P2 = q.get("M0", 0.0), q.get("P2", 0.0)
    q.setdefault("M2", 0.56 * M0)
    q.setdefault("M4", 0.31 * M0)
    q.setdefault("P4", 0.5 * P2)
    q.setdefault("P6", 0.1 * P2)
    return q


def J_values(n):
    Js = set()
    for t in terms(n):
        S, L = term_SL(t)
        J = abs(L - S)
        while J <= L + S + 1e-9:
            Js.add(J)
            J += 1
    return sorted(Js)


def J_block(N, J, p):
    """Hamiltonian block for total J of f^N: returns (list of terms, matrix)."""
    tb = _tables()
    hole = N > 7
    n = 14 - N if hole else N
    p = carnall_params(p)
    pp = dict(F0=0.0, F2=p["F2"], F4=p["F4"], F6=p["F6"], zeta=p["zeta"] * (-1 if hole else 1),
              M0=p["M0"], M2=p["M2"], M4=p["M4"], P2=p["P2"], P4=p["P4"], P6=p["P6"])
    for k in T_KEYS:
        pp[k] = p.get(k, 0.0) * (-1 if hole else 1)
    tbt = tb["three_body_t2switch0" if hole else "three_body_t2switch1"]
    ts = [t for t in terms(n) if abs(term_SL(t)[1] - term_SL(t)[0]) <= J + 1e-9 <= sum(term_SL(t)) + 1e-9]
    d = len(ts)
    H = np.zeros((d, d))
    for i, a in enumerate(ts):
        for j, b in enumerate(ts):
            v = _dot(_get(tb["electrostatic"], n, a, b), pp)
            v += _dot(_get(tb["spin_orbit"], n, a, b, J), pp)
            v += _dot(_get(tb["soo_ecso"], n, a, b, J), pp)
            v += _dot(_get(tb["spin_spin"], n, a, b, J), pp)
            v += _dot(_get(tbt, n, a, b), pp)
            if i == j and n > 1:
                S, L = term_SL(a)
                W, U = nk_labels(a)
                v += p.get("alpha", 0.0) * L * (L + 1) + p.get("beta", 0.0) * GG2U[U] \
                    + p.get("gamma", 0.0) * GSO7W[W]
            H[i, j] = v
    return ts, 0.5 * (H + H.T)


def jlabel(J):
    return f"{int(round(2 * J))}/2" if abs(J - round(J)) > 0.25 else f"{int(round(J))}"


class LSJResult:
    """Free-ion eigen-solution in the |tau S L J> basis; multiplets compatible with run_ln."""

    def __init__(self, N, multiplets):
        self.N, self.multiplets = N, multiplets
        self.n = 14 - N if N > 7 else N
        self.basis = None
        self.stark = None

    def table(self, emax=None):
        rows = [(m["label"], m["energy"], m["J"], m["purity"]) for m in self.multiplets]
        return [r for r in rows if emax is None or r[1] <= emax]

    # ---- intensities
    def U2(self, i, j, k):
        """|<i||U^(k)||j>|^2 between multiplets i, j (indices into self.multiplets)."""
        return reduced_U(self.n, self.multiplets[i], self.multiplets[j], k) ** 2

    def MD(self, i, j):
        return reduced_LS(self.multiplets[i], self.multiplets[j]) ** 2

    def reduced_elements(self, emax=None):
        m = self.multiplets
        idx = [i for i in range(len(m)) if emax is None or m[i]["energy"] <= emax]
        out = {}
        for i in idx:
            for j in idx:
                if m[j]["energy"] >= m[i]["energy"]:
                    continue
                d = {f"U{l}": self.U2(i, j, l) for l in (2, 4, 6)}
                d["MD"] = self.MD(i, j) if abs(m[i]["J"] - m[j]["J"]) <= 1 else 0.0
                out[(i, j)] = d
        return out


def reduced_U(n, mi, mj, k):
    tab = _tables()["Uk"]
    Ji, Jj = mi["J"], mj["J"]
    tot = 0.0
    for a, ca in zip(mi["terms"], mi["vec"]):
        Sa, La = term_SL(a)
        for b, cb in zip(mj["terms"], mj["vec"]):
            Sb, Lb = term_SL(b)
            if Sa != Sb or ca == 0 or cb == 0:
                continue
            r = _get(tab, n, a, b, k)
            if r is None:
                continue
            # tables store <a||U^k||b>; if only (b,a) present use the symmetric relation
            if f"{n}|{a}|{b}|{float(k)}" not in tab:
                r = r * (-1) ** (La - Lb)
            r = float(r)
            f = (-1) ** int(round(Sa + Lb + Ji + k)) * np.sqrt((2 * Ji + 1) * (2 * Jj + 1)) \
                * sixj(La, Ji, Sa, Jj, Lb, k)
            tot += ca * cb * f * r
    return tot


def reduced_LS(mi, mj):
    """<i||L+2S||j> (units of hbar) = <||J||> + <||S||>, diagonal in tau S L."""
    Ji, Jj = mi["J"], mj["J"]
    tot = 0.0
    idx = {t: c for t, c in zip(mj["terms"], mj["vec"])}
    for a, ca in zip(mi["terms"], mi["vec"]):
        cb = idx.get(a)
        if cb is None:
            continue
        S, L = term_SL(a)
        v = 0.0
        if abs(Ji - Jj) < 1e-9:
            v += np.sqrt(Ji * (Ji + 1) * (2 * Ji + 1))
        v += (-1) ** int(round(S + L + Ji + 1)) * np.sqrt((2 * Ji + 1) * (2 * Jj + 1)) \
            * sixj(S, Ji, L, Jj, S, 1) * np.sqrt(S * (S + 1) * (2 * S + 1))
        tot += ca * cb * v
    return tot


def solve_lsj(N, p, shift_ground=True):
    """Diagonalise the complete Carnall free-ion Hamiltonian, J block by J block.

    Each eigenvector is labelled by the SL term with the largest weight summed over its
    tau copies (Carnall convention: printed without tau, e.g. '2H11/2'); within a J block
    eigenvectors are taken in increasing energy and a label already used falls through to
    the next-largest SL, so labels are unique. Returns LSJResult (multiplets sorted by energy).
    """
    n = 14 - N if N > 7 else N
    mults = []
    for J in J_values(n):
        ts, H = J_block(N, J, p)
        E, V = np.linalg.eigh(H)
        used = set()
        bases = sorted({base_term(t) for t in ts})
        for c in range(len(E)):
            w = V[:, c] ** 2
            # Carnall convention: weight of an SL term summed over its tau (seniority/U) copies
            agg = {bt: float(sum(w[i] for i, t in enumerate(ts) if base_term(t) == bt)) for bt in bases}
            order_b = sorted(agg, key=lambda bt: -agg[bt])
            bt = next((b for b in order_b if b + jlabel(J) not in used), order_b[0])
            used.add(bt + jlabel(J))
            lab = bt + jlabel(J)
            comp = [i for i, t in enumerate(ts) if base_term(t) == bt]
            pick = max(comp, key=lambda i: w[i])
            order = np.argsort(-w)
            S, L = term_SL(ts[pick])
            mults.append(dict(label=lab, term=bt, tau_term=ts[pick], S=S, L=L, J=J,
                              energy=float(E[c]), purity=agg[bt], terms=ts, vec=V[:, c].copy(),
                              composition={ts[o]: float(w[o]) for o in order[:3]}))
    # labels unique per J already; resolve residual duplicates (same label, different J impossible)
    mults.sort(key=lambda m: m["energy"])
    seen = {}
    for m in mults:   # later copies of a label (another tau of the same SL) -> '2H9/2(2)'
        k = m["label"]
        seen[k] = seen.get(k, 0) + 1
        if seen[k] > 1:
            m["label"] = f"{k}({seen[k]})"
    if shift_ground:
        e0 = mults[0]["energy"]
        for m in mults:
            m["energy"] -= e0
    return LSJResult(N, mults)


def duplicate_labels(res, emax=None):
    labs = [m["label"] for m in res.multiplets if emax is None or m["energy"] <= emax]
    return len(labs) - len(set(labs))
