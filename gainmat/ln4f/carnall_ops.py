"""Higher-order Carnall effective operators in the 4f^N Slater-determinant basis.

Operators added to the free-ion Hamiltonian of fn_ci (Carnall et al. 1989 parametrisation):

  * three-body  sum_i T^i t_i, i = 2,3,4,6,7,8   (Judd 1966; electron picture, original t2)
  * spin-other-orbit (+ spin-spin)   sum_k M^k m_k, k = 0,2,4   (Judd, Crosswhite, Crosswhite 1968)
  * electrostatically-correlated spin-orbit   sum_k P^k p_k, k = 2,4,6   (same source)

Construction (no tables beyond f^2/f^3 are used):
  1. The f^2 (M, P) or f^3 (T) LS/LSJ matrix elements (carnall_ops_data, extracted from the qlanth
     tables of the published values) are turned into a matrix in the 91-dim (f^2) or 364-dim (f^3)
     determinant space using explicitly constructed coupled states. f^2 states are the coupled
     two-electron states [l x l]^L [s x s]^S; their relative signs are fixed against the tabulated
     spin-orbit matrix (verified: all f^2 spin-orbit elements reproduced). f^3 states are built from
     the Nielson-Koster CFPs, |f3 tau SL> ~ sum cfp [a^+ x |f2 SL_bar>]^{SL}; this fixes the
     seniority-split pairs 2D1/2D2 ... 2H1/2H2 (verified against the tabulated f^3 Coulomb and
     spin-orbit matrices, see tests).
  2. The f^2 (f^3) determinant matrix *is* the antisymmetrised two- (three-) body integral table;
     the operator is then applied as a genuine two-/three-body operator in the f^N space. For a pure
     k-body operator this is exactly Judd's CFP recursion n/(n-k) sum cfp cfp <f^{n-1}|.|f^{n-1}>.
  3. Beyond half filling Carnall/Crosswhite (and qlanth, which reproduces Carnall's computed levels)
     take the M^k and P^k matrices of f^N equal to those of f^{14-N} (no sign change). The true
     electron-picture two-body operator differs from that by a one-body l.s-like piece (an effective
     zeta renormalisation that the fitted zeta absorbs). We therefore apply M/P for N>7 through an
     explicit particle-hole map W (rotation covariant, checked on L, S) of the f^{14-N} operator.
     The t_i (i != 2) satisfy t_i(f^{14-N}) = -t_i(f^N) automatically for the electron-picture
     three-body operator, and t_2 is taken in the electron picture as in Carnall 1989/qlanth
     (t2Switch = 0 for N>7 together with T-sign flip == electron-picture t2).
"""
from __future__ import annotations
from functools import lru_cache
from itertools import combinations
import numpy as np
from sympy.physics.wigner import clebsch_gordan
from . import carnall_ops_data as D

L_F = 3
NSO = 14
LETTERS = "SPDFGHIKLMNOQ"


@lru_cache(maxsize=None)
def cg(j1, m1, j2, m2, j, m):
    from sympy import Rational as R
    f = lambda x: R(int(round(2 * x)), 2)
    return float(clebsch_gordan(f(j1), f(j2), f(j), f(m1), f(m2), f(m)))


def orb(m, ms2):
    """spin-orbital index for m_l=m, 2 m_s = ms2 (same convention as fn_ci)."""
    return 2 * (m + 3) + (0 if ms2 == 1 else 1)


def parse_term(t):
    """'2D1' -> (mult=2, L=2, tag='2D1')."""
    return int(t[0]), LETTERS.index(t[1])


def _norm(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else v


# ---------------------------------------------------------------- f^2 coupled states
def f2_state(S, L, MS, ML, basis):
    """|f2 S L MS ML> = sum CG(3 m1 3 m2|L ML) CG(1/2 s1 1/2 s2|S MS) a+_1 a+_2 |0>, normalised."""
    v = np.zeros(basis.dim)
    for m1 in range(-3, 4):
        m2 = ML - m1
        if abs(m2) > 3:
            continue
        cl = cg(3, m1, 3, m2, L, ML)
        if cl == 0:
            continue
        for s1 in (1, -1):
            s2 = 2 * MS - s1
            if s2 not in (1, -1):
                continue
            cs = cg(0.5, s1 / 2, 0.5, s2 / 2, S, MS)
            if cs == 0:
                continue
            p, q = orb(m1, s1), orb(m2, s2)
            if p == q:
                continue
            # a+_p a+_q |0>: = +det if p<q else -det
            sign = 1.0 if p < q else -1.0
            v[basis.index[(1 << p) | (1 << q)]] += sign * cl * cs
    return _norm(v)


F2_TERMS = ["3P", "3F", "3H", "1S", "1D", "1G", "1I"]


def _apply_create(basis_from, basis_to, v, p):
    """a+_p acting on vector v of basis_from -> vector in basis_to."""
    out = np.zeros(basis_to.dim)
    for j, d in enumerate(basis_from.dets):
        if v[j] == 0 or (d >> p) & 1:
            continue
        s = -1 if bin(d & ((1 << p) - 1)).count("1") % 2 else 1
        out[basis_to.index[d | (1 << p)]] += s * v[j]
    return out


@lru_cache(maxsize=None)
def f2_signs():
    """Per-term signs making the coupled f^2 states reproduce the tabulated spin-orbit matrix."""
    from .fn_ci import get_basis
    b = get_basis(2)
    so = b.angular()["so"]
    signs = {"3P": 1}
    # resolve signs through connected off-diagonal SO elements
    changed = True
    while changed:
        changed = False
        for key, val in D.SO2.items():
            a, c, J = key.split("|")
            if a == c or abs(val) < 1e-9:
                continue
            if (a in signs) == (c in signs):
                continue
            J = float(eval_frac(J))
            va = f2_lsj(a, J, J, b)
            vc = f2_lsj(c, J, J, b)
            mine = va @ so @ vc
            rel = np.sign(mine) * np.sign(val)
            if a in signs:
                signs[c] = int(signs[a] * rel)
            else:
                signs[a] = int(signs[c] * rel)
            changed = True
    return signs


def eval_frac(s):
    if "/" in s:
        a, b = s.split("/")
        return int(a) / int(b)
    return float(s)


def f2_lsj(term, J, MJ, b, sign=1):
    mult, L = parse_term(term)
    S = (mult - 1) / 2
    v = np.zeros(b.dim)
    for MS2 in range(-int(2 * S), int(2 * S) + 1, 2):
        MS = MS2 / 2
        ML = MJ - MS
        if abs(ML) > L + 1e-9 or abs(ML - round(ML)) > 1e-9:
            continue
        c = cg(S, MS, L, int(round(ML)), J, MJ)
        if c:
            v += c * f2_state(S, L, MS, int(round(ML)), b)
    return sign * v


def lsj_operator_f2(table):
    """91x91 determinant matrix of an f^2 J-scalar operator given by table 'SL|S'L'|J' -> value."""
    from .fn_ci import get_basis
    b = get_basis(2)
    sg = f2_signs()
    X = np.zeros((b.dim, b.dim))
    vec = {}
    for key, val in table.items():
        a, c, Js = key.split("|")
        J = eval_frac(Js)
        if abs(val) < 1e-15:
            continue
        for MJ2 in range(-int(2 * J), int(2 * J) + 1, 2):
            MJ = MJ2 / 2
            for t in (a, c):
                if (t, J, MJ) not in vec:
                    vec[(t, J, MJ)] = f2_lsj(t, J, MJ, b, sg[t])
            X += val * np.outer(vec[(a, J, MJ)], vec[(c, J, MJ)])
    return X


# ---------------------------------------------------------------- f^3 CFP states
def f3_state(term, MS, ML, b2, b3):
    """|f3 tau S L MS ML> from Nielson-Koster CFPs: sum_parent cfp [ |f2 parent> x f ]^{S L}."""
    mult, L = parse_term(term)
    S = (mult - 1) / 2
    sg = f2_signs()
    v = np.zeros(b3.dim)
    for par, c in D.CFP3[term].items():
        pm, pL = parse_term(par)
        pS = (pm - 1) / 2
        for pMS2 in range(-int(2 * pS), int(2 * pS) + 1, 2):
            pMS = pMS2 / 2
            s = MS - pMS
            if abs(abs(s) - 0.5) > 1e-9:
                continue
            cs = cg(pS, pMS, 0.5, s, S, MS)
            if cs == 0:
                continue
            for pML in range(-pL, pL + 1):
                m = ML - pML
                if abs(m) > 3:
                    continue
                cl = cg(pL, pML, 3, m, L, ML)
                if cl == 0:
                    continue
                par_v = sg[par] * f2_state(pS, pL, pMS, pML, b2)
                v += c * cs * cl * _apply_create(b2, b3, par_v, orb(m, int(round(2 * s))))
    return _norm(v)


def f3_lsj(term, J, MJ, b2, b3):
    mult, L = parse_term(term)
    S = (mult - 1) / 2
    v = np.zeros(b3.dim)
    for MS2 in range(-int(2 * S), int(2 * S) + 1, 2):
        MS = MS2 / 2
        ML = MJ - MS
        if abs(ML) > L + 1e-9:
            continue
        c = cg(S, MS, L, int(round(ML)), J, MJ)
        if c:
            v += c * f3_state(term, MS, int(round(ML)), b2, b3)
    return v


@lru_cache(maxsize=None)
def f3_states():
    """{(term, MS, ML): vector} for all f^3 terms."""
    from .fn_ci import get_basis
    b2, b3 = get_basis(2), get_basis(3)
    out = {}
    for term in D.CFP3:
        mult, L = parse_term(term)
        S = (mult - 1) / 2
        for MS2 in range(-int(2 * S), int(2 * S) + 1, 2):
            for ML in range(-L, L + 1):
                out[(term, MS2 / 2, ML)] = f3_state(term, MS2 / 2, ML, b2, b3)
    return out


TI = ("T2", "T3", "T4", "T6", "T7", "T8")


@lru_cache(maxsize=None)
def t_operator_f3(name):
    """364x364 determinant matrix of t_i in f^3 (SL-scalar, J-independent)."""
    from .fn_ci import get_basis
    b3 = get_basis(3)
    st = f3_states()
    X = np.zeros((b3.dim, b3.dim))
    for key, d in D.T3.items():
        val = d.get(name, 0.0)
        if abs(val) < 1e-15:
            continue
        a, c = key.split("|")
        mult, L = parse_term(a)
        S = (mult - 1) / 2
        for MS2 in range(-int(2 * S), int(2 * S) + 1, 2):
            for ML in range(-L, L + 1):
                X += val * np.outer(st[(a, MS2 / 2, ML)], st[(c, MS2 / 2, ML)])
    return X


# ---------------------------------------------------------------- generic k-body appliers
def _sign(det, p):
    return -1 if bin(det & ((1 << p) - 1)).count("1") % 2 else 1


def apply_kbody(basis, Xk, k):
    """Matrix in the N-electron basis of the k-body operator whose f^k determinant matrix is Xk:
    O = sum_{P,R} Xk[P,R] a+_{p1}..a+_{pk} a_{rk}..a_{r1}  (p1<..<pk, r1<..<rk)."""
    from .fn_ci import get_basis
    bk = get_basis(k)
    N = basis.N
    H = np.zeros((basis.dim, basis.dim))
    if N < k:
        return H
    occs = [tuple(p for p in range(NSO) if (d >> p) & 1) for d in bk.dets]
    by_r = {}
    nz = np.argwhere(np.abs(Xk) > 1e-13)
    for i, j in nz:
        by_r.setdefault(occs[j], []).append((occs[i], Xk[i, j]))
    for j, d in enumerate(basis.dets):
        occ = [p for p in range(NSO) if (d >> p) & 1]
        for R in combinations(occ, k):
            lst = by_r.get(R)
            if not lst:
                continue
            sg, dd = 1, d
            for r in R:                      # a_{r1} acts first
                sg *= _sign(dd, r); dd ^= (1 << r)
            for P, v in lst:
                ok, s2, d2 = True, sg, dd
                for p in reversed(P):        # a+_{pk} acts first
                    if (d2 >> p) & 1:
                        ok = False; break
                    s2 *= _sign(d2, p); d2 |= (1 << p)
                if ok:
                    H[basis.index[d2], j] += s2 * v
    return H


# ---------------------------------------------------------------- particle-hole map
def _tr(p):
    """time-reversed orbital (m, ms) -> (-m, -ms) and phase (-1)^(l-m+1/2-ms)."""
    m, ms2 = p // 2 - 3, (1 if p % 2 == 0 else -1)
    q = orb(-m, -ms2)
    ph = (-1) ** int(round(L_F - m + 0.5 - ms2 / 2))
    return q, ph


def hole_map(N):
    """W (dim_{14-N} x dim_N): N-electron det -> (14-N)-'electron' state of its holes, such that
    W X^{(14-N)} W^T (as an operator on f^N) is rotation-covariant (checked on L, S in tests)."""
    from .fn_ci import get_basis
    bN, bH = get_basis(N), get_basis(NSO - N)
    W = np.zeros((bH.dim, bN.dim))
    full = (1 << NSO) - 1
    for j, d in enumerate(bN.dets):
        holes = [p for p in range(NSO) if not (d >> p) & 1]
        # hole state  prod_{h in holes, ascending} a_h |full>  ~  prod tilde-creators |0>
        # sign of removing holes from the filled shell
        sg, dd = 1, full
        for h in holes:
            sg *= _sign(dd, h); dd ^= (1 << h)
        # map each hole h to electron in tr(h) with phase
        v = np.zeros(bH.dim); state = {0: float(sg)}
        for h in holes:
            q, ph = _tr(h)
            new = {}
            for dt, a in state.items():
                if (dt >> q) & 1:
                    continue
                new[dt | (1 << q)] = new.get(dt | (1 << q), 0) + a * ph * _sign(dt, q)
            state = new
        for dt, a in state.items():
            W[bH.index[dt], j] = a
    return W


# ---------------------------------------------------------------- public builders
@lru_cache(maxsize=None)
def _mp_tables():
    keys = ("M0", "M2", "M4", "P2", "P4", "P6")
    tabs = {k: {} for k in keys}
    for key, d in D.MP2.items():
        for k, v in d.items():
            tabs[k][key] = v
    ss = {k: {} for k in ("M0", "M2", "M4")}
    for key, d in D.SS2.items():
        for k, v in d.items():
            ss[k][key] = v
    return tabs, ss


@lru_cache(maxsize=None)
def f2_operator(name, spin_spin=True):
    tabs, ss = _mp_tables()
    X = lsj_operator_f2(tabs[name])
    if spin_spin and name in ss:
        X = X + lsj_operator_f2(ss[name])
    return X


def magnetic_operators(basis, spin_spin=True):
    """{'M0','M2','M4','P2','P4','P6'} matrices in the f^N determinant basis."""
    key = ("MP", spin_spin)
    if key in basis._cache:
        return basis._cache[key]
    from .fn_ci import get_basis
    N = basis.N
    out = {}
    names = ("M0", "M2", "M4", "P2", "P4", "P6")
    if N <= 7:
        for n in names:
            out[n] = apply_kbody(basis, f2_operator(n, spin_spin), 2)
    else:
        bH = get_basis(NSO - N)
        W = hole_map(N)
        for n in names:
            XH = apply_kbody(bH, f2_operator(n, spin_spin), 2)
            out[n] = W.T @ XH @ W
    basis._cache[key] = out
    return out


def three_body_operators(basis):
    """{'T2',...,'T8'} matrices (electron picture) in the f^N determinant basis."""
    if "T" in basis._cache:
        return basis._cache["T"]
    out = {n: apply_kbody(basis, t_operator_f3(n), 3) for n in TI}
    basis._cache["T"] = out
    return out
