"""Explicit many-electron 4f^N Hamiltonian in a Slater-determinant basis.

Basis: all C(14,N) determinants built from the 14 spin-orbitals |l=3, m_l, m_s>.
Electrons are treated explicitly for every N (no particle-hole transformation is
needed: for N>7 the full electron-basis Hamiltonian automatically reproduces the
hole picture, e.g. the sign flip of the spin-orbit term).

H = sum_k F^k f_k (Coulomb, k=2,4,6)  + zeta * sum_i l_i.s_i
    + alpha L(L+1) + beta G(G2) + gamma G(R7)  + sum_kq B^k_q C^k_q   (cm^-1)

Casimir operators are built from the many-electron unit tensors
U^(k) = sum_i u^(k)(i):   G(R7) = (1/5)*sum_{k=1,3,5}(2k+1) U^k.U^k,
G(G2) = (1/4)*sum_{k=1,5}(2k+1) U^k.U^k  (normalisation fixed so the eigenvalues
equal Racah's g(W)=[w1(w1+5)+w2(w2+3)+w3(w3+1)]/10 and g(U)=[u1^2+u1u2+u2^2+5u1+4u2]/12;
checked in tests). Carnall's higher-order operators (three-body T^i, spin-other-orbit +
spin-spin M^k, electrostatically-correlated spin-orbit P^k) are added via `extra=` (see
carnall_ops.py).  Labels: unique conventional SLJ names from a per-J assignment of the
(seniority-resolved) LS terms, see assign_labels().
"""
from __future__ import annotations
from functools import lru_cache
from itertools import combinations
import numpy as np
from sympy.physics.wigner import wigner_3j

L_F = 3
NSO = 14
# spin-orbital index p = 2*(m+3) + (0 for ms=+1/2, 1 for ms=-1/2)
ML = np.array([p // 2 - 3 for p in range(NSO)])
MS2 = np.array([1 if p % 2 == 0 else -1 for p in range(NSO)])  # 2*ms
LETTERS = "SPDFGHIKLMNOQRTUV"


@lru_cache(maxsize=None)
def w3j(a, b, c, d, e, f):
    return float(wigner_3j(a, b, c, d, e, f))


def ck(k, m1, m2):
    m1, m2 = int(m1), int(m2)
    """c^k(l m1, l m2) = (-1)^m1 (2l+1) (l k l;0 0 0)(l k l;-m1 m1-m2 m2)."""
    l = L_F
    return (-1) ** m1 * (2 * l + 1) * w3j(l, k, l, 0, 0, 0) * w3j(l, k, l, -m1, m1 - m2, m2)


def unit_tensor_1e(k, q):
    """14x14 one-electron matrix of u^(k)_q (<l||u^k||l>=1), spin-diagonal."""
    l = L_F
    M = np.zeros((NSO, NSO))
    for p in range(NSO):
        for r in range(NSO):
            if MS2[p] != MS2[r]:
                continue
            m, mp = ML[p], ML[r]
            if m != q + mp:
                continue
            M[p, r] = (-1) ** int(l - m) * w3j(l, k, l, int(-m), q, int(mp))
    return M


def ck_tensor_1e(k, q):
    """one-electron C^(k)_q: <l||C^k||l> = (-1)^l (2l+1) (l k l;000)."""
    return (-1) ** L_F * (2 * L_F + 1) * w3j(L_F, k, L_F, 0, 0, 0) * unit_tensor_1e(k, q)


def angmom_1e():
    """lz, l+, sz, s+ as 14x14 matrices."""
    lz = np.diag(ML.astype(float))
    sz = np.diag(MS2 / 2.0)
    lp = np.zeros((NSO, NSO)); sp = np.zeros((NSO, NSO))
    for r in range(NSO):
        m = ML[r]
        if m < L_F:
            p = r + 2
            lp[p, r] = np.sqrt(L_F * (L_F + 1) - m * (m + 1))
        if MS2[r] == -1:
            sp[r - 1, r] = 1.0
    return lz, lp, sz, sp


class FnBasis:
    """Determinant basis and second-quantised operator builders for 4f^N."""

    def __init__(self, N):
        if not 0 < N < NSO:
            raise ValueError("N must be 1..13")
        self.N = N
        self.dets = [sum(1 << p for p in occ) for occ in combinations(range(NSO), N)]
        self.index = {d: i for i, d in enumerate(self.dets)}
        self.dim = len(self.dets)
        self._cache = {}

    @staticmethod
    def _sign(det, p):
        # number of occupied orbitals below p
        return -1 if bin(det & ((1 << p) - 1)).count("1") % 2 else 1

    def one_body(self, o):
        """Matrix of sum_pq o[p,q] a+_p a_q in the N-electron basis."""
        o = np.asarray(o)
        H = np.zeros((self.dim, self.dim), dtype=o.dtype)
        nz = list(zip(*np.nonzero(o)))
        for j, d in enumerate(self.dets):
            for p, q in nz:
                if not (d >> q) & 1:
                    continue
                s = self._sign(d, q)
                d1 = d ^ (1 << q)
                if (d1 >> p) & 1:
                    continue
                s *= self._sign(d1, p)
                i = self.index[d1 | (1 << p)]
                H[i, j] += s * o[p, q]
        return H

    def coulomb(self):
        """Matrices f_2, f_4, f_6 (coefficients of F^2, F^4, F^6)."""
        if "coul" in self._cache:
            return self._cache["coul"]
        out = {}
        # pair lists grouped by (ML, MS)
        pairs = list(combinations(range(NSO), 2))
        groups = {}
        for (p, q) in pairs:
            groups.setdefault((ML[p] + ML[q], MS2[p] + MS2[q]), []).append((p, q))

        def g(k, a, b, c, d):  # <ab|cd>
            if MS2[a] != MS2[c] or MS2[b] != MS2[d] or ML[a] + ML[b] != ML[c] + ML[d]:
                return 0.0
            return ck(k, ML[a], ML[c]) * ck(k, ML[d], ML[b])

        for k in (2, 4, 6):
            H = np.zeros((self.dim, self.dim))
            V = {}
            for key, plist in groups.items():
                for (p, q) in plist:
                    for (r, s) in plist:
                        v = g(k, p, q, r, s) - g(k, p, q, s, r)
                        if abs(v) > 1e-14:
                            V[(p, q, r, s)] = v
            Vby = {}
            for (p, q, r, s), v in V.items():
                Vby.setdefault((r, s), []).append((p, q, v))
            for j, d in enumerate(self.dets):
                occ = [p for p in range(NSO) if (d >> p) & 1]
                for r, s in combinations(occ, 2):  # r<s ; operator a+p a+q a_s a_r (a_r acts first)
                    sg = self._sign(d, r); d1 = d ^ (1 << r)
                    sg *= self._sign(d1, s); d2 = d1 ^ (1 << s)
                    for p, q, v in Vby.get((r, s), ()):
                        if (d2 >> q) & 1 or (d2 >> p) & 1:
                            continue
                        s2 = sg * self._sign(d2, q); d3 = d2 | (1 << q)
                        s2 *= self._sign(d3, p); d4 = d3 | (1 << p)
                        H[self.index[d4], j] += s2 * v
            out[k] = H
        self._cache["coul"] = out
        return out

    def angular(self):
        """Many-electron Lz, L+, Sz, S+, and L^2, S^2, J^2, Jz, and l.s."""
        if "ang" in self._cache:
            return self._cache["ang"]
        lz, lp, sz, sp = angmom_1e()
        Lz, Lp, Sz, Sp = (self.one_body(x) for x in (lz, lp, sz, sp))
        Lm, Sm = Lp.T, Sp.T
        L2 = Lm @ Lp + Lz @ Lz + Lz
        S2 = Sm @ Sp + Sz @ Sz + Sz
        Jz, Jp = Lz + Sz, Lp + Sp
        J2 = Jp.T @ Jp + Jz @ Jz + Jz
        ls1 = lz @ sz + 0.5 * (lp @ sp.T + lp.T @ sp)
        LS = self.one_body(ls1)
        d = dict(Lz=Lz, Lp=Lp, Sz=Sz, Sp=Sp, L2=L2, S2=S2, J2=J2, Jz=Jz, so=LS)
        self._cache["ang"] = d
        return d

    def U(self, k, q):
        key = ("U", k, q)
        if key not in self._cache:
            self._cache[key] = self.one_body(unit_tensor_1e(k, q))
        return self._cache[key]

    def UdotU(self, k):
        return sum((-1) ** q * self.U(k, q) @ self.U(k, -q) for q in range(-k, k + 1))

    def casimirs(self):
        if "cas" not in self._cache:
            uu = {k: self.UdotU(k) for k in (1, 3, 5)}
            GR7 = 0.2 * (3 * uu[1] + 7 * uu[3] + 11 * uu[5])
            GG2 = 0.25 * (3 * uu[1] + 11 * uu[5])
            self._cache["cas"] = (GG2, GR7)
        return self._cache["cas"]

    def crystal_field(self, bkq):
        """bkq: {(k,q): B} for q>=0 (Wybourne normalisation, complex allowed);
        B^k_{-q} = (-1)^q conj(B^k_q) is added automatically."""
        o = np.zeros((NSO, NSO), dtype=complex)
        for (k, q), B in bkq.items():
            if q < 0:
                continue
            o += B * ck_tensor_1e(k, q)
            if q > 0:
                o += (-1) ** q * np.conj(B) * ck_tensor_1e(k, -q)
        H = self.one_body(o)
        if np.allclose(H.imag, 0):
            H = H.real
        return H

    def ls_projectors(self):
        """{(2S+1, L): projector} onto simultaneous L^2,S^2 eigenspaces."""
        if "lsp" in self._cache:
            return self._cache["lsp"]
        a = self.angular()
        w, v = np.linalg.eigh(a["L2"] + 1000.0 * a["S2"])
        out = {}
        s2 = np.sum(v * (a["S2"] @ v), axis=0)
        l2 = np.sum(v * (a["L2"] @ v), axis=0)
        Sv = np.round((-1 + np.sqrt(1 + 4 * s2)) / 2 * 2) / 2
        Lv = np.round((-1 + np.sqrt(1 + 4 * l2)) / 2).astype(int)
        for key in sorted(set(zip(Sv, Lv))):
            m = (Sv == key[0]) & (Lv == key[1])
            out[(int(round(2 * key[0] + 1)), int(key[1]))] = v[:, m]
        self._cache["lsp"] = out
        return out


    def term_projectors(self):
        """[(mult, L, tau, nslots, P)]: LS eigenspaces split further by the Casimir labels
        (W from G(R7), U from G(G2)) -> one entry per (tau S L) class, ordered by seniority
        within an SL (tau = 1, 2, ...). nslots > 1 if W and U do not separate repeated terms."""
        if "tp" in self._cache:
            return self._cache["tp"]
        GG2, GR7 = self.casimirs()
        out = []
        for (mult, L), P in self.ls_projectors().items():
            M = P.T @ (GR7 + 1e-3 * GG2) @ P
            w, u = np.linalg.eigh((M + M.T) / 2)
            Q = P @ u
            deg = mult * (2 * L + 1)
            i, tau = 0, 0
            while i < len(w):
                j = i
                while j < len(w) and abs(w[j] - w[i]) < 1e-6:
                    j += 1
                tau += 1
                out.append((mult, L, tau, (j - i) // deg, Q[:, i:j]))
                i = j
        self._cache["tp"] = out
        return out


@lru_cache(maxsize=None)
def get_basis(N):
    return FnBasis(N)


def jlabel(J):
    return f"{int(round(2 * J))}/2" if abs(J - round(J)) > 0.25 else f"{int(round(J))}"


def term_label(mult, L):
    return f"{mult}{LETTERS[L]}"


def _qn(x):  # from x(x+1)
    return (-1 + np.sqrt(1 + 4 * max(x, 0))) / 2


class FnResult:
    """Eigen-decomposition of a 4f^N Hamiltonian plus labelled multiplets."""

    def __init__(self, basis, E, V, multiplets, stark=None):
        self.basis, self.E, self.V = basis, E, V
        self.multiplets = multiplets      # list of dicts (free-ion J multiplets)
        self.stark = stark                # list of dicts (when CF given)

    def table(self, emax=None):
        rows = [(m["label"], m["energy"], m["J"], m["purity"]) for m in self.multiplets]
        if emax is not None:
            rows = [r for r in rows if r[1] <= emax]
        return rows


HIGHER = ("T2", "T3", "T4", "T6", "T7", "T8", "M0", "M2", "M4", "P2", "P4", "P6")


def build_hamiltonian(N, F2, F4, F6, zeta, alpha=0.0, beta=0.0, gamma=0.0, bkq=None,
                      extra=None, spin_spin=True):
    """extra: optional dict with any of T2..T8, M0, M2, M4, P2, P4, P6 (Carnall operators,
    see carnall_ops)."""
    b = get_basis(N)
    c = b.coulomb()
    a = b.angular()
    H = F2 * c[2] + F4 * c[4] + F6 * c[6] + zeta * a["so"] + alpha * a["L2"]
    if beta or gamma:
        GG2, GR7 = b.casimirs()
        H = H + beta * GG2 + gamma * GR7
    if extra:
        from . import carnall_ops as CO
        if any(extra.get(k, 0.0) for k in CO.TI):
            T = CO.three_body_operators(b)
            for k in CO.TI:
                if extra.get(k, 0.0):
                    H = H + extra[k] * T[k]
        if any(extra.get(k, 0.0) for k in ("M0", "M2", "M4", "P2", "P4", "P6")):
            MP = CO.magnetic_operators(b, spin_spin=spin_spin)
            for k in ("M0", "M2", "M4", "P2", "P4", "P6"):
                if extra.get(k, 0.0):
                    H = H + extra[k] * MP[k]
    if bkq:
        H = H + b.crystal_field(bkq)
    return H


def _group_multiplets(b, E, V, tol=1e-4):
    """Group degenerate eigenvectors into J multiplets (free-ion case)."""
    a = b.angular()
    groups, i = [], 0
    while i < len(E):
        j = i + 1
        while j < len(E) and abs(E[j] - E[i]) < tol * max(1.0, abs(E[i])) + 1e-6:
            j += 1
        Vg = V[:, i:j]
        J2g = Vg.T @ a["J2"] @ Vg
        wj, uj = np.linalg.eigh((J2g + J2g.T) / 2)
        Jv = np.round(np.array([_qn(x) for x in wj]) * 2) / 2
        for J in sorted(set(Jv)):
            sel = Jv == J
            groups.append((float(np.mean(E[i:j])), float(J), Vg @ uj[:, sel]))
        i = j
    return groups


def _label(b, Vg):
    P = b.ls_projectors()
    w = {key: float(np.sum((Pk.T @ Vg) ** 2)) / Vg.shape[1] for key, Pk in P.items()}
    best = max(w, key=w.get)
    return best, w[best], w


def assign_labels(b, groups):
    """Conventional, unique SLJ labels.

    For each J, the multiplets of that J and the (tau S L) terms that can couple to J are equal
    in number. Going up in energy, each multiplet takes the SL term (summed over seniority) with
    the largest squared component among those not yet used up for this J ("largest LS component, each label
    used once, in order", as in Carnall's tables). E.g. Nd 17241 -> 2G7/2, 19060 -> 4G7/2
    instead of two "4G7/2".
    The name is the bare SLJ (e.g. 2H9/2) for the lowest-energy multiplet carrying that SLJ,
    and 2H(2)9/2 etc. for the higher-lying multiplet of a repeated SL term. Er3+ 28,000 cm-1
    (4G 0.42, 2G1+2G2 0.49) is thus 2G7/2, as in Carnall's table."""
    TP = b.term_projectors()
    out = [None] * len(groups)
    byJ = {}
    for i, (e, J, Vg) in enumerate(groups):
        byJ.setdefault(J, []).append(i)
    for J, idx in byJ.items():
        # SL classes allowed for this J; capacity = number of (tau S L) terms with that SL
        cls = {}
        for (mult, L, tau, ns, P) in TP:
            S = (mult - 1) / 2
            if abs(S - L) - 1e-9 <= J <= S + L + 1e-9:
                c = cls.setdefault((mult, L), dict(cap=0, P=[]))
                c["cap"] += ns
                c["P"].append(P)
        keys = list(cls)
        Wt = np.zeros((len(idx), len(keys)))
        for a, i in enumerate(idx):
            Vg = groups[i][2]
            for c, k in enumerate(keys):
                Wt[a, c] = sum(float(np.sum(np.abs(P.T @ Vg) ** 2)) for P in cls[k]["P"]) / Vg.shape[1]
        cap = np.array([cls[k]["cap"] for k in keys])
        # greedy in energy order: each multiplet takes its largest LS component whose SL still
        # has unused capacity (repeated SL terms -> second occurrence becomes e.g. 2H(2)9/2)
        order = sorted(range(len(idx)), key=lambda a: groups[idx[a]][0])
        for a in order:
            cc = int(np.argmax(np.where(cap > 0, Wt[a], -1.0)))
            cap[cc] -= 1
            mult, L = keys[cc]
            out[idx[a]] = (mult, L, 0, Wt[a, cc])
    return out


def solve(N, F2, F4, F6, zeta, alpha=0.0, beta=0.0, gamma=0.0, bkq=None, shift_ground=True,
          extra=None, spin_spin=True):
    """Diagonalise and label. Returns FnResult with J-multiplets (barycentres when
    a crystal field is given) and, if bkq is given, Stark levels. extra: Carnall higher-order
    parameters T2..T8, M0,M2,M4, P2,P4,P6 (see carnall_ops)."""
    b = get_basis(N)
    H0 = build_hamiltonian(N, F2, F4, F6, zeta, alpha, beta, gamma, extra=extra,
                           spin_spin=spin_spin)
    E0, V0 = np.linalg.eigh(H0)
    groups = _group_multiplets(b, E0, V0)
    labs = assign_labels(b, groups)
    mults = []
    for (e, J, Vg), (mult, L, tau, wt) in zip(groups, labs):
        (m1, L1), pur, comp = _label(b, Vg)
        mults.append(dict(label=term_label(mult, L) + jlabel(J), term=term_label(mult, L),
                          S=(mult - 1) / 2, L=L, J=J, energy=e, purity=wt, tau=tau,
                          largest=term_label(m1, L1), largest_weight=pur, vecs=Vg))
    stark = None
    if bkq:
        H = H0 + b.crystal_field(bkq)
        E, V = np.linalg.eigh(H)
        stark = []
        W = np.array([[np.sum(np.abs(m["vecs"].conj().T @ V[:, n]) ** 2) for m in mults]
                      for n in range(len(E))])
        for n in range(len(E)):
            im = int(np.argmax(W[n]))
            stark.append(dict(energy=float(E[n]), parent=im, weight=float(W[n, im])))
        # barycentres of each parent multiplet from the Stark levels
        for im, m in enumerate(mults):
            w = W[:, im]
            m["energy_free"] = m["energy"]
            m["energy"] = float(np.sum(w * E) / np.sum(w))
        E0, V0 = E, V
    mults_by_id = list(mults)
    eg = min(m["energy"] for m in mults)
    if shift_ground:
        for m in mults:
            m["energy"] -= eg
            if "energy_free" in m:
                m["energy_free"] -= eg
        if stark:
            for s in stark:
                s["energy"] -= eg
    mults.sort(key=lambda m: m["energy"])
    # unique names: bare SLJ for the lowest occurrence, 2H(2)9/2 for a higher repeated SL term
    seen = {}
    for m in mults:
        base = m["label"]
        k = seen.get(base, 0) + 1
        seen[base] = k
        if k > 1:
            m["label"] = f"{m['term']}({k}){jlabel(m['J'])}"
    if stark:
        for st in stark:                      # parent index -> final label
            st["parent"] = mults_by_id[st["parent"]]["label"]
    return FnResult(b, E0, V0, mults, stark)


def solve_params(N, p, bkq=None, casimir=True, higher=True, spin_spin=True):
    """p: dict with F2,F4,F6,zeta,alpha[,beta,gamma][,T2..T8,M0,M2,M4,P2,P4,P6].
    higher=False drops the Carnall three-body/magnetic operators."""
    extra = {k: p[k] for k in HIGHER if k in p} if higher else None
    return solve(N, p["F2"], p["F4"], p["F6"], p["zeta"], p.get("alpha", 0.0),
                 p.get("beta", 0.0) if casimir else 0.0,
                 p.get("gamma", 0.0) if casimir else 0.0, bkq=bkq, extra=extra,
                 spin_spin=spin_spin)


def ls_terms(N):
    """LS terms of 4f^N (no spin-orbit): list of (label, multiplicity count)."""
    b = get_basis(N)
    P = b.ls_projectors()
    out = {}
    for (mult, L), Pk in P.items():
        n = Pk.shape[1] // (mult * (2 * L + 1))
        out[term_label(mult, L)] = n
    return out


def hund_ground(N, zeta_sign_test=True):
    """Ground multiplet label from a calculation with generic Coulomb + SO."""
    r = solve(N, 80000.0, 57000.0, 40000.0, 1000.0 + 140 * N)
    return r.multiplets[0]["label"]
