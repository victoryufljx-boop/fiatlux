"""Per-ion 4f^N models for the glass forward model.

IonModel builds the complete Carnall free-ion Hamiltonian once (determinant basis, with
T^i, M^k, P^k), diagonalises it, and keeps the free-ion eigenvectors below `emax_basis`
as a truncated basis.  Every site's crystal field is then a small matrix in that basis,
so per-site diagonalisation is cheap even for f9/f10, and J-mixing between the kept
multiplets is retained.

Also precomputed in the truncated basis: U^(2,4,6)_q (electric dipole, Judd-Ofelt) and
L+2S (magnetic dipole), so Stark-resolved line strengths come straight from site
eigenvectors.
"""
from __future__ import annotations

import numpy as np

from ..ln4f.fn_ci import get_basis, build_hamiltonian, ck_tensor_1e, solve_params, HIGHER
from ..ln4f.judd_ofelt import md_ops
from ..refs import ln_refs as R

# Dy3+ and Yb3+ are not in ln_refs.FREE_ION; values from Carnall, Goodman, Rajnak, Rana,
# J. Chem. Phys. 90, 3443 (1989), Table 1 (LaF3), as transcribed in the qlanth project
# (github.com/zia-lab/qlanth, data/carnall.py, commit 14563fd; MIT licence).
# M2=0.56 M0, M4=0.31 M0, P4=0.5 P2, P6=0.1 P2 (Carnall Table footnote convention).
_EXTRA_IONS = {
    "Dy3+": dict(N=9, F2=91903.0, F4=64372.0, F6=49386.0, zeta=1913.0, alpha=18.02,
                 beta=-633.4, gamma=1790.0, T2=329.0, T3=36.0, T4=127.0, T6=-314.0,
                 T7=404.0, T8=315.0, M0=3.39, P2=719.0),
    "Yb3+": dict(N=13, F2=0.0, F4=0.0, F6=0.0, zeta=2928.0, alpha=0.0, beta=0.0, gamma=0.0),
}
for _p in _EXTRA_IONS.values():
    if "M0" in _p:
        _p.update(M2=0.56 * _p["M0"], M4=0.31 * _p["M0"])
    if "P2" in _p:
        _p.update(P4=0.5 * _p["P2"], P6=0.1 * _p["P2"])
    _p["src"] = "Carnall et al. 1989 (LaF3) via qlanth data/carnall.py"


def ion_params(ion: str) -> dict:
    if ion in R.FREE_ION:
        return dict(R.FREE_ION[ion])
    return dict(_EXTRA_IONS[ion])


IONS = list(R.FREE_ION) + list(_EXTRA_IONS)
_MODELS: dict = {}


def get_ion(ion: str, emax_basis: float = 32000.0) -> "IonModel":
    key = (ion, emax_basis)
    if key not in _MODELS:
        _MODELS[key] = IonModel(ion, emax_basis)
    return _MODELS[key]


class IonModel:
    def __init__(self, ion: str, emax_basis: float = 32000.0):
        self.ion = ion
        p = ion_params(ion)
        self.N = N = int(p["N"])
        b = get_basis(N)
        extra = {k: p[k] for k in HIGHER if k in p} or None
        H0 = build_hamiltonian(N, p["F2"], p["F4"], p["F6"], p["zeta"], p.get("alpha", 0.0),
                               p.get("beta", 0.0), p.get("gamma", 0.0), extra=extra)
        res = solve_params(N, p)                     # labelled free-ion multiplets
        E, V = np.linalg.eigh(H0)
        E = E - E.min()
        keep = E <= emax_basis
        self.E_free = E[keep]
        self.V = V[:, keep]                          # truncated basis (dim x m)
        # multiplet membership of every kept free-ion eigenvector
        self.multiplets = []
        owner = -np.ones(self.V.shape[1], int)
        for im, m in enumerate(sorted(res.multiplets, key=lambda m: m["energy"])):
            ov = np.sum(np.abs(m["vecs"].conj().T @ self.V) ** 2, axis=0)
            idx = np.where(ov > 0.5)[0]
            if len(idx) == 0:
                continue
            owner[idx] = len(self.multiplets)
            self.multiplets.append(dict(label=m["label"], J=m["J"], energy=float(m["energy"]),
                                        idx=idx))
        self.owner = owner
        Vt = self.V
        proj = lambda O: Vt.conj().T @ O @ Vt
        # crystal-field operators C^k_q (k = 2,4,6), all q, in the truncated basis
        self.Ckq = {(k, q): proj(b.one_body(ck_tensor_1e(k, q)))
                    for k in (2, 4, 6) for q in range(-k, k + 1)}
        self.Uq = {lam: [proj(b.U(lam, q)) for q in range(-lam, lam + 1)] for lam in (2, 4, 6)}
        self.MDq = [proj(O) for O in md_ops(b)]
        self.H0t = np.diag(self.E_free).astype(complex)

    # ---------------------------------------------------------------- per-site
    def site_hamiltonian(self, bkq: dict) -> np.ndarray:
        """bkq: {(k,q): B} q >= 0, Wybourne normalisation (cm^-1)."""
        H = self.H0t.copy()
        for (k, q), B in bkq.items():
            if q < 0:
                continue
            H += B * self.Ckq[(k, q)]
            if q > 0:
                H += (-1) ** q * np.conj(B) * self.Ckq[(k, -q)]
        return 0.5 * (H + H.conj().T)

    def stark(self, bkq: dict):
        """Eigenvalues (cm^-1, ground = 0), eigenvectors and parent multiplet index."""
        E, W = np.linalg.eigh(self.site_hamiltonian(bkq))
        E = E - E.min()
        wt = np.zeros((len(self.multiplets), len(E)))
        for im, m in enumerate(self.multiplets):
            wt[im] = np.sum(np.abs(W[m["idx"], :]) ** 2, axis=0)
        parent = np.argmax(wt, axis=0)
        return E, W, parent

    def mult_index(self, label: str) -> int:
        for i, m in enumerate(self.multiplets):
            if m["label"] == label:
                return i
        raise KeyError(f"{label} not among kept multiplets of {self.ion}: "
                       f"{[m['label'] for m in self.multiplets]}")
