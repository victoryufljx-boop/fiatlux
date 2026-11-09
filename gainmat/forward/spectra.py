"""Glass lineshapes from an ensemble of dopant sites.

For each site:
  1. crystal field from the Newman superposition model (SPM):
        B^k_q = sum_L  Bbar_k(R_L) * conj(C^k_q(theta_L, phi_L)),
        Bbar_k(R) = Bbar_k(R0) * (R0 / R)^t_k          (Wybourne normalisation)
  2. diagonalise the free-ion + crystal-field Hamiltonian (IonModel, truncated basis);
  3. Stark-resolved line strengths between every state a (upper multiplet) and b (lower):
        S_ED(a,b) = e^2 sum_lambda Omega_lambda sum_q |<a|U^lambda_q|b>|^2
        S_MD(a,b) = mu_B^2 sum_q |<a|(L+2S)_q|b>|^2
     (summed over a,b in two multiplets these reproduce the Judd-Ofelt multiplet line
     strength exactly, so the integrated cross sections are JO-normalised);
  4. per state pair  A_ab = 64 pi^4 nu^3 / (3 h) [chi_ED S_ED + chi_MD S_MD]
     and the Fuechtbauer-Ladenburg / Einstein relation
        sigma_ab(nu) = A_ab / (8 pi c n^2 nu^2) * L(nu - nu_ab)
     with a homogeneous Lorentzian L.  Emission weights upper states by their Boltzmann
     population inside the upper multiplet, absorption weights lower states likewise;
     this makes sigma_a and sigma_e obey McCumber reciprocity by construction.
The inhomogeneous width is not a parameter: it comes from the spread of site geometries.
"""
from __future__ import annotations

import numpy as np
from scipy.special import sph_harm_y

from .ions import get_ion
from ..ln4f.judd_ofelt import E_ESU, H_CGS, C_CGS, MU_B_ESU, chi_ed, chi_md

KB_CM = 0.6950348  # Boltzmann constant in cm^-1 / K


# ----------------------------------------------------------------------- crystal field
def ckq(k, q, theta, phi):
    """Racah-normalised C^k_q = sqrt(4 pi/(2k+1)) Y_kq."""
    return np.sqrt(4 * np.pi / (2 * k + 1)) * sph_harm_y(k, q, theta, phi)


def spm_bkq(ligands, spm: dict) -> dict:
    """ligands: iterable of (element, x, y, z) relative to the dopant (Angstrom).
    spm: {element: dict(R0=.., B2=.., B4=.., B6=.., t2=.., t4=.., t6=..)}.
    Returns {(k,q): B} for q >= 0."""
    out = {(k, q): 0j for k in (2, 4, 6) for q in range(0, k + 1)}
    for el, x, y, z in ligands:
        p = spm.get(el)
        if p is None:
            continue
        r = np.sqrt(x * x + y * y + z * z)
        th, ph = np.arccos(np.clip(z / r, -1, 1)), np.arctan2(y, x)
        for k in (2, 4, 6):
            bk = p[f"B{k}"] * (p["R0"] / r) ** p[f"t{k}"]
            for q in range(0, k + 1):
                out[(k, q)] += bk * np.conj(ckq(k, q, th, ph))
    return out


# --------------------------------------------------------------- stochastic site model
_POLY = {}


def _reference_polyhedron(cn: int) -> np.ndarray:
    """Unit vectors of a near-ideal coordination polyhedron (Thomson-like points
    from a deterministic repulsion relaxation; octahedron for 6, square antiprism-like for 8)."""
    if cn in _POLY:
        return _POLY[cn]
    rng = np.random.default_rng(1234 + cn)
    v = rng.normal(size=(cn, 3))
    v /= np.linalg.norm(v, axis=1)[:, None]
    for _ in range(3000):                        # Coulomb repulsion on the sphere
        d = v[:, None, :] - v[None, :, :]
        r = np.linalg.norm(d, axis=2) + np.eye(cn)
        f = np.sum(d / r[..., None] ** 3, axis=1)
        v += 0.02 * f
        v /= np.linalg.norm(v, axis=1)[:, None]
    _POLY[cn] = v
    return v


def random_sites(n_sites: int, ligand: str, cn_mean: float, r_mean: float, r_std: float,
                 ang_std_deg: float = 12.0, seed: int = 0):
    """Stochastic glass sites: coordination number drawn around cn_mean (two nearest
    integers, weighted to give the mean), ligands on a jittered reference polyhedron with
    angular noise ang_std_deg and radial noise r_std, random overall orientation."""
    rng = np.random.default_rng(seed)
    lo = int(np.floor(cn_mean)); w_hi = cn_mean - lo
    sites = []
    for _ in range(n_sites):
        cn = lo + (1 if rng.random() < w_hi else 0)
        v = _reference_polyhedron(cn).copy()
        v += rng.normal(scale=np.radians(ang_std_deg), size=v.shape)
        v /= np.linalg.norm(v, axis=1)[:, None]
        q, _ = np.linalg.qr(rng.normal(size=(3, 3)))          # random rotation
        v = v @ q.T
        r = rng.normal(r_mean, r_std, size=cn)
        sites.append([(ligand, *(ri * vi)) for ri, vi in zip(r, v)])
    return sites


# ------------------------------------------------------------------------- spectra
def _bin(arr, nu, x, w):
    """Add weights w at positions x onto the uniform grid nu by linear interpolation
    (arr then holds integrated strength per bin; divided by the bin width later)."""
    dn = nu[1] - nu[0]
    f = (x - nu[0]) / dn
    i = np.floor(f).astype(int); t = f - i
    ok = (i >= 0) & (i < len(nu) - 1)
    np.add.at(arr, i[ok], w[ok] * (1 - t[ok]) / dn)
    np.add.at(arr, i[ok] + 1, w[ok] * t[ok] / dn)

def nu_grid(lmin_nm=400.0, lmax_nm=5000.0, dnu=2.0):
    return np.arange(1e7 / lmax_nm, 1e7 / lmin_nm, dnu)


def ensemble_spectra(ion: str, sites, spm: dict, omegas, n_index: float, pairs=None,
                     T: float = 295.0, gamma_h: float = 20.0, nu=None, emax: float = 27000.0):
    """Site-averaged absorption and emission cross sections (cm^2) on grid nu (cm^-1)
    for multiplet pairs (upper_label, lower_label).  Returns dict with nu, pairs ->
    (sigma_abs[lower->upper], sigma_em[upper->lower]), mean Stark energies per multiplet,
    and per-multiplet total Stark splitting statistics."""
    m = get_ion(ion)
    nu = nu_grid() if nu is None else nu
    O = np.asarray(omegas, float) * 1e-20
    mults = [i for i, mm in enumerate(m.multiplets) if mm["energy"] <= emax]
    if pairs is None:
        pairs = [(m.multiplets[i]["label"], m.multiplets[j]["label"])
                 for i in mults for j in mults if j < i]
    pidx = [(m.mult_index(u), m.mult_index(l)) for u, l in pairs]
    sa = {p: np.zeros_like(nu) for p in pairs}
    se = {p: np.zeros_like(nu) for p in pairs}
    bary = np.zeros(len(m.multiplets)); split = [[] for _ in m.multiplets]
    xe, xm = chi_ed(n_index), chi_md(n_index)
    hg = gamma_h / 2.0
    for lig in sites:
        E, W, par = m.stark(spm_bkq(lig, spm))
        Uw = {lam: [W.conj().T @ Uq @ W for Uq in m.Uq[lam]] for lam in (2, 4, 6)}
        Mw = [W.conj().T @ Mq @ W for Mq in m.MDq]
        S_ed = sum(O[i] * sum(np.abs(X) ** 2 for X in Uw[lam])
                   for i, lam in enumerate((2, 4, 6))) * E_ESU ** 2
        S_md = sum(np.abs(X) ** 2 for X in Mw) * MU_B_ESU ** 2
        for im in range(len(m.multiplets)):
            e = E[par == im]
            if len(e):
                bary[im] += e.mean(); split[im].append(e.max() - e.min())
        for (u, l), (iu, il) in zip(pairs, pidx):
            a_idx = np.where(par == iu)[0]; b_idx = np.where(par == il)[0]
            if len(a_idx) == 0 or len(b_idx) == 0:
                continue
            Ea, Eb = E[a_idx], E[b_idx]
            fa = np.exp(-(Ea - Ea.min()) / (KB_CM * T)); fa /= fa.sum()
            fb = np.exp(-(Eb - Eb.min()) / (KB_CM * T)); fb /= fb.sum()
            dnu = Ea[:, None] - Eb[None, :]                       # transition energies
            ok = dnu > 50
            A = (64 * np.pi ** 4 * dnu ** 3 / (3 * H_CGS)) * (
                xe * S_ed[np.ix_(a_idx, b_idx)] + xm * S_md[np.ix_(a_idx, b_idx)])
            pref = np.where(ok, A / (8 * np.pi * C_CGS * n_index ** 2 * np.maximum(dnu, 1) ** 2), 0)
            # stick spectra on the grid (linear binning); Lorentzian convolution at the end
            w = pref.ravel(); d = dnu.ravel()
            sel = (w > 0) & (d > nu[0]) & (d < nu[-1])
            if np.any(sel):
                _bin(se[(u, l)], nu, d[sel], (fa[:, None] * pref).ravel()[sel])
                _bin(sa[(u, l)], nu, d[sel], (fb[None, :] * pref).ravel()[sel])
    ns = len(sites)
    dn = nu[1] - nu[0]
    ker_x = np.arange(-int(40 * gamma_h / dn) - 1, int(40 * gamma_h / dn) + 2) * dn
    ker = (hg / np.pi) / (ker_x ** 2 + hg * hg)          # per cm^-1; integrates to ~1
    for p in pairs:
        sa[p] = np.convolve(sa[p] / ns, ker * dn, mode="same")
        se[p] = np.convolve(se[p] / ns, ker * dn, mode="same")
    bary /= ns
    return dict(nu=nu, sigma_abs=sa, sigma_em=se, pairs=pairs,
                bary={mm["label"]: float(bary[i]) for i, mm in enumerate(m.multiplets)},
                split={mm["label"]: (float(np.mean(split[i])) if split[i] else None)
                       for i, mm in enumerate(m.multiplets)})


def peak(nu, s):
    i = int(np.argmax(s))
    return 1e7 / nu[i], float(s[i])


def fwhm_nm(nu, s):
    i = int(np.argmax(s)); h = s[i] / 2
    above = np.where(s >= h)[0]
    return 1e7 / nu[above.min()] - 1e7 / nu[above.max()]
