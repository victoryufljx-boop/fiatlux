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
     with a homogeneous lineshape L (Gaussian of FWHM gamma_h by default; Lorentzian optional,
     whose slow wings over-extend band edges when multiplied by high doping).  Emission weights upper states by their Boltzmann
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


def _sideband_kernel(dn, hw_max, T, n_pts=None):
    """One-phonon vibronic sideband kernels on a grid of spacing dn (cm^-1).
    Phonon density of states g(e) ~ e^2 (Debye-like) below hw_max, zero above.
    Absorption: phonon creation at nu0 + e weighted n(e)+1, annihilation at nu0 - e weighted n(e).
    Emission: the mirror image.  Normalised so the T = 0 sideband has unit weight; detailed
    balance between the two kernels keeps McCumber reciprocity."""
    m = int(np.ceil(hw_max / dn)) + 1
    e = np.arange(1, m + 1) * dn
    g = np.where(e <= hw_max, e ** 2, 0.0)
    g = g / g.sum()
    nb = 1.0 / np.expm1(e / (KB_CM * T)) if T > 0 else np.zeros_like(e)
    k_abs = np.zeros(2 * m + 1); k_em = np.zeros(2 * m + 1)
    k_abs[m + 1:] = g * (nb + 1); k_abs[:m] = (g * nb)[::-1]
    k_em[:m] = (g * (nb + 1))[::-1]; k_em[m + 1:] = g * nb
    return k_abs, k_em


# ------------------------------------------------- Stark-resolved electric-dipole intensities
# Judd-Ofelt with an explicit odd-rank crystal field (Goerller-Walrand & Binnemans, Handbook
# Phys. Chem. Rare Earths 25, ch. 167, eq. for <A|D_rho|B>):
#   <a|D_rho|b> = sum_{lambda,t,p} (-1)^{p+rho} (2 lambda+1) B_{lambda t p}
#                 ( 1 lambda t ; rho, -(p+rho), p ) <a|U^lambda_{p+rho}|b>
# with the intensity parameters from the same site geometry by superposition:
#   B_{lambda t p} = a_{lambda t} sum_L (R0/R_L)^{tau_t} conj(C^t_p(theta_L, phi_L)).
# Summed over complete multiplets this reduces exactly to Judd-Ofelt with the site's
#   Omega_lambda = (2 lambda+1) sum_{t,p} |B_{lambda t p}|^2 / (2t+1),
# and the a_{lambda t} are renormalised so the ensemble-mean Omega equal the measured Omega:
# band integrals are unchanged; what changes is how intensity is shared among Stark components
# (the isotropic form gives every component the same Omega-weighted |U|^2 share).
_W3J = {}


def _w3j(j1, j2, j3, m1, m2, m3):
    k = (j1, j2, j3, m1, m2, m3)
    if k not in _W3J:
        from sympy.physics.wigner import wigner_3j
        _W3J[k] = float(wigner_3j(j1, j2, j3, m1, m2, m3))
    return _W3J[k]


def intensity_B(ligands, ratio=None, tau=None, R0=2.3):
    """{(lam, t, p): B} for one site, intrinsic amplitude 1 for t = lam-1, ratio[lam] for t = lam+1."""
    ratio = ratio or {2: 1.0, 4: 1.0, 6: 1.0}
    out = {}
    for lam in (2, 4, 6):
        for t in (lam - 1, lam + 1):
            a = 1.0 if t == lam - 1 else ratio.get(lam, 1.0)
            ta = (tau or {}).get(t, t + 1.0)
            for p in range(-t, t + 1):
                out[(lam, t, p)] = 0j
            for el, x, y, z in ligands:
                r = np.sqrt(x * x + y * y + z * z)
                th, ph = np.arccos(np.clip(z / r, -1, 1)), np.arctan2(y, x)
                f = a * (R0 / r) ** ta
                for p in range(-t, t + 1):
                    out[(lam, t, p)] += f * np.conj(ckq(t, p, th, ph))
    return out


def site_omegas(B):
    om = {}
    for lam in (2, 4, 6):
        om[lam] = (2 * lam + 1) * sum(abs(v) ** 2 / (2 * t + 1) for (l, t, p), v in B.items() if l == lam)
    return om


def aniso_S(Uw, B, scale):
    """S_ED matrix (cm^2-equivalent, before e^2) between all states from site parameters B, with
    per-lambda scale factors (renormalisation to measured Omega)."""
    M = [None, None, None]
    for ir, rho in enumerate((-1, 0, 1)):
        acc = None
        for lam in (2, 4, 6):
            for q in range(-lam, lam + 1):
                p = q - rho
                c = 0j
                for t in (lam - 1, lam + 1):
                    if abs(p) > t:
                        continue
                    w = _w3j(1, lam, t, rho, -q, p)
                    if w == 0.0:
                        continue
                    c += (-1) ** (p + rho) * (2 * lam + 1) * B[(lam, t, p)] * w
                if c == 0:
                    continue
                term = (c * np.sqrt(scale[lam])) * Uw[lam][q + lam]
                acc = term if acc is None else acc + term
        M[ir] = acc
    return sum(np.abs(m) ** 2 for m in M)


def _apply_level_anchor(m, diags, anchor):
    """See ensemble_spectra(level_anchor=...)."""
    kramers = (m.N % 2 == 1) if hasattr(m, "N") else None
    out = [(E.copy(), W, par) for E, W, par in diags]
    for label, a in anchor.items():
        im = m.mult_index(label)
        rel = []
        for E, W, par in diags:
            e = np.sort(E[par == im]); rel.append(e - e[0])
        rel = np.array(rel)
        n = rel.shape[1]
        kr = kramers if kramers is not None else (n % 2 == 0 and np.allclose(rel[:, 0::2], rel[:, 1::2], atol=1e-3))
        mean_lv = rel.mean(0)[::2] if kr else rel.mean(0)
        meas = list(a.get("levels") or [])
        xs, ys = [0.0], [0.0]
        for i, v in enumerate(meas[1:], start=1):
            if v is not None and i < len(mean_lv) and mean_lv[i] > xs[-1]:
                xs.append(float(mean_lv[i])); ys.append(float(v))
        xs, ys = np.array(xs), np.array(ys)
        slope = ys[-1] / xs[-1] if len(xs) > 1 else 1.0

        def fmap(x):
            return np.where(x <= xs[-1], np.interp(x, xs, ys), ys[-1] + slope * (x - xs[-1]))
        lows = []
        for k, (E, W, par) in enumerate(out):
            sel = par == im; e = E[sel]; e0 = e.min()
            E[sel] = e0 + fmap(e - e0)
            lows.append(e0 - E[par == 0].min())
        if a.get("E0") is not None:
            d = float(a["E0"]) - float(np.mean(lows))
            for E, W, par in out:
                E[par == im] += d
    return out


def ensemble_spectra(ion: str, sites, spm: dict, omegas, n_index: float, pairs=None,
                     T: float = 295.0, gamma_h: float = 20.0, nu=None, emax: float = 27000.0,
                     lineshape: str = "gauss", neph_std: float = 0.0, sideband_S: float = 0.0,
                     hw_max: float = 1000.0, seed: int = 0, intensity=None, energy_scale: float = 1.0,
                     level_anchor=None):
    """level_anchor: {multiplet label: dict(levels=[0, 425, None, 865], E0=None)}: measured Stark
    levels (cm^-1 above the multiplet's lowest level; one entry per Kramers doublet for odd-electron
    ions; None = not measured) and optionally E0, the measured lowest level of the multiplet above
    the lowest ground level.  Each site's levels of that multiplet go through one monotone map,
    piecewise linear through (site-mean model level, measured level) and proportional beyond the
    last measured level, so the ensemble mean matches the measurement while the site-to-site
    scatter and the crystal-field eigenvectors (intensities) are kept.
    intensity: None (isotropic Judd-Ofelt per Stark component) or dict(ratio={2:..,4:..,6:..},
    tau=None, R0=2.3) for Stark-resolved intensities from the odd-rank crystal field (aniso_S).
    neph_std: per-site fractional scatter of all multiplet energies (first-order
    nephelauxetic variation, E -> E (1 + eps), eps ~ N(0, neph_std)); an inhomogeneous
    mechanism the crystal-field sites alone miss.  sideband_S: relative strength of the
    one-phonon vibronic sideband (fraction of each line's intensity moved into it at T = 0)."""
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
    e_lo = [[] for _ in m.multiplets]; e_hi = [[] for _ in m.multiplets]
    site_edges = []                       # per site {multiplet index: (lowest, highest) Stark level}
    xe, xm = chi_ed(n_index), chi_md(n_index)
    hg = gamma_h / 2.0
    e_free = np.array([mm["energy"] for mm in m.multiplets])
    rng_n = np.random.default_rng(10007 + seed)
    if intensity is not None and intensity != "uniform":
        Bs = [intensity_B(lig, intensity.get("ratio"), intensity.get("tau"), intensity.get("R0", 2.3)) for lig in sites]
        oms = [site_omegas(b) for b in Bs]
        iscale = {lam: O[i] / max(np.mean([o[lam] for o in oms]), 1e-300) for i, lam in enumerate((2, 4, 6))}
    diags = [m.stark(spm_bkq(lig, spm)) for lig in sites]
    if level_anchor:
        diags = _apply_level_anchor(m, diags, level_anchor)
    for isite, lig in enumerate(sites):
        E, W, par = diags[isite]
        if neph_std > 0:
            E = E + rng_n.normal(0.0, neph_std) * e_free[par]
        if energy_scale != 1.0:                       # mean shift of all multiplets, E -> E (1 + d)
            E = E + (energy_scale - 1.0) * e_free[par]
        Uw = {lam: [W.conj().T @ Uq @ W for Uq in m.Uq[lam]] for lam in (2, 4, 6)}
        Mw = [W.conj().T @ Mq @ W for Mq in m.MDq]
        if intensity is None or intensity == "uniform":
            S_ed = sum(O[i] * sum(np.abs(X) ** 2 for X in Uw[lam])
                       for i, lam in enumerate((2, 4, 6))) * E_ESU ** 2
            if intensity == "uniform":       # every Stark component of a multiplet pair equally strong;
                S_ed = S_ed.copy()           # the pair's total line strength (hence tau_rad) is unchanged
                for ia_ in range(len(m.multiplets)):
                    a_i = np.where(par == ia_)[0]
                    for ib_ in range(ia_):
                        b_i = np.where(par == ib_)[0]
                        if len(a_i) and len(b_i):
                            blk = S_ed[np.ix_(a_i, b_i)]; S_ed[np.ix_(a_i, b_i)] = blk.mean()
                            S_ed[np.ix_(b_i, a_i)] = blk.mean()
        else:
            S_ed = aniso_S(Uw, Bs[isite], iscale) * E_ESU ** 2
        S_md = sum(np.abs(X) ** 2 for X in Mw) * MU_B_ESU ** 2
        se_site = {}
        site_edges.append(se_site)
        for im in range(len(m.multiplets)):
            e = E[par == im]
            if len(e):
                bary[im] += e.mean(); split[im].append(e.max() - e.min())
                e_lo[im].append(e.min()); e_hi[im].append(e.max())
                se_site[m.multiplets[im]["label"]] = (float(e.min()), float(e.max()))
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
    if lineshape == "lorentz":
        ker = (hg / np.pi) / (ker_x ** 2 + hg * hg)      # per cm^-1; integrates to ~1
    else:   # Gaussian of FWHM gamma_h: far wings of RT 4f lines fall faster than Lorentzian
        sg = gamma_h / 2.3548
        ker = np.exp(-0.5 * (ker_x / sg) ** 2) / (sg * np.sqrt(2 * np.pi))
    if sideband_S > 0:
        k_abs, k_em = _sideband_kernel(dn, hw_max, T)
        for p in pairs:
            sa[p] = (1 - sideband_S) * sa[p] + sideband_S * np.convolve(sa[p], k_abs, mode="same")
            se[p] = (1 - sideband_S) * se[p] + sideband_S * np.convolve(se[p], k_em, mode="same")
    for p in pairs:
        sa[p] = np.convolve(sa[p] / ns, ker * dn, mode="same")
        se[p] = np.convolve(se[p] / ns, ker * dn, mode="same")
    bary /= ns
    return dict(nu=nu, sigma_abs=sa, sigma_em=se, pairs=pairs,
                bary={mm["label"]: float(bary[i]) for i, mm in enumerate(m.multiplets)},
                split={mm["label"]: (float(np.mean(split[i])) if split[i] else None)
                       for i, mm in enumerate(m.multiplets)},
                edges={mm["label"]: (float(np.mean(e_lo[i])), float(np.mean(e_hi[i])))
                       for i, mm in enumerate(m.multiplets) if e_lo[i]},
                site_edges=site_edges)


def peak(nu, s):
    i = int(np.argmax(s))
    return 1e7 / nu[i], float(s[i])


def fwhm_nm(nu, s):
    i = int(np.argmax(s)); h = s[i] / 2
    above = np.where(s >= h)[0]
    return 1e7 / nu[above.min()] - 1e7 / nu[above.max()]
