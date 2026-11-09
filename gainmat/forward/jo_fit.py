"""Judd-Ofelt parameters from measured absorption bands, using the forward model's own
band integrals (so the fit is consistent with how the model turns Omega into cross sections).

The integrated absorption cross section of every band is linear in Omega (electric dipole) plus a
constant (magnetic dipole):  I_b = sum_lambda Omega_lambda K_b,lambda + M_b.  K and M are obtained
from three unit-Omega runs and one zero-Omega run of the site-ensemble spectra; Omega comes from
non-negative least squares on the measured integrals, weighted by their relative uncertainty.

    om, info = fit_omegas("Tm3+", "silica", {("3F4","3H6"): (lam_nm, sigma_cm2), ...})
"""
from __future__ import annotations

import numpy as np


def _integral(nu, s):
    o = np.argsort(nu)
    return float(np.trapezoid(np.asarray(s)[o], np.asarray(nu)[o]))


def band_kernels(ion, family, pairs, n_sites=16, seed=0, overrides=None):
    from .glass import run
    K = {p: np.zeros(3) for p in pairs}; M = {p: 0.0 for p in pairs}
    base = dict(overrides or {}); base.pop("measured_absorption", None)
    r0 = run(ion, family, n_sites=n_sites, pumps=[], seed=seed, overrides=dict(base, omegas=(0.0, 0.0, 0.0)))
    for p in pairs:
        M[p] = _integral(r0["_spec"]["nu"], r0["_spec"]["sigma_abs"][p])
    for i in range(3):
        om = [0.0, 0.0, 0.0]; om[i] = 1.0
        r = run(ion, family, n_sites=n_sites, pumps=[], seed=seed, overrides=dict(base, omegas=tuple(om)))
        for p in pairs:
            K[p][i] = _integral(r["_spec"]["nu"], r["_spec"]["sigma_abs"][p]) - M[p]
    return K, M


def fit_omegas(ion, family, bands: dict, rel_sigma=0.1, n_sites=16, overrides=None):
    """bands: {(upper, lower): (lam_nm, sigma_abs_cm2)} complete measured absorption bands."""
    from scipy.optimize import nnls
    pairs = [tuple(p) for p in bands]
    K, M = band_kernels(ion, family, pairs, n_sites, overrides=overrides)
    A, y, rows = [], [], []
    for p, (lam, sig) in bands.items():
        p = tuple(p)
        lam = np.asarray(lam, float); sig = np.asarray(sig, float); ok = np.isfinite(sig)
        meas = _integral(1e7 / lam[ok], sig[ok])
        w = 1.0 / (rel_sigma * meas)
        A.append(K[p] * w); y.append((meas - M[p]) * w); rows.append((p, meas))
    A = np.array(A); y = np.array(y)
    om, _ = nnls(A, y)
    info = []
    for (p, meas) in rows:
        model = float(K[p] @ om + M[p])
        info.append(dict(band=f"{p[0]}<-{p[1]}", measured_integral=meas, fitted_integral=model,
                         rel_err=(model - meas) / meas))
    rank = int(np.linalg.matrix_rank(A, tol=1e-3 * np.abs(A).max()))
    return tuple(float(x) for x in om), dict(bands=info, n_bands=len(rows), rank=rank,
                                             note="rank < 3: some Omega are unconstrained by these bands (NNLS sets them to 0)")
