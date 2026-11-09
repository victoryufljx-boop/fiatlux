"""Steady-state rate equations and small-signal gain spectra for a doped glass.

Levels are the 4f^N multiplets (thermalised internally).  Rates:
  radiative  A_ij      Judd-Ofelt (ED + MD) with the site-averaged level energies;
  nonradiative W_i     multiphonon decay to the next lower multiplet (energy-gap law,
                       family parameters or the host-agnostic phonon-order form);
  pump       R         sigma(lambda_p) I / (h nu) on every multiplet pair: ground-state and
                       excited-state absorption plus stimulated emission at the pump
                       wavelength(s); an optional second pump gives dual-wavelength pumping.
Steady state  dN/dt = 0,  sum N = N_tot.  Small-signal gain
  g(nu) = sum_{i>j} [ N_i sigma_em^{i->j}(nu) - N_j sigma_abs^{j->i}(nu) ]
covers ground-state absorption, stimulated emission and excited-state absorption of the
signal at once.  The gain window is the envelope over pump intensities.
"""
from __future__ import annotations

import numpy as np

from .ions import get_ion
from ..ln4f.judd_ofelt import E_ESU, H_CGS, MU_B_ESU, chi_ed, chi_md
from ..ln4f.nonrad import w_nr

HC = 6.62607015e-34 * 2.99792458e8   # J m
DB = 10 / np.log(10)                 # 4.343


def multiplet_rates(ion, omegas, n_index, energies: dict, emax=27000.0):
    """Radiative rate matrix A[i, j] (s^-1) between kept multiplets, from Judd-Ofelt
    reduced elements computed in the free-ion basis and the given (site-averaged) energies."""
    m = get_ion(ion)
    idx = [i for i, mm in enumerate(m.multiplets) if energies[mm["label"]] <= emax]
    O = np.asarray(omegas, float) * 1e-20
    xe, xm = chi_ed(n_index), chi_md(n_index)
    A = np.zeros((len(m.multiplets),) * 2)
    for i in idx:
        for j in idx:
            Ei, Ej = energies[m.multiplets[i]["label"]], energies[m.multiplets[j]["label"]]
            if Ei <= Ej:
                continue
            ii, jj = m.multiplets[i]["idx"], m.multiplets[j]["idx"]
            u = [sum(np.sum(np.abs(Uq[np.ix_(ii, jj)]) ** 2) for Uq in m.Uq[lam]) for lam in (2, 4, 6)]
            md = sum(np.sum(np.abs(Mq[np.ix_(ii, jj)]) ** 2) for Mq in m.MDq)
            g = 2 * m.multiplets[i]["J"] + 1
            sig = Ei - Ej
            pref = 64 * np.pi ** 4 * sig ** 3 / (3 * H_CGS * g)
            A[i, j] = pref * (xe * E_ESU ** 2 * float(np.dot(O, u)) + xm * MU_B_ESU ** 2 * md)
    return A


def nonradiative(ion, energies: dict, hw_max: float, family: str | None = None, emax=27000.0,
                 edges: dict | None = None, site_edges: list | None = None, T: float | None = None):
    """W[i] (s^-1): multiphonon decay of each multiplet to the next lower one.  The gap is the
    Stark-edge gap (lowest Stark level of the upper multiplet minus highest of the lower one,
    site-averaged) when `edges` is given, as in the energy-gap-law literature; otherwise the
    barycentre gap.  With `site_edges` (one {label: (lo, hi)} per site) the rate is averaged
    over sites, <W> = <B exp(-alpha dE_site)>: the exponential law weights small-gap sites
    heavily, so this exceeds W(<dE>) by about exp(alpha^2 var(dE) / 2)."""
    m = get_ion(ion)
    mode = "family" if family else "phonon"
    E = np.array([energies[mm["label"]] for mm in m.multiplets])
    order = np.argsort(E)
    W = np.zeros(len(E)); tgt = -np.ones(len(E), int)
    for k in range(1, len(order)):
        i, j = order[k], order[k - 1]
        if E[i] > emax:
            continue
        li, lj = m.multiplets[i]["label"], m.multiplets[j]["label"]
        if site_edges and all(li in s and lj in s for s in site_edges):
            W[i] = float(np.mean([w_nr(max(s[li][0] - s[lj][1], 1.0), hw_max=hw_max, family=family,
                                       model=mode, T=T) for s in site_edges]))
        else:
            if edges and li in edges and lj in edges:
                gap = max(edges[li][0] - edges[lj][1], 1.0)
            else:
                gap = E[i] - E[j]
            W[i] = w_nr(gap, hw_max=hw_max, family=family, model=mode, T=T)
        tgt[i] = j
    return W, tgt


def steady_state(A, W, tgt, pump_terms, active):
    """Solve dN/dt = 0 for normalised populations over `active` levels.
    pump_terms: iterable of (from_level, to_level, rate s^-1) for optical pumping (absorption
    and stimulated emission at every pump wavelength, on every transition)."""
    n = len(A)
    M = np.zeros((n, n))
    for i in active:
        for j in active:
            if A[i, j] > 0:
                M[j, i] += A[i, j]; M[i, i] -= A[i, j]
        if W[i] > 0 and tgt[i] in active:
            M[tgt[i], i] += W[i]; M[i, i] -= W[i]
    for a, b, r in pump_terms:
        if r > 0 and a in active and b in active:
            M[b, a] += r; M[a, a] -= r
    act = list(active)
    Ms = M[np.ix_(act, act)].copy()
    Ms[0, :] = 1.0                                   # normalisation replaces one equation
    rhs = np.zeros(len(act)); rhs[0] = 1.0
    x = np.linalg.lstsq(Ms, rhs, rcond=None)[0]
    N = np.zeros(n); N[act] = np.clip(x, 0, None)
    return N / N.sum()


def _pump_terms(spec, labels, pump_nm, I_W_cm2):
    """(from, to, rate) for a pump of intensity I at pump_nm acting on every multiplet pair:
    absorption l -> u at sigma_abs and stimulated emission u -> l at sigma_em.  This covers
    ground-state absorption, pump excited-state absorption and the 'second pump' of
    dual-wavelength schemes in one rule."""
    nu, nup = spec["nu"], 1e7 / pump_nm
    phi = I_W_cm2 / (HC / (pump_nm * 1e-9))           # photons / cm^2 / s
    out = []
    for (u, l) in spec["pairs"]:
        iu, il = labels.index(u), labels.index(l)
        sa = float(np.interp(nup, nu, spec["sigma_abs"][(u, l)], left=0, right=0))
        se = float(np.interp(nup, nu, spec["sigma_em"][(u, l)], left=0, right=0))
        if sa > 0:
            out.append((il, iu, sa * phi))
        if se > 0:
            out.append((iu, il, se * phi))
    return out


def gain_spectra(ion, spec: dict, omegas, n_index, hw_max, family=None, N_tot_cm3=1e19,
                 pump_nm=980.0, intensities_W_cm2=np.logspace(2, 7, 11), emax=27000.0,
                 pump2_nm=None, intensities2_W_cm2=None, T=None):
    """Gain g(nu) in dB/m for each pump intensity, using ensemble spectra `spec`
    (from spectra.ensemble_spectra with all pairs among multiplets below emax).
    pump2_nm: optional second pump (dual-wavelength pumping); the runs then cover the grid
    intensities_W_cm2 x intensities2_W_cm2."""
    m = get_ion(ion)
    labels = [mm["label"] for mm in m.multiplets]
    energies = {k: (spec["bary"][k] if spec["bary"].get(k) else m.multiplets[i]["energy"])
                for i, k in enumerate(labels)}
    A = multiplet_rates(ion, omegas, n_index, energies, emax)
    W, tgt = nonradiative(ion, energies, hw_max, family, emax, edges=spec.get("edges"),
                          site_edges=spec.get("site_edges"), T=T)
    active = [i for i, k in enumerate(labels) if energies[k] <= emax]
    nu = spec["nu"]; nup = 1e7 / pump_nm
    g0 = labels.index(min(labels, key=lambda k: energies[k]))
    # which multiplet does the (first) pump land in?  the GSA pair with the largest sigma at nu_p
    best, pt = 0.0, None
    for (u, l) in spec["pairs"]:
        if l == labels[g0]:
            v = np.interp(nup, nu, spec["sigma_abs"][(u, l)])
            if v > best:
                best, pt = v, labels.index(u)
    if pt is None or best <= 0:
        raise ValueError(f"no ground-state absorption at {pump_nm} nm for {ion}")
    grid = [(I, 0.0) for I in intensities_W_cm2]
    if pump2_nm is not None:
        I2s = intensities_W_cm2 if intensities2_W_cm2 is None else intensities2_W_cm2
        grid = [(I, J) for I in intensities_W_cm2 for J in I2s]
    out = []
    for I, J in grid:
        terms = _pump_terms(spec, labels, pump_nm, I)
        if pump2_nm is not None and J > 0:
            terms += _pump_terms(spec, labels, pump2_nm, J)
        N = steady_state(A, W, tgt, terms, active)
        g = np.zeros_like(nu)
        for (u, l) in spec["pairs"]:
            iu, il = labels.index(u), labels.index(l)
            if N[iu] > 1e-9:
                g += N[iu] * spec["sigma_em"][(u, l)]
            if N[il] > 1e-9:
                g -= N[il] * spec["sigma_abs"][(u, l)]
        out.append(dict(I=I, I2=J, N={labels[i]: float(N[i]) for i in active},
                        g_dB_m=DB * g * N_tot_cm3 * 100.0))
    tau = {labels[i]: 1.0 / (A[i].sum() + W[i]) for i in active if A[i].sum() + W[i] > 0}
    eta = {labels[i]: A[i].sum() / (A[i].sum() + W[i]) for i in active if A[i].sum() + W[i] > 0}
    return dict(nu=nu, runs=out, pump_level=labels[pt], sigma_pump=best, tau=tau, eta=eta,
                A_tot={labels[i]: float(A[i].sum()) for i in active},
                W={labels[i]: float(W[i]) for i in active})


def envelope(gs: dict):
    return np.max(np.array([r["g_dB_m"] for r in gs["runs"]]), axis=0)


def windows(nu, g_dB_m, loss_dB_m=0.0, threshold_dB_m=0.1, min_width_nm=2.0):
    """Contiguous wavelength ranges (nm) where gain - loss >= threshold."""
    ok = (g_dB_m - loss_dB_m) >= threshold_dB_m
    lam = 1e7 / nu
    out = []
    i = 0
    while i < len(ok):
        if ok[i]:
            j = i
            while j + 1 < len(ok) and ok[j + 1]:
                j += 1
            lo, hi = sorted((lam[i], lam[j]))
            if hi - lo >= min_width_nm:
                k = i + int(np.argmax(g_dB_m[i:j + 1]))
                out.append(dict(lo_nm=round(lo, 1), hi_nm=round(hi, 1),
                                peak_nm=round(lam[k], 1), peak_dB_m=round(float(g_dB_m[k]), 3)))
            i = j + 1
        else:
            i += 1
    return sorted(out, key=lambda w: w["lo_nm"])
