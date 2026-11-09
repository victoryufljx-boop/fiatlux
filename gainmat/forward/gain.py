"""Steady-state rate equations and small-signal gain spectra for a doped glass.

Levels are the 4f^N multiplets (thermalised internally).  Rates:
  radiative  A_ij      Judd-Ofelt (ED + MD) with the site-averaged level energies;
  nonradiative W_i     multiphonon decay to the next lower multiplet (energy-gap law,
                       family parameters or the host-agnostic phonon-order form);
  pump       R         sigma_abs(lambda_p) I / (h nu) from the ground multiplet to the
                       multiplet the pump lands in, plus stimulated emission back.
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


def nonradiative(ion, energies: dict, hw_max: float, family: str | None = None, emax=27000.0):
    """W[i] (s^-1): multiphonon decay of each multiplet to the next lower one."""
    m = get_ion(ion)
    E = np.array([energies[mm["label"]] for mm in m.multiplets])
    order = np.argsort(E)
    W = np.zeros(len(E)); tgt = -np.ones(len(E), int)
    for k in range(1, len(order)):
        i, j = order[k], order[k - 1]
        if E[i] > emax:
            continue
        W[i] = w_nr(E[i] - E[j], hw_max=hw_max, family=family,
                    model="family" if family else "phonon")
        tgt[i] = j
    return W, tgt


def steady_state(A, W, tgt, pump_from, pump_to, R_up, R_down, active):
    """Solve dN/dt = 0 for normalised populations over `active` levels."""
    n = len(A)
    M = np.zeros((n, n))
    for i in active:
        for j in active:
            if A[i, j] > 0:
                M[j, i] += A[i, j]; M[i, i] -= A[i, j]
        if W[i] > 0 and tgt[i] in active:
            M[tgt[i], i] += W[i]; M[i, i] -= W[i]
    M[pump_to, pump_from] += R_up; M[pump_from, pump_from] -= R_up
    M[pump_from, pump_to] += R_down; M[pump_to, pump_to] -= R_down
    act = list(active)
    Ms = M[np.ix_(act, act)].copy()
    Ms[0, :] = 1.0                                   # normalisation replaces one equation
    rhs = np.zeros(len(act)); rhs[0] = 1.0
    x = np.linalg.lstsq(Ms, rhs, rcond=None)[0]
    N = np.zeros(n); N[act] = np.clip(x, 0, None)
    return N / N.sum()


def gain_spectra(ion, spec: dict, omegas, n_index, hw_max, family=None, N_tot_cm3=1e19,
                 pump_nm=980.0, intensities_W_cm2=np.logspace(2, 7, 11), emax=27000.0):
    """Gain g(nu) in dB/m for each pump intensity, using ensemble spectra `spec`
    (from spectra.ensemble_spectra with all pairs among multiplets below emax)."""
    m = get_ion(ion)
    labels = [mm["label"] for mm in m.multiplets]
    energies = {k: (spec["bary"][k] if spec["bary"].get(k) else m.multiplets[i]["energy"])
                for i, k in enumerate(labels)}
    A = multiplet_rates(ion, omegas, n_index, energies, emax)
    W, tgt = nonradiative(ion, energies, hw_max, family, emax)
    active = [i for i, k in enumerate(labels) if energies[k] <= emax]
    nu = spec["nu"]; nup = 1e7 / pump_nm
    g0 = labels.index(min(labels, key=lambda k: energies[k]))
    # which multiplet does the pump land in?  the GSA pair with the largest sigma at nu_p
    best, pt = 0.0, None
    for (u, l) in spec["pairs"]:
        if l == labels[g0]:
            v = np.interp(nup, nu, spec["sigma_abs"][(u, l)])
            if v > best:
                best, pt = v, labels.index(u)
    if pt is None or best <= 0:
        raise ValueError(f"no ground-state absorption at {pump_nm} nm for {ion}")
    s_up = best
    s_dn = float(np.interp(nup, nu, spec["sigma_em"][(labels[pt], labels[g0])]))
    hnu = HC / (pump_nm * 1e-9)
    out = []
    for I in intensities_W_cm2:
        phi = I / hnu                                 # photons / cm^2 / s
        N = steady_state(A, W, tgt, g0, pt, s_up * phi, s_dn * phi, active)
        g = np.zeros_like(nu)
        for (u, l) in spec["pairs"]:
            iu, il = labels.index(u), labels.index(l)
            if N[iu] > 1e-9:
                g += N[iu] * spec["sigma_em"][(u, l)]
            if N[il] > 1e-9:
                g -= N[il] * spec["sigma_abs"][(u, l)]
        out.append(dict(I=I, N={labels[i]: float(N[i]) for i in active},
                        g_dB_m=DB * g * N_tot_cm3 * 100.0))
    tau = {labels[i]: 1.0 / (A[i].sum() + W[i]) for i in active if A[i].sum() + W[i] > 0}
    eta = {labels[i]: A[i].sum() / (A[i].sum() + W[i]) for i in active if A[i].sum() + W[i] > 0}
    return dict(nu=nu, runs=out, pump_level=labels[pt], sigma_pump=s_up, tau=tau, eta=eta,
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
