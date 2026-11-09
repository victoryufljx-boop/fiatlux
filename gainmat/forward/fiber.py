"""Fibre-level model: pump and signal propagating along a doped single-mode fibre.

The local physics is the same as gain.py (multiplet rate equations with Judd-Ofelt radiative
rates, Stark-edge multiphonon decay, optional ion-ion energy transfer, site-ensemble
cross sections).  Here the optical fields are guided:

    I_k(z)    = Gamma_k P_k(z) / A_core                 (mode intensity seen by the dopant)
    dP_k/dz   = +/- [Gamma_k g(lambda_k; N(z)) - alpha_bg] P_k
    g(lambda) = N_tot sum_{u>l} [ n_u sigma_em^{ul}(lambda) - n_l sigma_abs^{ul}(lambda) ]

with the populations n(z) solved locally from every field at that point (pump(s) and signal),
so pump depletion, signal saturation, ground-state absorption and excited-state absorption
are all included.  Overlap Gamma from the Gaussian (Marcuse) fit to the LP01 mode of a
step-index core with a uniformly doped core.  Amplified spontaneous emission is NOT
included (small-signal / low-ASE regime: fine for gain vs pump power at moderate gain,
optimistic at >25-30 dB).

Usage:  amplifier("Pr3+", "ZBLAN", L_m=8, core_um=4, NA=0.13, N_cm3=..., pumps=[(1017, 0.18)],
                  signal_nm=1310, P_signal_W=5e-7)
"""
from __future__ import annotations

import numpy as np

from .gain import _pump_terms, multiplet_rates, nonradiative, steady_state, DB
from .ions import get_ion


def mode_overlap(core_um: float, NA: float, lam_nm: float):
    """(Gamma, w_um, V) for the LP01 mode, Gaussian approximation (Marcuse 1977).  Above
    cutoff (V > 2.405) the fundamental-mode formula is still used: this describes a pump
    launched into the fundamental mode (single-mode fibres) and is approximate otherwise."""
    a = core_um / 2
    V = 2 * np.pi * a * NA / (lam_nm * 1e-3)
    w = a * (0.65 + 1.619 / V ** 1.5 + 2.879 / V ** 6)
    return 1 - np.exp(-2 * a * a / (w * w)), w, V


class LocalModel:
    """Rates and cross sections of one doped glass, ready for repeated local solves."""

    def __init__(self, ion, spec, omegas, n_index, hw_max, family=None, emax=27000.0,
                 energy_transfer=None, quench=None):
        m = get_ion(ion)
        self.labels = [mm["label"] for mm in m.multiplets]
        self.spec = spec
        en = {k: (spec["bary"][k] if spec["bary"].get(k) else m.multiplets[i]["energy"])
              for i, k in enumerate(self.labels)}
        self.A = multiplet_rates(ion, omegas, n_index, en, emax)
        self.W, self.tgt = nonradiative(ion, en, hw_max, family, emax, edges=spec.get("edges"),
                                        site_edges=spec.get("site_edges"))
        if quench:                                  # extra first-order loss of a level (s^-1)
            for lab, w in quench.items():
                i = self.labels.index(lab)
                self.W[i] += w
                if self.tgt[i] < 0:
                    self.tgt[i] = self.labels.index(min(self.labels, key=lambda k: en[k]))
        self.active = [i for i, k in enumerate(self.labels) if en[k] <= emax]
        self.et = energy_transfer or []

    def populations(self, fields, N_tot_cm3, n0=None):
        """fields: list of (lambda_nm, intensity W/cm^2).  Returns fractional populations."""
        terms = []
        for lam, I in fields:
            if I > 0:
                terms += _pump_terms(self.spec, self.labels, lam, I)
        if not self.et:
            return steady_state(self.A, self.W, self.tgt, terms, self.active)
        return solve_with_energy_transfer(self, terms, N_tot_cm3, n0)

    def gain_cm(self, lam_nm, n, N_tot_cm3):
        """Local gain coefficient (1/cm) at lam_nm for populations n."""
        nu, x = self.spec["nu"], 1e7 / lam_nm
        g = 0.0
        for (u, l) in self.spec["pairs"]:
            iu, il = self.labels.index(u), self.labels.index(l)
            se = np.interp(x, nu, self.spec["sigma_em"][(u, l)], left=0, right=0)
            sa = np.interp(x, nu, self.spec["sigma_abs"][(u, l)], left=0, right=0)
            g += n[iu] * se - n[il] * sa
        return g * N_tot_cm3


def solve_with_energy_transfer(lm, terms, N_tot, n0=None, iters=60, tol=1e-9):
    """Steady state with ion-ion energy transfer (quadratic in populations).
    Each process (dA, dA2, aB, aB2, k) moves a donor dA -> dA2 and an acceptor aB -> aB2 at
    rate k N_tot n_dA n_aB per ion (k in cm^3/s).  Solved by successive linearisation: the
    donor-acceptor product is split as k N_tot n_aB (frozen) times n_dA, plus the matching
    acceptor term, which keeps the system linear and converges in a few tens of steps."""
    L = lm.labels
    n = steady_state(lm.A, lm.W, lm.tgt, terms, lm.active) if n0 is None else np.array(n0, float)
    for _ in range(iters):
        ext = list(terms)
        for (dA, dA2, aB, aB2, k) in lm.et:
            i, i2, j, j2 = L.index(dA), L.index(dA2), L.index(aB), L.index(aB2)
            r = k * N_tot * n[i] * n[j]               # events per ion per second
            if r <= 0:
                continue
            # donor: linear in n_i with coefficient k N n_j ; acceptor: linear in n_j with k N n_i
            ext.append((i, i2, k * N_tot * n[j]))
            ext.append((j, j2, k * N_tot * n[i]))
        n_new = steady_state(lm.A, lm.W, lm.tgt, ext, lm.active)
        if np.max(np.abs(n_new - n)) < tol:
            n = n_new
            break
        n = 0.5 * n + 0.5 * n_new
    return n


def amplifier(lm: LocalModel, L_m, core_um, NA, N_cm3, pumps, signal_nm, P_signal_W=1e-6,
              alpha_bg_dB_m=0.0, nz=200, backward_pump=False):
    """Single-pass amplifier.  pumps: list of (lambda_nm, P_W) launched at z = 0 (or z = L when
    backward_pump, solved by shooting on the pump powers).  Returns dict with z, P_s(z),
    P_p(z), net gain dB and internal (no background loss) gain dB."""
    A_core = np.pi * (core_um * 1e-4 / 2) ** 2
    Gs = mode_overlap(core_um, NA, signal_nm)[0]
    Gp = [mode_overlap(core_um, NA, lam)[0] for lam, _ in pumps]
    a_bg = alpha_bg_dB_m / DB / 100.0               # 1/cm
    dz = L_m * 100.0 / nz

    def run(pp0):
        Ps, Pp, n = P_signal_W, list(pp0), None
        zs, PS, PP, G = [0.0], [Ps], [list(Pp)], 0.0
        for k in range(nz):
            fields = [(lam, g * P / A_core) for (lam, _), g, P in zip(pumps, Gp, Pp)]
            fields.append((signal_nm, Gs * Ps / A_core))
            n = lm.populations(fields, N_cm3, n0=n)
            gs = Gs * lm.gain_cm(signal_nm, n, N_cm3)
            Ps *= np.exp((gs - a_bg) * dz); G += gs * dz
            sgn = -1.0 if backward_pump else 1.0
            Pp = [P * np.exp(sgn * (g * lm.gain_cm(lam, n, N_cm3) - a_bg) * dz)
                  for (lam, _), g, P in zip(pumps, Gp, Pp)]
            zs.append((k + 1) * dz / 100); PS.append(Ps); PP.append(list(Pp))
        return np.array(zs), np.array(PS), np.array(PP), G

    if not backward_pump:
        z, PS, PP, G = run([P for _, P in pumps])
    else:                                            # shoot: find P(0) such that P(L) = launched
        guess = [P for _, P in pumps]
        for _ in range(30):
            z, PS, PP, G = run(guess)
            ratio = [P / max(PP[-1][i], 1e-30) for i, (_, P) in enumerate(pumps)]
            if all(abs(r - 1) < 1e-3 for r in ratio):
                break
            guess = [g * r for g, r in zip(guess, ratio)]
    return dict(z_m=z, P_signal_W=PS, P_pump_W=PP,
                gain_dB=float(10 * np.log10(PS[-1] / PS[0])),
                internal_gain_dB=float(DB * G),
                Gamma_signal=float(Gs), Gamma_pump=[float(g) for g in Gp])


def local_model_for(ion, family, n_sites=30, seed=0, overrides=None, energy_transfer=None,
                    quench=None, N_cm3=None):
    """Convenience: build the site-ensemble spectra for (ion, family) via glass.run and wrap
    them in a LocalModel."""
    from .glass import run, FAMILIES, EMAX, omegas_for
    import copy
    r = run(ion, family, n_sites=n_sites, seed=seed, overrides=overrides,
            pumps=[])  # spectra only
    fam = copy.deepcopy(FAMILIES[family]); fam.update(overrides or {})
    om = fam.get("omegas") or omegas_for(ion, family)[0]
    return LocalModel(ion, r["_spec"], om, fam["n"], fam["hw_max"], fam.get("gap_family"),
                      emax=fam.get("emax", EMAX[ion]), energy_transfer=energy_transfer,
                      quench=quench)


def laser(lm: LocalModel, L_m, core_um, NA, N_cm3, pumps, laser_nm, R1=0.99, R2=0.04,
          alpha_bg_dB_m=0.0, nz=120, pump_overlap=None, sweeps=200, tol=1e-4):
    """Steady-state Fabry-Perot fibre laser at a fixed wavelength (Rigrod-type boundary
    conditions, both laser directions and the pump(s) resolved along z, populations solved
    locally from all fields).  pumps: list of (lambda_nm, P_W) launched at z = 0 through the
    R1 mirror.  pump_overlap: override Gamma for cladding pumping (= A_core / A_clad).
    Returns output power through R2 and the threshold estimate.
    Solved by alternating forward / backward sweeps with damping."""
    A_core = np.pi * (core_um * 1e-4 / 2) ** 2
    Gl = mode_overlap(core_um, NA, laser_nm)[0]
    Gp = [pump_overlap if pump_overlap else mode_overlap(core_um, NA, lam)[0] for lam, _ in pumps]
    a_bg = alpha_bg_dB_m / DB / 100.0
    dz = L_m * 100.0 / nz
    Pf = np.full(nz + 1, 1e-3); Pb = np.full(nz + 1, 1e-3)
    pops = [None] * (nz + 1)
    Pp = np.array([[P for _, P in pumps]] * (nz + 1), float)
    out_prev = None
    for it in range(sweeps):
        # local populations and gains with the current fields
        gl = np.zeros(nz + 1); gp = np.zeros((nz + 1, len(pumps)))
        for k in range(nz + 1):
            fields = [(lam, g * Pp[k, i] / A_core) for i, ((lam, _), g) in enumerate(zip(pumps, Gp))]
            fields.append((laser_nm, Gl * (Pf[k] + Pb[k]) / A_core))
            pops[k] = lm.populations(fields, N_cm3, n0=pops[k])
            gl[k] = Gl * lm.gain_cm(laser_nm, pops[k], N_cm3) - a_bg
            gp[k] = [g * lm.gain_cm(lam, pops[k], N_cm3) - a_bg for (lam, _), g in zip(pumps, Gp)]
        # pump: forward propagation from z = 0
        newPp = np.empty_like(Pp); newPp[0] = [P for _, P in pumps]
        for k in range(nz):
            newPp[k + 1] = newPp[k] * np.exp(0.5 * (gp[k] + gp[k + 1]) * dz)
        # laser: round trip; total single-pass log gain
        G = np.sum(0.5 * (gl[1:] + gl[:-1])) * dz
        rt = np.log(R1 * R2) + 2 * G                      # log round-trip gain
        # cumulative gain profile: Pf(z) = Pf(0) exp(int_0^z g), Pb(z) = Pb(L) exp(int_z^L g)
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (gl[1:] + gl[:-1]) * dz)])
        shape_f = np.exp(cum); shape_b = np.exp(G - cum)
        # adjust the circulating power multiplicatively toward rt = 0 (gain clamping)
        scale = np.exp(np.clip(rt, -2, 2) * 0.5)
        Pmax = 5.0 * sum(P for _, P in pumps) + 1e-6       # energy conservation bound
        Pf0 = min(max(Pf[0] * scale, 1e-12), Pmax)
        if not np.isfinite(Pf0):
            Pf0 = 1e-3
        newPf = np.minimum(Pf0 * shape_f, Pmax)
        PbL = newPf[-1] * R2
        newPb = np.minimum(PbL * shape_b, Pmax)
        newPf[0] = R1 * newPb[0] if R1 * newPb[0] > 0 else newPf[0]
        Pf = 0.5 * Pf + 0.5 * newPf; Pb = 0.5 * Pb + 0.5 * newPb; Pp = 0.5 * Pp + 0.5 * newPp
        out = (1 - R2) * Pf[-1]
        if rt < 0 and Pf.max() < 1e-9:
            out = 0.0
            break
        if out_prev is not None and abs(out - out_prev) <= tol * max(out, 1e-6) and abs(rt) < 1e-3:
            break
        out_prev = out
    converged = abs(rt) < 1e-3 or (rt < 0 and out == 0.0)
    method = "sweeps"
    if not converged:
        # Fallback: the damped sweeps can stall (e.g. when the pump field and the laser field
        # compete through ESA).  Solve the gain-clamping condition directly: for a given
        # circulating power Pf(0) relax pump + populations with the laser field fixed, return the
        # log round-trip gain, and bracket its zero (rt falls as the laser saturates the gain).
        ctx = dict(pops=list(pops), Pp=Pp.copy())

        def rt_at(Pf0, n_inner=25):
            pops_l, Pp_l = ctx["pops"], ctx["Pp"].copy()
            Pf_l = np.full(nz + 1, Pf0); Pb_l = np.full(nz + 1, Pf0 * R2)
            for _ in range(n_inner):
                gl_l = np.zeros(nz + 1); gp_l = np.zeros((nz + 1, len(pumps)))
                for k in range(nz + 1):
                    fields = [(lam, g * Pp_l[k, i] / A_core) for i, ((lam, _), g) in enumerate(zip(pumps, Gp))]
                    fields.append((laser_nm, Gl * (Pf_l[k] + Pb_l[k]) / A_core))
                    pops_l[k] = lm.populations(fields, N_cm3, n0=pops_l[k])
                    gl_l[k] = Gl * lm.gain_cm(laser_nm, pops_l[k], N_cm3) - a_bg
                    gp_l[k] = [g * lm.gain_cm(lam, pops_l[k], N_cm3) - a_bg for (lam, _), g in zip(pumps, Gp)]
                newPp = np.empty_like(Pp_l); newPp[0] = [P for _, P in pumps]
                for k in range(nz):
                    newPp[k + 1] = newPp[k] * np.exp(0.5 * (gp_l[k] + gp_l[k + 1]) * dz)
                G_l = np.sum(0.5 * (gl_l[1:] + gl_l[:-1])) * dz
                cum = np.concatenate([[0.0], np.cumsum(0.5 * (gl_l[1:] + gl_l[:-1]) * dz)])
                # Pf(0) fixed; Pb(L) = R2 Pf(L)
                nPf = Pf0 * np.exp(cum); nPb = nPf[-1] * R2 * np.exp(G_l - cum)
                if np.max(np.abs(newPp - Pp_l)) < 1e-6 * max(1.0, Pp_l.max()) and np.allclose(nPf, Pf_l, rtol=1e-4):
                    Pp_l, Pf_l, Pb_l = newPp, nPf, nPb
                    break
                Pp_l = 0.5 * Pp_l + 0.5 * newPp
                Pf_l = np.minimum(0.5 * Pf_l + 0.5 * nPf, Pmax); Pb_l = np.minimum(0.5 * Pb_l + 0.5 * nPb, Pmax)
            ctx["pops"], ctx["Pp"] = pops_l, Pp_l
            return float(np.log(R1 * R2) + 2 * G_l), Pf_l, Pp_l

        lo, hi = 1e-7, Pmax
        r_lo = rt_at(lo)[0]
        if r_lo <= 0:
            out, rt = 0.0, r_lo
        else:
            for _ in range(40):
                mid = np.sqrt(lo * hi)
                r_mid, Pf_m, Pp_m = rt_at(mid)
                if r_mid > 0:
                    lo = mid
                else:
                    hi = mid
                if hi / lo < 1.01:
                    break
            r_mid, Pf_m, Pp_m = rt_at(np.sqrt(lo * hi))
            out, rt, Pp = (1 - R2) * Pf_m[-1], r_mid, Pp_m
        method = "bisection"
    return dict(P_out_W=float(out if rt > -1e-2 else 0.0), log_roundtrip=float(rt),
                P_pump_residual_W=[float(x) for x in Pp[-1]], iterations=it + 1, method=method,
                Gamma_laser=float(Gl), Gamma_pump=[float(g) for g in Gp])


def slope_and_threshold(lm, L_m, core_um, NA, N_cm3, pump_nm, P_list, laser_nm, **kw):
    """Laser output vs launched pump; linear fit above threshold -> (slope, threshold_W, rows)."""
    rows = []
    for P in P_list:
        r = laser(lm, L_m, core_um, NA, N_cm3, [(pump_nm, P)], laser_nm, **kw)
        rows.append((P, r["P_out_W"]))
    x = np.array([p for p, o in rows if o > 1e-4]); y = np.array([o for p, o in rows if o > 1e-4])
    if len(x) < 2:
        return None, None, rows
    s, c = np.polyfit(x, y, 1)
    return float(s), float(-c / s), rows


# ---------------------------------------------------------------------------- co-doped glasses
class CoDopedModel:
    """Several dopant species in one glass, coupled by inter-ion energy transfer.

    members: {name: (LocalModel, N_cm3)}
    et: list of (donor_species, d_from, d_to, acceptor_species, a_from, a_to, k_cm3_s); a donor
        ion in d_from transfers to an acceptor ion in a_from at rate k N_acceptor n_a (per donor
        ion), both ions changing level.  Same-species processes (A == B) are allowed.
    Same interface as LocalModel for amplifier()/laser(): populations(fields) -> {name: n},
    gain_cm(lam, pops) -> sum over species."""

    def __init__(self, members: dict, et=None):
        self.members = members
        self.et = et or []

    def populations(self, fields, N_cm3=None, n0=None, iters=80, tol=1e-9):
        base = {}
        for name, (lm, N) in self.members.items():
            terms = []
            for lam, I in fields:
                if I > 0:
                    terms += _pump_terms(lm.spec, lm.labels, lam, I)
            base[name] = terms
        n = dict(n0) if n0 else {name: steady_state(lm.A, lm.W, lm.tgt, base[name], lm.active)
                                 for name, (lm, N) in self.members.items()}
        for _ in range(iters):
            ext = {k: list(v) for k, v in base.items()}
            for (sa, d0, d1, sb, a0, a1, k) in self.et:
                lA, NA = self.members[sa]; lB, NB = self.members[sb]
                i, i2 = lA.labels.index(d0), lA.labels.index(d1)
                j, j2 = lB.labels.index(a0), lB.labels.index(a1)
                ext[sa].append((i, i2, k * NB * n[sb][j]))          # donor, linear in n_A,i
                ext[sb].append((j, j2, k * NA * n[sa][i]))          # acceptor, linear in n_B,j
            new = {name: steady_state(lm.A, lm.W, lm.tgt, ext[name], lm.active)
                   for name, (lm, N) in self.members.items()}
            dmax = max(np.max(np.abs(new[k] - n[k])) for k in n)
            n = {k: 0.5 * n[k] + 0.5 * new[k] for k in n}
            if dmax < tol:
                break
        return n

    def gain_cm(self, lam_nm, pops, N_cm3=None):
        return sum(lm.gain_cm(lam_nm, pops[name], N) for name, (lm, N) in self.members.items())

    def lifetimes(self, level_species, level, eps=1e-6):
        """Effective small-signal decay rate of one level including energy transfer to the
        other species in their unexcited (ground-state) populations.  Returns seconds."""
        lm, N = self.members[level_species]
        i = lm.labels.index(level)
        rate = lm.A[i].sum() + lm.W[i]
        for (sa, d0, d1, sb, a0, a1, k) in self.et:
            if sa == level_species and d0 == level:
                lB, NB = self.members[sb]
                g = lB.labels.index(min(lB.labels, key=lambda x: lB.spec["bary"].get(x) or 1e9))
                if lB.labels.index(a0) == g:
                    rate += k * NB
        return 1.0 / rate


# ------------------------------------------------------------------- amplifier with ASE
def _xs_matrix(lm, lams_nm):
    """sigma_abs/sigma_em of every multiplet pair at the given wavelengths: (n_pairs, n_lam)."""
    nu, x = lm.spec["nu"], 1e7 / np.asarray(lams_nm, float)
    P = lm.spec["pairs"]
    Sa = np.array([np.interp(x, nu, lm.spec["sigma_abs"][p], left=0, right=0) for p in P])
    Se = np.array([np.interp(x, nu, lm.spec["sigma_em"][p], left=0, right=0) for p in P])
    iu = np.array([lm.labels.index(u) for u, l in P]); il = np.array([lm.labels.index(l) for u, l in P])
    return Sa, Se, iu, il


def amplifier_ase(lm: LocalModel, L_m, core_um, NA, N_cm3, pumps, signals, ase_band_nm=(1450, 1650),
                  ase_step_nm=2.0, alpha_bg_dB_m=0.0, nz=150, n_iter=8, backward_pump=False):
    """Amplifier with forward and backward amplified spontaneous emission.
    pumps: [(lam_nm, P_W)] launched at z = 0 (or at z = L if backward_pump);
    signals: [(lam_nm, P_W)] launched at z = 0.
    ASE in bins of ase_step_nm over ase_band_nm, both directions, spontaneous source
    2 h nu dnu Gamma N_tot sum_u n_u sigma_em(nu) per unit length (two polarisations).
    Solved by forward/backward relaxation sweeps.  Returns gains, ASE spectra, noise figure
    (signal-spontaneous beat limit: NF = (1 + 2 P_ase_fwd(lam_s) / (h nu dnu G)) / G ... , i.e.
    NF = 1/G + 2 n_sp (1 - 1/G) with n_sp from the forward ASE density at the signal)."""
    H, C = 6.62607015e-34, 2.99792458e8
    A_core = np.pi * (core_um * 1e-4 / 2) ** 2
    a_bg = alpha_bg_dB_m / DB / 100.0
    lam_ase = np.arange(ase_band_nm[0], ase_band_nm[1] + 1e-9, ase_step_nm)
    lam_all = np.concatenate([[l for l, _ in pumps], [l for l, _ in signals], lam_ase, lam_ase])
    npmp, nsig, nb = len(pumps), len(signals), len(lam_ase)
    Gam = np.array([mode_overlap(core_um, NA, l)[0] for l in lam_all])
    Sa, Se, iu, il = _xs_matrix(lm, lam_all)
    hnu = H * C / (lam_all * 1e-9)                                   # J
    dnu = C * (ase_step_nm * 1e-9) / (lam_ase * 1e-9) ** 2           # Hz per bin
    dz = L_m * 100.0 / nz
    direction = np.concatenate([[-1.0 if backward_pump else 1.0] * npmp, [1.0] * nsig, [1.0] * nb, [-1.0] * nb])
    P = np.zeros((nz + 1, len(lam_all)))
    P[:, :npmp] = [p for _, p in pumps]; P[:, npmp:npmp + nsig] = [p for _, p in signals]
    n_lev = len(lm.labels)
    pops = np.zeros((nz + 1, n_lev))

    def local(Prow):
        phi = Gam * Prow / A_core / hnu                               # photons/cm^2/s per field
        up = Sa @ phi; dn = Se @ phi
        terms = [(il[k], iu[k], up[k]) for k in range(len(iu))] + [(iu[k], il[k], dn[k]) for k in range(len(iu))]
        n = steady_state(lm.A, lm.W, lm.tgt, terms, lm.active)
        g = Gam * N_cm3 * ((n[iu] @ Se) - (n[il] @ Sa))                  # 1/cm per field
        spont = Gam * N_cm3 * (n[iu] @ Se)                             # 1/cm, emission part
        return n, g, spont

    for it in range(n_iter):
        for k in range(nz + 1):
            pops[k], _, _ = local(P[k])
        # forward-travelling fields
        newP = P.copy()
        fw = direction > 0; bw = ~fw
        start = P[0].copy()
        start[npmp + nsig:npmp + nsig + nb] = 0.0                     # forward ASE starts at 0
        newP[0, fw] = start[fw]
        for k in range(nz):
            _, g, sp = local(np.where(fw, newP[k], P[k]))
            src = np.zeros(len(lam_all)); src[npmp + nsig:npmp + nsig + nb] = 2 * hnu[npmp + nsig:npmp + nsig + nb] * dnu
            newP[k + 1, fw] = (newP[k, fw] * np.exp((g[fw] - a_bg) * dz)
                               + (src * sp)[fw] * dz)
        # backward-travelling fields
        endrow = newP[nz].copy()
        endrow[npmp + nsig + nb:] = 0.0
        if backward_pump:
            endrow[:npmp] = [p for _, p in pumps]
        newP[nz, bw] = endrow[bw]
        for k in range(nz, 0, -1):
            _, g, sp = local(np.where(bw, newP[k], newP[k]))
            src = np.zeros(len(lam_all)); src[npmp + nsig + nb:] = 2 * hnu[npmp + nsig + nb:] * dnu
            newP[k - 1, bw] = (newP[k, bw] * np.exp((g[bw] - a_bg) * dz) + (src * sp)[bw] * dz)
        P = newP
    sig = P[:, npmp:npmp + nsig]
    G = sig[-1] / sig[0]
    ase_f = P[-1, npmp + nsig:npmp + nsig + nb]; ase_b = P[0, npmp + nsig + nb:]
    nf = []
    for s, (lam_s, _) in enumerate(signals):
        pa = np.interp(lam_s, lam_ase, ase_f) / (ase_step_nm)        # W per nm at output
        dnu_s = C * 1e-9 / (lam_s * 1e-9) ** 2                         # Hz per nm
        nsp_term = pa / (H * C / (lam_s * 1e-9) * dnu_s)              # = 2 n_sp (G-1) per 2 pol -> photons
        nf.append(float(10 * np.log10((1 + nsp_term) / G[s])))
    return dict(gain_dB=[float(10 * np.log10(g)) for g in G], noise_figure_dB=nf,
                ase_fwd_W=float(ase_f.sum()), ase_bwd_W=float(ase_b.sum()),
                lam_ase_nm=lam_ase, ase_fwd_spectrum_W=ase_f, ase_bwd_spectrum_W=ase_b,
                pump_out_W=[float(x) for x in P[-1, :npmp]] if not backward_pump else [float(x) for x in P[0, :npmp]])
