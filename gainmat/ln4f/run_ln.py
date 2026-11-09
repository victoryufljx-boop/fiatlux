"""Driver for the lanthanide tier: run_dopant() -> schema.DopantProps, and validate()."""
from __future__ import annotations
import numpy as np
from ..schema import DopantProps, Level, Transition
from ..refs import ln_refs as R
from .fn_ci import solve_params
from .lsj import solve_lsj
from .judd_ofelt import rates, find
from .nonrad import w_nr, quantum_eff, FAMILIES

_CACHE = {}


SOLVER = "det"   # default free-ion solver: "det" (determinant CI) or "lsj" (complete Carnall, |tau S L J>)


def free_ion(ion, casimir=True, bkq=None, solver=None):
    """solver='det': determinant CI (F^k, zeta, alpha, beta, gamma [, B^k_q]).
    solver='lsj': complete Carnall effective Hamiltonian in the coupled |tau S L J> basis
    (adds T^i, M^k, P^k; free ion only, so bkq forces 'det')."""
    solver = solver or SOLVER
    if bkq:
        solver = "det"
    key = (ion, casimir, None if bkq is None else tuple(sorted(bkq.items())), solver)
    if key not in _CACHE:
        if solver == "lsj":
            p = dict(R.FREE_ION_CARNALL[ion])
            if not casimir:
                p["beta"] = p["gamma"] = 0.0
            _CACHE[key] = solve_lsj(p["N"], p)
        else:
            p = R.FREE_ION[ion]
            _CACHE[key] = solve_params(p["N"], p, bkq=bkq, casimir=casimir)
    return _CACHE[key]


def score_barycentres(ion, solver="lsj", emax=30000.0):
    """Computed vs experimental (Carnall 1989 App. II/VI) barycentres, ground-barycentre
    referenced. Returns dict(rows, rms, max, n, duplicates)."""
    res = free_ion(ion, solver=solver)
    exp = R.EXP_BARY_CARNALL1989[ion]
    calc = {}
    for m in res.multiplets:
        calc.setdefault(m["label"], m["energy"])
    g = next(iter(exp))
    rows = []
    for lab, e in exp.items():
        if e - exp[g] > emax or lab not in calc:
            continue
        c = calc[lab] - calc[g]
        rows.append(dict(label=lab, calc=round(c, 1), exp=round(e - exp[g], 1), diff=round(c - (e - exp[g]), 1)))
    d = np.array([r["diff"] for r in rows[1:]])
    labs = [m["label"] for m in res.multiplets if m["energy"] <= emax + 3000]
    return dict(rows=rows, rms=float(np.sqrt(np.mean(d ** 2))), max=float(np.max(np.abs(d))),
                n=len(d), duplicates=len(labs) - len(set(labs)))


def _host_info(ion, host_key):
    el = ion.rstrip("3+")
    jo = R.JUDD_OFELT.get(f"{el}:{host_key}") or R.JUDD_OFELT.get(host_key)
    h = R.HOSTS.get(host_key, {})
    return jo, h


def run_dopant(ion, host_key, n_index=None, phonon_max_cm1=None, omegas=None, emax=27000.0,
               family=None, nr_model="auto", bkq=None, solver=None):
    """Levels, JO radiative rates, branching, tau_rad, multiphonon W_nr and eta.

    ion: "Er3+"; host_key: "silica"|"ZBLAN"|"YAG"|"LiYF4" or "Er:ZBLAN"-style key.
    omegas: (O2,O4,O6) 1e-20 cm^2 (defaults to refs table). phonon_max_cm1 feeds the
    host-agnostic gap law (if family is None and nr_model == 'auto', the literature
    family law is used when the host family is tabulated, else the phonon-order form).
    """
    res = free_ion(ion, bkq=bkq, solver=solver)
    jo, h = _host_info(ion, host_key)
    if omegas is None:
        if jo is None:
            raise KeyError(f"no Judd-Ofelt parameters for {ion}:{host_key}; pass omegas")
        omegas = (jo["O2"], jo["O4"], jo["O6"])
    n_index = n_index or (jo or {}).get("n") or h.get("n", 1.5)
    phonon_max_cm1 = phonon_max_cm1 or h.get("hw")
    family = family if family is not None else h.get("family")
    if nr_model == "phonon":
        family = None
    red = res.reduced_elements(emax) if hasattr(res, "reduced_elements") else None
    trans, tot, tau = rates(res, omegas, n_index, emax=emax, red=red)
    m = res.multiplets
    idx = [i for i in range(len(m)) if m[i]["energy"] <= emax]
    out = DopantProps()
    out.levels = [Level(label=m[i]["label"], energy_cm1=float(m[i]["energy"]), J=float(m[i]["J"]))
                  for i in idx]
    W, eta = {}, {}
    for k, i in enumerate(idx):
        if k == 0:
            continue
        gap = m[i]["energy"] - m[idx[k - 1]]["energy"]
        W[i] = w_nr(gap, phonon_max_cm1, family if family in FAMILIES else None,
                    model="family" if family in FAMILIES else "phonon")
        eta[m[i]["label"]] = quantum_eff(tot.get(i, 0.0), W[i])
    for t in trans:
        nxt = idx[idx.index(t["i"]) - 1] == t["j"]
        out.transitions.append(Transition(upper=t["upper"], lower=t["lower"],
                                          wavelength_nm=float(t["wavelength_nm"]),
                                          A_ed_s1=float(t["A_ed"]), A_md_s1=float(t["A_md"]),
                                          branching=float(t["branching"]),
                                          W_nr_s1=float(W.get(t["i"], 0.0)) if nxt else 0.0))
    out.tau_rad_s = {m[i]["label"]: float(tau[i]) for i in tau}
    out.quantum_eff = eta
    out.judd_ofelt = dict(O2=omegas[0], O4=omegas[1], O6=omegas[2])
    out.method = dict(tier="ln4f", hamiltonian=("4f^N |tau S L J>: complete Carnall (F^k,zeta,alpha,beta,gamma,T^i,M^k,P^k)"
                                   if hasattr(res, "reduced_elements") else
                                   "4f^N determinant CI: F2,F4,F6,zeta,alpha,beta,gamma"
                                   + (", B^k_q" if bkq else "")), params=R.FREE_ION[ion]["src"],
                      n_index=n_index, phonon_max_cm1=phonon_max_cm1,
                      nr_model=("family:" + family) if family in FAMILIES else "phonon-order",
                      W_nr_s1={m[i]["label"]: float(W[i]) for i in W})
    return out


def level_table(ion, casimir=True, emax=25000.0):
    res = free_ion(ion, casimir=casimir)
    exp = R.EXP_LEVELS.get(ion, {})
    rows, seen = [], set()
    for mm in res.multiplets:
        lab = mm["label"]
        if lab in seen or mm["energy"] > emax + 2000:
            continue
        seen.add(lab)
        e = exp.get(lab)
        rows.append(dict(label=lab, calc=round(mm["energy"], 1), exp=e,
                         diff=None if e is None else round(mm["energy"] - e, 1)))
    d = [r["diff"] for r in rows if r["diff"] is not None and r["exp"] <= emax]
    rms = float(np.sqrt(np.mean(np.square(d)))) if d else None
    return rows, rms


def _tr(p, up, lo):
    for t in p.transitions:
        if t.upper == up and t.lower == lo:
            return t
    raise KeyError((up, lo))


def validate(verbose=True):
    out = {}
    # ---- levels
    for ion in ("Er3+", "Nd3+"):
        rows, rms = level_table(ion)
        rows0, rms0 = level_table(ion, casimir=False)
        out[f"levels_{ion}"] = dict(rows=rows, rms=rms, rms_no_casimir=rms0)
    # ---- Nd:YAG
    nd = run_dopant("Nd3+", "YAG")
    br = {lo: _tr(nd, "4F3/2", lo).branching for lo in ("4I9/2", "4I11/2", "4I13/2", "4I15/2")}
    out["Nd:YAG"] = dict(branching_4F3_2=br, tau_rad_4F3_2_us=nd.tau_rad_s["4F3/2"] * 1e6,
                         eta_4F3_2=nd.quantum_eff["4F3/2"],
                         tau_exp_us=R.LIFETIMES[("Nd:YAG", "4F3/2")]["tau"] * 1e6,
                         branching_exp=R.BRANCHING[("Nd:YAG", "4F3/2")]["values"],
                         lambda_4F3_2_4I11_2_nm=_tr(nd, "4F3/2", "4I11/2").wavelength_nm)
    # ---- Er: silica vs ZBLAN
    for host in ("silica", "ZBLAN"):
        er = run_dopant("Er3+", host)
        t15 = _tr(er, "4I13/2", "4I15/2"); t27 = _tr(er, "4I11/2", "4I13/2")
        A11 = 1 / er.tau_rad_s["4I11/2"]
        W11 = er.method["W_nr_s1"]["4I11/2"]; W13 = er.method["W_nr_s1"]["4I13/2"]
        out[f"Er:{host}"] = dict(
            lambda_1p5_nm=t15.wavelength_nm, tau_rad_4I13_2_ms=er.tau_rad_s["4I13/2"] * 1e3,
            W_nr_4I13_2=W13, eta_4I13_2=er.quantum_eff["4I13/2"],
            tau_eff_4I13_2_ms=1e3 / (1 / er.tau_rad_s["4I13/2"] + W13),
            lambda_2p7_nm=t27.wavelength_nm, A_2p7_s1=t27.A_ed_s1 + t27.A_md_s1,
            beta_2p7=t27.branching, A_tot_4I11_2=A11, W_nr_4I11_2=W11,
            eta_4I11_2=er.quantum_eff["4I11/2"], tau_eff_4I11_2_ms=1e3 / (A11 + W11),
            tau_exp_4I13_2_ms=R.LIFETIMES[(f"Er:{host}", "4I13/2")]["tau"] * 1e3,
            tau_exp_4I11_2_ms=R.LIFETIMES[(f"Er:{host}", "4I11/2")]["tau"] * 1e3,
            phonon_max_cm1=er.method["phonon_max_cm1"], nr_model=er.method["nr_model"])
    if verbose:
        print_report(out)
    return out


def print_report(out):
    for ion in ("Er3+", "Nd3+"):
        d = out[f"levels_{ion}"]
        print(f"\n{ion} multiplets (LaF3 Carnall params)   calc   exp   diff")
        for r in d["rows"]:
            print(f"  {r['label']:8s} {r['calc']:9.0f} {r['exp'] if r['exp'] else '':>7} "
                  f"{'' if r['diff'] is None else round(r['diff']):>6}")
        print(f"  RMS (<25000 cm-1) = {d['rms']:.0f} cm-1 ; without beta/gamma = {d['rms_no_casimir']:.0f}")
    n = out["Nd:YAG"]
    print("\nNd:YAG 4F3/2: branching calc / exp:",
          {k: f"{v:.3f}/{n['branching_exp'][k]}" for k, v in n["branching_4F3_2"].items()})
    print(f"  tau_rad = {n['tau_rad_4F3_2_us']:.0f} us (exp fluorescence {n['tau_exp_us']:.0f} us), "
          f"eta = {n['eta_4F3_2']:.3f}, lambda = {n['lambda_4F3_2_4I11_2_nm']:.0f} nm")
    for host in ("silica", "ZBLAN"):
        e = out[f"Er:{host}"]
        print(f"\nEr:{host} (hw_max={e['phonon_max_cm1']} cm-1, {e['nr_model']})")
        print(f"  4I13/2->4I15/2 {e['lambda_1p5_nm']:.0f} nm: tau_rad={e['tau_rad_4I13_2_ms']:.2f} ms, "
              f"W_nr={e['W_nr_4I13_2']:.2e}/s, tau_eff={e['tau_eff_4I13_2_ms']:.2f} ms "
              f"(exp {e['tau_exp_4I13_2_ms']} ms)")
        print(f"  4I11/2->4I13/2 {e['lambda_2p7_nm']:.0f} nm: A={e['A_2p7_s1']:.1f}/s, beta={e['beta_2p7']:.3f}, "
              f"A_tot={e['A_tot_4I11_2']:.0f}/s, W_nr={e['W_nr_4I11_2']:.2e}/s, eta={e['eta_4I11_2']:.2e}, "
              f"tau_eff={e['tau_eff_4I11_2_ms']:.4f} ms (exp {e['tau_exp_4I11_2_ms']} ms)")


if __name__ == "__main__":
    validate()
