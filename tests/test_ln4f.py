import numpy as np
import pytest

pytestmark = pytest.mark.filterwarnings("ignore::pytest.PytestUnknownMarkWarning")
from gainmat.ln4f.fn_ci import solve, ls_terms, get_basis, term_label
from gainmat.ln4f.judd_ofelt import U2, LS2, find
from gainmat.ln4f.run_ln import free_ion, run_dopant, level_table
from gainmat.ln4f.nonrad import w_nr, w_phonon_order

HUND = {1: "2F5/2", 2: "3H4", 3: "4I9/2", 4: "5I4", 5: "6H5/2", 6: "7F0", 7: "8S7/2",
        8: "7F6", 9: "6H15/2", 10: "5I8", 11: "4I15/2", 12: "3H6", 13: "2F7/2"}
NJ = {1: 2, 2: 13, 3: 41, 4: 107, 5: 198, 6: 295, 7: 327}
NTERMS = {1: 1, 2: 7, 3: 17, 4: 47, 5: 73, 6: 119, 7: 119}   # counted with repetition
FAST = [1, 2, 3, 4, 10, 11, 12, 13]


@pytest.mark.parametrize("N", FAST)
def test_hund_and_counts(N):
    r = solve(N, 80000.0, 57000.0, 40000.0, 1000.0)
    assert r.multiplets[0]["label"] == HUND[N]
    assert len(r.multiplets) == NJ[min(N, 14 - N)]
    assert sum(ls_terms(N).values()) == NTERMS[min(N, 14 - N)]


@pytest.mark.slow  # ~2.5 min total; deselect with -m "not slow"
@pytest.mark.parametrize("N", [5, 6, 7, 8, 9])
def test_hund_half_filled(N):
    r = solve(N, 80000.0, 57000.0, 40000.0, 1000.0)
    assert r.multiplets[0]["label"] == HUND[N]


def test_f2_terms():
    assert set(ls_terms(2)) == {"1S", "1D", "1G", "1I", "3P", "3F", "3H"}


def test_f1_spin_orbit():
    r = solve(1, 0, 0, 0, 1000.0)
    assert r.multiplets[1]["energy"] - r.multiplets[0]["energy"] == pytest.approx(3500.0)
    r = solve(13, 0, 0, 0, 1000.0)          # hole: inverted
    assert r.multiplets[0]["label"] == "2F7/2"


def test_casimir_eigenvalues_f2():
    b = get_basis(2); GG2, GR7 = b.casimirs(); P = b.ls_projectors()
    exp = {"1S": (0, 0), "3P": (1, 1), "3H": (1, 1), "3F": (0.5, 1), "1D": (7 / 6, 1.4)}
    for (m, L), Pk in P.items():
        lab = term_label(m, L)
        if lab in exp:
            assert np.diag(Pk.T @ GG2 @ Pk).mean() == pytest.approx(exp[lab][0], abs=1e-9)
            assert np.diag(Pk.T @ GR7 @ Pk).mean() == pytest.approx(exp[lab][1], abs=1e-9)


def test_er_U_table():
    r = free_ion("Er3+", solver="det"); b = r.basis; m = r.multiplets
    i, g = find(r, "4I13/2"), find(r, "4I15/2")
    u = [U2(b, m[i]["vecs"], m[g]["vecs"], l) for l in (2, 4, 6)]
    assert u == pytest.approx([0.0195, 0.1173, 1.4316], abs=1e-3)


def test_nd_U_table():
    r = free_ion("Nd3+", solver="det"); b = r.basis; m = r.multiplets
    i = find(r, "4F3/2")
    for lo, ref in [("4I9/2", [0, 0.2293, 0.0548]), ("4I11/2", [0, 0.1423, 0.4083])]:
        j = find(r, lo)
        u = [U2(b, m[i]["vecs"], m[j]["vecs"], l) for l in (2, 4, 6)]
        assert u == pytest.approx(ref, abs=5e-3)


def test_md_sum_rule():
    # sum over J' of |<J||L+2S||J'>|^2 within a pure LS term equals (2J+1)<|(L+2S)^2|>; check f1
    r = solve(1, 0, 0, 0, 1000.0); b = r.basis; m = r.multiplets
    tot = sum(LS2(b, m[0]["vecs"], mm["vecs"]) for mm in m)
    # (L+2S)^2 = L^2 + 4 L.S + 4 S^2 ; for 2F5/2: 12 + 4*(-2) + 3 = 7 per state, x6 states
    assert tot == pytest.approx(42.0)


@pytest.mark.parametrize("ion,tol", [("Er3+", 450.0), ("Nd3+", 350.0)])
def test_barycentres(ion, tol):
    rows, rms = level_table(ion)
    assert rms < tol
    first = {r["label"]: r for r in rows}
    assert abs(first["4I13/2" if ion == "Er3+" else "4I11/2"]["diff"]) < 150


def test_ndyag_branching():
    p = run_dopant("Nd3+", "YAG")
    br = {t.lower: t.branching for t in p.transitions if t.upper == "4F3/2"}
    assert br["4I11/2"] > br["4I9/2"] > br["4I13/2"] > br["4I15/2"]
    assert 150e-6 < p.tau_rad_s["4F3/2"] < 450e-6


def test_er_2p7_silica_vs_zblan():
    s = run_dopant("Er3+", "silica"); z = run_dopant("Er3+", "ZBLAN")
    assert s.quantum_eff["4I11/2"] < 0.05 < 0.5 < z.quantum_eff["4I11/2"]
    assert s.quantum_eff["4I13/2"] > 0.99
    assert 5e-3 < s.tau_rad_s["4I13/2"] < 20e-3


def test_phonon_order_monotone():
    assert w_phonon_order(3600, 1100) > w_phonon_order(3600, 580) > w_phonon_order(3600, 350)
    assert w_nr(3600, 1100, "silicate") > 1e3


# ---------------- complete-Carnall |tau S L J> solver (lsj.py) ----------------
from gainmat.ln4f.run_ln import score_barycentres
from gainmat.ln4f.lsj import solve_lsj, duplicate_labels
from gainmat.refs import ln_refs as RR


def test_lsj_counts_and_hund():
    for N, ion in [(3, "Nd3+"), (11, "Er3+")]:
        r = solve_lsj(N, RR.FREE_ION_CARNALL[ion])
        assert len(r.multiplets) == 41
        assert r.multiplets[0]["label"] == HUND[N]
        assert duplicate_labels(r) == 0


def test_lsj_f2_matches_det_solver():
    # same physics (F^k, zeta, alpha, beta, gamma only) must give identical levels in both bases
    p = dict(F2=68878, F4=50347, F6=32901, zeta=751.7, alpha=16.23, beta=-566.6, gamma=1371, M0=0, P2=0)
    a = sorted(m["energy"] for m in solve_lsj(2, p).multiplets)
    b = sorted(m["energy"] for m in solve(2, 68878, 50347, 32901, 751.7, 16.23, -566.6, 1371).multiplets)
    assert np.allclose(a, b, atol=1e-6)


@pytest.mark.parametrize("ion,tol", [("Er3+", 40.0), ("Nd3+", 60.0)])
def test_lsj_barycentres_carnall1989(ion, tol):
    s = score_barycentres(ion, solver="lsj")
    assert s["rms"] < tol and s["duplicates"] == 0


def test_lsj_U_tables():
    r = free_ion("Er3+", solver="lsj"); L = {m["label"]: i for i, m in enumerate(r.multiplets)}
    assert [r.U2(L["4I13/2"], L["4I15/2"], k) for k in (2, 4, 6)] == pytest.approx([0.0195, 0.1173, 1.4316], abs=2e-3)
    r = free_ion("Nd3+", solver="lsj"); L = {m["label"]: i for i, m in enumerate(r.multiplets)}
    assert [r.U2(L["4F3/2"], L["4I11/2"], k) for k in (2, 4, 6)] == pytest.approx([0, 0.1423, 0.4083], abs=3e-3)


def test_lsj_md_f1():
    r = solve_lsj(1, dict(F2=0, F4=0, F6=0, zeta=1000.0, M0=0, P2=0))
    # <2F5/2||L+2S||2F5/2>^2 = J(J+1)(2J+1) g^2, g=6/7
    assert r.MD(0, 0) == pytest.approx(2.5 * 3.5 * 6 * (6 / 7) ** 2)


def test_run_dopant_lsj():
    p = run_dopant("Er3+", "ZBLAN", solver="lsj")
    assert 5e-3 < p.tau_rad_s["4I13/2"] < 20e-3
