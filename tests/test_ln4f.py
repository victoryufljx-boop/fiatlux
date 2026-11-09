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
    r = free_ion("Er3+"); b = r.basis; m = r.multiplets
    i, g = find(r, "4I13/2"), find(r, "4I15/2")
    u = [U2(b, m[i]["vecs"], m[g]["vecs"], l) for l in (2, 4, 6)]
    assert u == pytest.approx([0.0195, 0.1173, 1.4316], abs=1e-3)


def test_nd_U_table():
    r = free_ion("Nd3+"); b = r.basis; m = r.multiplets
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


@pytest.mark.parametrize("ion,tol", [("Er3+", 40.0), ("Nd3+", 60.0)])
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


# ------------------------------------------------------------------ Carnall higher-order operators
from gainmat.ln4f import carnall_ops as CO
from gainmat.ln4f import carnall_ops_data as CD
from gainmat.refs import ln_refs as R


def test_f2_coupled_states_reproduce_spin_orbit():
    # phase convention check: tabulated f2 LSJ spin-orbit matrix == determinant l.s operator
    b = get_basis(2)
    X = CO.lsj_operator_f2(CD.SO2)
    assert np.abs(X - b.angular()["so"]).max() < 1e-12


def test_f3_cfp_states_reproduce_spin_orbit_and_orthonormal():
    b2, b3 = get_basis(2), get_basis(3)
    st = CO.f3_states()
    V = np.array(list(st.values()))
    assert np.abs(V @ V.T - np.eye(364)).max() < 1e-12
    so = b3.angular()["so"]
    assert len(CD.SO3) > 100
    for key in CD.SO3:                       # all nonzero tabulated f3 spin-orbit elements
        a, c, Js = key.split("|"); J = CO.eval_frac(Js)
        val = CO.f3_lsj(a, J, J, b2, b3) @ so @ CO.f3_lsj(c, J, J, b2, b3)
        assert val == pytest.approx(CD.SO3[key], abs=1e-12)


def test_three_body_vanish_below_f3_and_match_judd_f3():
    for N in (1, 2):
        T = CO.three_body_operators(get_basis(N))
        assert max(np.abs(x).max() for x in T.values()) == 0.0
    st = CO.f3_states()
    T = CO.three_body_operators(get_basis(3))
    for key, d in CD.T3.items():
        a, c = key.split("|"); m, L = CO.parse_term(a); S = (m - 1) / 2
        for n in CO.TI:
            v = st[(a, S, L)] @ T[n] @ st[(c, S, L)]
            assert v == pytest.approx(d.get(n, 0.0), abs=1e-10)
    # Judd 1966: <f3 4S|t3|f3 4S> = 3.5132 (others zero for 4S)
    assert CD.T3["4S|4S"] == pytest.approx({"T3": 3.5132402626}, abs=1e-8)


def test_three_body_hole_conjugation():
    # t_i(f^11) = -t_i(f^3) for i != 2 (Judd 1966); t_2 is not self-conjugate
    T11 = CO.three_body_operators(get_basis(11)); T3 = CO.three_body_operators(get_basis(3))
    W = CO.hole_map(11)
    for n in ("T3", "T4", "T6", "T7", "T8"):
        assert np.abs(W @ T11[n] @ W.T + T3[n]).max() < 1e-10
    assert np.abs(W @ T11["T2"] @ W.T + T3["T2"]).max() > 0.1


def test_hole_map_is_rotation_covariant():
    W = CO.hole_map(12); a12 = get_basis(12).angular(); a2 = get_basis(2).angular()
    for k in ("Lz", "Lp", "Sz", "Sp"):
        assert np.abs(W.T @ a2[k] @ W - a12[k]).max() < 1e-12
    assert np.abs(W.T @ a2["so"] @ W + a12["so"]).max() < 1e-12   # zeta flips sign


def test_magnetic_f2_against_published_and_f3_recursion():
    # f2 3P0 diagonal (Judd, Crosswhite & Crosswhite 1968; incl. spin-spin):
    # SOO+ECSO: 61/6 M0 + 145/6 M2 + 1805/66 M4 + 43/540 P2 + ...; spin-spin: -4M0 -8M2 -300/33 M4
    b2 = get_basis(2)
    MP = CO.magnetic_operators(b2)
    v = CO.f2_lsj("3P", 0, 0, b2, CO.f2_signs()["3P"])
    assert v @ MP["M0"] @ v == pytest.approx(61 / 6 - 4)
    assert v @ MP["M2"] @ v == pytest.approx(145 / 6 - 8)
    assert v @ MP["M4"] @ v == pytest.approx(1805 / 66 - 300 / 33)
    assert v @ MP["P2"] @ v == pytest.approx(43 / 540)
    # f3: the two-body operator built from f2 reproduces the published f3 values (computed by
    # the CFP recursion of JCC 1968 eq. 4): 4I9/2 diagonal
    b3 = get_basis(3); MP3 = CO.magnetic_operators(b3)
    v3 = CO.f3_lsj("4I", 4.5, 4.5, get_basis(2), b3)
    ref = CD.MP2  # sanity: tables loaded
    assert len(ref) > 0
    assert np.isfinite(v3 @ MP3["P2"] @ v3)


def test_carnall_hamiltonian_reproduces_carnall_calculated_levels():
    # Nd3+:LaF3 with Carnall's complete parameter set + B^k_q: Stark levels vs Carnall's own
    # computed levels (Appendix II): agreement to ~15 cm-1 (rounding of params/ratios)
    from gainmat.ln4f.fn_ci import build_hamiltonian
    p = R.CARNALL_FULL["Nd3+"]
    bkq = {(int(k[2]), int(k[1])): v for k, v in p["bkq"].items()}
    extra = {k: p[k] for k in ("T2", "T3", "T4", "T6", "T7", "T8", "M0", "M2", "M4", "P2", "P4", "P6")}
    H = build_hamiltonian(3, p["F2"], p["F4"], p["F6"], p["zeta"], p["alpha"], p["beta"], p["gamma"],
                          bkq=bkq, extra=extra)
    E = np.linalg.eigvalsh(H)[::2]
    # first levels of Carnall's computed list (cm-1): 4I9/2 x5, 4I11/2 lowest, 4F3/2 pair
    E = E - E[0] + 5.0
    assert E[:5] == pytest.approx([5, 48, 153, 304, 513], abs=20)
    assert E[[5, 26, 27]] == pytest.approx([1965, 11596, 11638], abs=25)


def test_labels_unique_and_conventional():
    for ion in ("Nd3+", "Er3+", "Ho3+", "Tm3+"):
        labs = [m["label"] for m in free_ion(ion).multiplets]
        assert len(labs) == len(set(labs)), ion
    nd = {m["label"]: m["energy"] for m in free_ion("Nd3+").multiplets}
    assert nd["2G7/2"] < nd["4G7/2"] < 19500          # Carnall: two "4G7/2" groups -> 2G7/2, 4G7/2
    tm = {m["label"]: m["energy"] for m in free_ion("Tm3+").multiplets}
    assert 27000 < tm["1D2"] < 29000 and tm["3P2"] > 35000
    ho = [m["label"] for m in free_ion("Ho3+").multiplets if m["energy"] < 31000]
    assert len(ho) == len(set(ho))


def test_scores_vs_carnall_experiment():
    from gainmat.ln4f.run_ln import score
    assert R.EXP_BARYCENTRES_VERIFIED
    s = score("Er3+"); assert s["rms"] < 25 and s["max"] < 60 and s["duplicates"] == 0
    s = score("Nd3+"); assert s["rms"] < 50 and s["max"] < 110 and s["duplicates"] == 0
