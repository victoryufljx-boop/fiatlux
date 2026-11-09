"""Reference data for the glass-site tier: superposition-model (SPM) intrinsic
parameters, rare-earth first-shell geometry in glass families, and room-T
homogeneous linewidths.

Conventions (same as the other refs modules):
    value    : number (or None when no sourced number could be obtained)
    src      : citation; "via X" means the number was read in X, which quotes
               the primary work
    verified : True  -> the number was read in the cited source (or in a page /
                        paper quoting it, named in `src`) during this build
               False -> NOT checked; None-valued placeholders are always False
Never add numbers without a source; never flip `verified` without reading it.

STATUS (build of 2026-11-08, 35-min web pass):
  * SPM: the Newman & Ng handbook tables (Crystal Field Handbook, CUP 2000,
    ch. 5 "Superposition model", pp. 83-119, doi:10.1017/CBO9780511524295.007)
    and Rep. Prog. Phys. 52, 699 (1989) were NOT accessible (paywalled). No
    Ln-O / Ln-F Bbar_k(R0), t_k table for Er/Nd/Pr/Tm/Ho/Dy could be read.
    Only one numeric SPM set was found (Yb3+-O2-, see below). All other ions
    carry None placeholders so code fails loudly instead of using invented
    numbers. Fill from the handbook Table 5.x before calibration.
  * Qualitative facts that were read: intrinsic parameters are positive for
    negatively charged ligands; t_k (point charge) = k+1 but fitted power
    laws are "generally higher than electrostatic power laws"
    (arXiv:1809.07397, Sec. 2.2, citing Newman & Ng); Bbar_k are defined in
    Wybourne normalisation there:
      B^k_q = sum_L Bbar_k(R0) (-1)^q C^(k)_{-q}(theta_L, phi_L) (R0/R_L)^t_k.
"""
from __future__ import annotations


def _v(value, src, verified=False, unit=None, note=None):
    d = dict(value=value, src=src, verified=verified)
    if unit:
        d["unit"] = unit
    if note:
        d["note"] = note
    return d


_NN2000 = ("D.J. Newman & B. Ng, Crystal Field Handbook (CUP 2000), ch. 5 "
           "'Superposition model', pp. 83-119 -- NOT accessed in this build")


def _missing_spm(ion, ligand):
    """Placeholder SPM entry: no sourced numbers available."""
    return dict(
        ion=ion, ligand=ligand, normalisation="Wybourne",
        Bbar2=_v(None, _NN2000, unit="cm^-1"),
        Bbar4=_v(None, _NN2000, unit="cm^-1"),
        Bbar6=_v(None, _NN2000, unit="cm^-1"),
        t2=_v(None, _NN2000), t4=_v(None, _NN2000), t6=_v(None, _NN2000),
        R0_A=_v(None, _NN2000, unit="Angstrom"),
        note="no numeric SPM set found in accessible sources; fill from "
             "Newman & Ng 2000 Table 5.x / Rep. Prog. Phys. 52, 699 (1989)",
    )


# ---------------------------------------------------------------------------
# 1. SPM intrinsic parameters, keyed (ion, ligand)
# ---------------------------------------------------------------------------
_DONG2019 = ("H.N. Dong & R. Zhang, Rev. Mex. Fis. 65, 1-4 (2019) (Yb3+ in "
             "YAl3(BO3)4, [YbO6]9- cluster); values stated as taken from "
             "H.N. Dong, J. Wang, X. Shuai, S.Y. Wu, Spectrochim. Acta A 70, "
             "7-10 (2008) for Yb3+ in Bi4Ge3O12")

SPM = {
    ("Yb3+", "O"): dict(
        ion="Yb3+", ligand="O",
        normalisation="UNCLEAR -- source writes A_k(R0); sign of A2 negative, "
                      "so likely not plain Wybourne Bbar_k; check before use",
        Bbar2=_v(-522, _DONG2019, True, "cm^-1", "quoted as A2(R0) ~ -522"),
        Bbar4=_v(66.3, _DONG2019, True, "cm^-1", "quoted as A4(R0) ~ 66.3"),
        Bbar6=_v(4.1, _DONG2019, True, "cm^-1", "quoted as A6(R0) ~ 4.1"),
        t2=_v(3.5, _DONG2019, True), t4=_v(6, _DONG2019, True),
        t6=_v(6, _DONG2019, True),
        R0_A=_v(None, _DONG2019, False, "Angstrom",
                "R0 not stated; 'usually the average metal-ligand bond "
                "length'; host Y-O distance R_H ~ 2.302 A"),
        note="EPR-oriented fit (g-factors), not an optical Stark-level SPM "
             "fit; A_k values are of the Stevens-like A_k<r^k>-free form used "
             "by that group. Treat as indicative only.",
    ),
}
for _ion in ("Er3+", "Nd3+", "Pr3+", "Tm3+", "Ho3+", "Dy3+"):
    for _lig in ("O", "F"):
        SPM[(_ion, _lig)] = _missing_spm(_ion, _lig)
SPM[("Yb3+", "F")] = _missing_spm("Yb3+", "F")

# SPM fits that exist in the literature (abstracts read) but whose numbers
# were not accessible -- the places to look first:
SPM_LEADS = [
    "N. Magnani, G. Amoretti, A. Baraldi, R. Capelletti, Eur. Phys. J. B 29, "
    "79-84 (2002): CF + SPM analysis of Er, Dy, Nd in BaY2F8 (F ligands); "
    "SPM parameters said to be similar to RE3+ in LaF3.",
    "A. Baraldi et al., Phys. Rev. B (2005), Er3+ in BaY2F8, CF levels of 10 "
    "manifolds + Newman SPM analysis.",
    "D.J. Newman, Aust. J. Phys. 31, 79 (1978): empirical orbit-lattice / "
    "SPM parameters for F, Cl, O ligands.",
    "L.-H. Xie & Y.Y. Yeung, Appl. Magn. Reson. 44, 917 (2013): SPM fitted to "
    "levels of Yb3+ in alpha-Al2O3.",
    "J.-H. Li et al., Chin. Phys. B 33, 117601 (2024): SPM Bbar_k fitted for "
    "Yb:GdScO3 (values in full text only).",
    "R.-P. Chai et al., Chin. Phys. B 30, 037601 (2021): SPM for Er3+ "
    "(ErFO4) centres in CdF2/CaF2 (full text blocked in this build).",
]

# Crystal-field strength / Stark-splitting cross-checks
CF_CHECKS = {
    "Nv_Nd_YAG": _v(4215, "Y. Fan et al., Acta Phys. Sin. 73, 044207 (2024), "
                    "Nv from Bkq of Burdick et al.", True, "cm^-1"),
    "Nv_Nd_GdScO3": _v(2709, "Y. Fan et al., Acta Phys. Sin. 73, 044207 (2024)",
                       True, "cm^-1"),
    "Nv_Nd_YAP": _v(3406, "Y. Fan et al., Acta Phys. Sin. 73, 044207 (2024)",
                    True, "cm^-1"),
    "Er_4I15_2_total_splitting_glass": _v(
        (335, 400), "S. Zemon et al., J. Appl. Phys. 69, 6799 (1991) "
        "(FLN of Er in fluorozirconate, fluorophosphate, phosphate, silicate, "
        "silica glasses), abstract", True, "cm^-1",
        "adjacent sublevel spacings 10-140 cm^-1"),
    "Er_4I13_2_total_splitting_glass": _v(None, "not found", False, "cm^-1"),
    "Yb_2F7_2_levels_aluminosilicate": _v(
        (0, 400, 760, 1210), "M.A. Melkumov et al., Preprint No. 5, GPI RAS "
        "(2004), arXiv:1502.02885, Table 5 (short version Quantum Electron. "
        "34, 843 (2004))", True, "cm^-1", "total 2F7/2 splitting 1210 cm^-1"),
    "Yb_2F5_2_levels_aluminosilicate": _v(
        (10245, 10917, 10940), "Melkumov et al. 2004, arXiv:1502.02885, "
        "Table 5", True, "cm^-1", "total 2F5/2 splitting 695 cm^-1"),
    "Yb_2F7_2_levels_phosphosilicate": _v(
        (0, 260, 440, 740), "Melkumov et al. 2004, arXiv:1502.02885, Table 5",
        True, "cm^-1", "total 740 cm^-1"),
    "Yb_2F5_2_levels_phosphosilicate": _v(
        (10260, 10520, 10930), "Melkumov et al. 2004, arXiv:1502.02885, "
        "Table 5", True, "cm^-1", "total 670 cm^-1"),
    "Ho_5I8_total_splitting_silica": _v(
        330, "S.D. Jackson & Y. Li, Electron. Lett. 40, 1474 (2004) "
        "doi:10.1049/el:20046979 -- inferred by authors from 2019-2163 nm "
        "tuning range", True, "cm^-1", "approximate"),
    "Er_ZBLAN_4F9_2_to_4I9_2_dE_over_E": _v(
        0.15, "O. Henderson-Sapir, S.D. Jackson, D.J. Ottaway, Opt. Lett. 41, "
        "1676 (2016) (arXiv:1602.06011), Stark-splitting estimate", True,
        None, "Delta E_L / E_ZL ~ 0.15 -> potential tuning ~0.53 um"),
}

# ---------------------------------------------------------------------------
# 2. First-coordination-shell geometry of RE in glass families
# ---------------------------------------------------------------------------
_PHW = ("P.M. Peters & S.N. Houde-Walter, J. Non-Cryst. Solids (1998) "
        "doi:10.1016/S0022-3093(98)00733-9, 'Local structure of Er3+ in "
        "multicomponent glasses' (EXAFS), abstract")
_KARA = ("M. Karabulut, G.K. Marasinghe, E. Metwalli, A.K. Wittenauer, "
         "R.K. Brow, C.H. Booth, D.K. Shuh, Phys. Rev. B 65, 104206 (2002) "
         "(EXAFS, Nd and Er phosphate glasses), abstract")
_MP_VIA = ("M.A. Marcus & A. Polman, J. Non-Cryst. Solids 136, 260 (1991) "
           "(EXAFS Er:SiO2, 0.3 and 1.6 at.%), via arXiv:2101.11067")


def _site(ligand, cn, r_A, r_spread_A, src, verified, note=None):
    d = dict(ligand=ligand, CN=cn, R_mean_A=r_A, R_spread_A=r_spread_A,
             src=src, verified=verified)
    if note:
        d["note"] = note
    return d


SITES = {
    "silica": [
        _site("O", 6.0, 2.28, 0.09, _MP_VIA, True,
              "R_spread is rms relative displacement (EXAFS sigma); same CN as "
              "crystalline Er2O3 (Er-O 2.229 and 2.317 A)"),
        _site("O", (1.2, 3.75), (2.03, 2.05), None,
              "Paul et al. (0.073-0.075 mol% Er2O3: CN 3.75, 2.05 A) and "
              "Haruna et al. (200 ppm Er fiber: CN 1.2, 2.03 A), via "
              "arXiv:2101.11067", True,
              "anomalous low-concentration EXAFS results; likely defect-bound "
              "Er -- do not use as default"),
    ],
    "aluminosilicate": [
        _site("O", 6.3, (2.21, 2.25), None, _PHW, True,
              "CN 6.3 in the aluminosilicate host; Er-O 2.21-2.25 A across "
              "all hosts studied"),
        _site("O", ">=6", (2.32, 2.36), None,
              "Paul et al.; Haruna et al. (Al or Al+Ge co-doped silica), via "
              "arXiv:2101.11067", True),
    ],
    "phosphate": [
        _site("O", (6.4, 9.0), (2.23, 2.29), None, _KARA, True,
              "Er-O shortens 2.29(1)->2.23(1) A and CN drops 9.0(10) "
              "(ultraphosphate, x<0.15) -> 6.4(9) (metaphosphate, x~0.25); "
              "Nd-O 2.40(1)->2.37(1) A"),
        _site("O", (6.5, 7.0), None, None,
              "G. Mountjoy, J. Non-Cryst. Solids 353, 2029 (2007), abstract",
              True, "0.18<=x<=0.28, large RE; slightly lower for small RE; "
                    "7.5-8 for x<=0.15"),
    ],
    "fluoride_multicomponent": [
        _site("F", 10.0, (2.21, 2.25), None, _PHW, True,
              "CN 10.0 in the fluoride host (Er); the 2.21-2.25 A range is "
              "quoted for Er-O(F) across hosts in the abstract -- check full "
              "text for the Er-F value"),
    ],
    "ZBLAN": [
        _site("F", 11, (2.38, 2.95), None,
              "M. Braglia et al., J. Synchrotron Rad. 6, 509 (1999) (Pr in "
              "ZBLAN, XAFS)", True,
              "NOT a glass fit: PrF3 reference shells 3F@2.38, 4F@2.44, "
              "2F@2.59, 2F@2.95 A were held fixed, only Debye-Waller varied"),
    ],
    "germanate": [_site("O", None, None, None, "not found", False)],
    "tellurite": [_site("O", None, None, None, "not found", False)],
    "fluoroindate": [_site("F", None, None, None, "not found", False)],
    "chalcogenide": [_site("S", None, None, None, "not found", False)],
}

# ---------------------------------------------------------------------------
# 3. Homogeneous linewidths
# ---------------------------------------------------------------------------
HOMOG_WIDTH_cm1 = {
    ("Er3+", "aluminosilicate", "4I13/2->4I15/2", 300): _v(
        49, "Spectral gain hole-burning at 1.53 um in erbium-doped fiber "
        "amplifiers, IEEE Photon. Technol. Lett. 2(4) (1990), AT&T Bell Labs "
        "(OSTI 6610918; authors not listed on the page read)", True, "cm^-1",
        "source gives 11.5 nm at 1.531 um, RT; 49 cm^-1 is our conversion"),
    ("Er3+", "silicate_fiber", "4I13/2->4I15/2", 2): _v(
        (8e6, 10e6), "M.U. Staudt et al., arXiv:quant-ph/0603192 (SHB)", True,
        "Hz", "plateau below 2 K; scales ~T^1.4 above; low-T only"),
    ("Er3+", "fluorozirconate", "4I13/2->4I15/2", 1.6): _v(
        106e6, "quoted in Staudt et al., arXiv:quant-ph/0603192", True, "Hz",
        "low-T only"),
    ("Er3+", "silica/fluoride", "4I13/2->4I15/2", 300): _v(
        None, "Zemon et al., J. Appl. Phys. 69, 6799 (1991): homogeneous "
        "broadening 'a reasonable approximation' for 300 K band (qualitative)",
        False, "cm^-1"),
}
