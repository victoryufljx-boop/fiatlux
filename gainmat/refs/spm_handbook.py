"""Superposition-model (SPM) intrinsic parameters for trivalent lanthanides,
extracted from D. J. Newman & B. Ng (eds.), *Crystal Field Handbook*
(Cambridge University Press, 2000), ch. 5 "Superposition model" (Newman & Ng),
plus a few related passages (ch. 1, ch. 9).

Page citations are BOOK page numbers ("p. 95"); in the 302-page PDF that was
read, PDF page = book page + 12. Only numbers were extracted; no text or PDF
was copied into the repository.

CONVENTION USED BY OUR CODE (Wybourne)
--------------------------------------
    B^k_q = sum_L Bbar_k(R_L) * conj(C^k_q(theta_L, phi_L)),
    Bbar_k(R) = Bbar_k(R0) * (R0/R)**t_k,        k = 2, 4, 6.

This is exactly the handbook's eq. (5.4)/(5.5) (pp. 86-87): Bbar_k is the
Wybourne B^k_0 of ONE ligand placed on the z axis, and the handbook's
coordination factors g_{k,q}(theta, phi) (Table 5.1, p. 86) are the real
(cos q*phi) combinations of C^k_q, i.e. identical to our C^k_q for q = 0 and
equivalent for real crystal-field parameters. No prefactor is needed for
values the handbook labels Bbar_k.

STEVENS-NORMALISED INTRINSIC PARAMETERS (Abar_k) -> WYBOURNE (Bbar_k)
--------------------------------------------------------------------
Many handbook tables (5.4, 5.5, 5.7, 5.9, 5.14 and several values in the text)
give Abar_k, the single-ligand STEVENS parameter A_k0<r^k> (the <r^k> is
already absorbed; there is NO Stevens theta_k/alpha,beta,gamma factor and no
separate <r^k> to multiply). Eq. (5.4) with Table 5.1 / Table 2.2 (p. 32,
lambda_kq = A_kq<r^k>/B^k_q; lambda_20 = 1/2, lambda_40 = 1/8, lambda_60 = 1/16)
gives, for the axial single-ligand parameters:

    Bbar_2 = 2  * Abar_2
    Bbar_4 = 8  * Abar_4
    Bbar_6 = 16 * Abar_6          (STEVENS_TO_WYBOURNE below)

Checks in the book: p. 97 "B4 = 8A4" and "B2 = 2A2" (Fe2+ garnet); eq. (5.30)
(p. 116) in both normalisations is consistent with these factors; and
Table 5.14 Er:YAG Abar4(R1) = 110.6, Abar6(R1) = 25.2 -> 884.8, 403.2, which is
Table 5.2's YAG entry (885, 403); Er:YGG 95.2, 23.6 -> 761.6, 377.6 = Table
5.2's YGG entry (762, 378); Table 5.7(ii) Cl values x8/x16 reproduce Table 5.6.

WHAT THE HANDBOOK SAYS ABOUT TRENDS / TYPICAL VALUES (quoted as numbers)
-----------------------------------------------------------------------
* Intrinsic parameters of all ranks are POSITIVE for ionic (negatively charged)
  ligands: F-, Cl-, Br-, O2-, S2- (pp. 85, 87, 91, 118).
* With Wybourne normalisation a given ionic ligand has Bbar_2, Bbar_4, Bbar_6
  of similar magnitude (p. 118).
* Er3+ values (Table 5.2) are "typical of trivalent lanthanides, although the
  parameters are somewhat larger for the lighter lanthanides"; Table 5.2 "can
  be used as a first approximation to the intrinsic parameters for all
  trivalent lanthanide ions, especially those with more than half-filled
  shells" (p. 91). -> used here for Tm, Yb (no direct data).
* Trivalent actinides: Bbar_4, Bbar_6 about twice the lanthanide values (p. 94).
* Spectrochemical series (rank-4 strength): I- < Br- < S2- < Cl- < F- < O2-
  (p. 90). O2- and F- expected similar, but O2- values vary strongly between
  hosts (real ionicity) (p. 96). S2- much smaller than O2-, and with a
  different Bbar4/Bbar6 ratio (p. 96).
* F- in cubic fluorites: Abar6/Abar4 = 0.28-0.32 for all Ln3+ and hosts
  (i.e. Bbar6/Bbar4 = 0.56-0.64), except Eu3+:PbF2 (suspect) (pp. 92-93);
  this implies t4 ~= t6 for F-.
* For most negatively charged ligands mu = Bbar6/Bbar4 lies in 0.67-1.0
  (ch. 9, p. 179).
* Power laws: point-charge values t2, t4, t6 = 3, 5, 7, but empirical t_k
  "are not generally in agreement" (p. 87). Empirically t4 > t6 "usually"
  (pp. 94, 112, 118). Specific numbers in the book:
    - ab initio Pr3+-F- (Newman & Curtis): t2 = 3.4(6), t4 = 5.2(4),
      t6 = 5.1(6) (ch. 1, p. 25)                       <- used for F- here
    - Cl-, Table 5.7 (Cs2NaLnCl6 vs LnCl3): t4 = 12-16, t6 = 5-9 (p. 95)
    - Cl-, Pr/Nd:LaCl3 under pressure [TGH93]: Abar4(2.95 A) = 30(4),
      t4 = 6(2), Abar6 = 18(2), t6 = 5.5(2); [GdSH89] Pr: Abar6(2.95) = 16(2),
      t6 = 7.5(1.5), Abar4 = 38(2), t4 = 7(2) (p. 114)
    - O2-, Nd:YAG (two O shells): t4 = 12.2, t6 = 6.2 (p. 112)
    - O2-, zircon structures rank 2: t2 > 5 required, t2 = 7.3 best, "t2 = 7 +- 1
      is more realistic"; Bbar2 = 725(25) cm^-1 for Ln3+ (p. 108)
    - tetravalent actinides, Cl-/Br- (very approximate): t2 = 7, t4 = 11,
      t6 = 8 (p. 96)
  The handbook does NOT give a generic "t4 ~ 5-6, t6 ~ 6 for oxides" rule.
  Caution (p. 106): t_k deduced from one crystal's two ligand shells are
  unreliable; and (p. 93) t_k cannot be extracted from the fluorite series
  because local relaxation around the dopant is unknown (nominal R only).

RANK 2 CAVEATS (pp. 108, 112, 115-119)
-------------------------------------
* Bbar_2 is the least reliable: usually only one B^2_0 per site, near-neighbour
  contributions cancel strongly (tiny angle errors flip the sign), and
  long-range electrostatic contributions break the superposition assumption.
* Table 5.14 garnet Abar_2 values are "unphysical" (negative for gallium
  garnets, too large for aluminium garnets) (p. 112); allowing for local
  distortion, Newman & Edgar [NE76] get Abar2 = 573 (Nd:YAG), 451 (Er:YAG)
  cm^-1 (nearest neighbours), i.e. Bbar2 = 1146, 902.
* Recommended fallback when no empirical Bbar_2 (p. 119): eq. (5.30), a LOWER
  BOUND for the contact (covalency+overlap) part,
      Bbar2_c = (5/11) * (2*Bbar4 - (7/13)*Bbar6)       (bbar2_contact()).
  In LaCl3 the observed Bbar2 exceeds 3/4 Bbar4, i.e. it contains significant
  electrostatic contributions (p. 116).

DATA-QUALITY NOTES
------------------
* Table 5.2 lists R = 2.42 A for LaCl3 (same as LaF3). The real La-Cl
  distance is ~2.95 A and the handbook itself uses R0 = 2.95 A for
  Pr:LaCl3 (p. 114). Cl- values from Tables 5.2/5.6 are therefore given here
  with R0 = 2.95 A (flagged).
* Table 5.2 gives YVO4 R = 2.24 A whereas Table 5.11 gives YVO4 R1 = 2.29 A;
  the Table 5.2 YVO4 Bbar4 = 441 is the lower bound 483 > Bbar4 > 441 (p. 103).
* Table 5.4 R values are nominal host spacings, not the relaxed Ln-F distance.
* Ho2+ (Table 5.4) is DIVALENT and is not used for Ho3+.

SPM_FOR_MODEL RECIPE
--------------------
Direct handbook values where tabulated (src starts "direct"), otherwise
derived as flagged in `src` ("interp"/"adopted"/"scaled"):
  O2- : Nd, Er, Dy from the YAG garnet fits (Table 5.14 x8/x16, nearest shell,
        R0 = R1(YAG) = 2.303 A, Table 5.13); t4, t6 from the two YAG shells via
        eq. (5.6) (Nd quoted p. 112; Er, Dy computed here); Bbar2 from [NE76]
        for Nd, Er (p. 112), t2 = 7 (p. 108). Eu direct from Table 5.8 (YGG,
        2.34 A) - note it is far below the Nd/Dy/Er garnet trend for Bbar4.
        Others linear in Z within the YAG set; Pr = Nd; Tm, Yb = Er (p. 91).
  F-  : Gd, Dy, Er from CaF2 (Table 5.4 x8/x16, R0 = 2.366 A); Eu from CdF2
        (2.333 A); Pr from LaF3 (Yeung & Reid, 2.44 A, p. 93); Nd, Sm, Tb, Ho
        linear in Z after moving neighbours to R0 = 2.366 A with the power law;
        Tm, Yb = Er. t_k = ab initio Pr-F values (p. 25). Bbar2: only Er:LaF3
        480(60) at 2.42 A (Table 5.2); other ions scaled by the Cl- Bbar2 ion
        trend of Table 5.6 (assumption, flagged).
  Cl- : Table 5.6 (LaCl3, Bbar_k direct, R0 = 2.95 A flagged) for Pr, Nd, Eu,
        Tb, Ho, Er; Sm from Table 5.7(ii) x8/x16 (Bbar2 interp Pm/Eu); Gd, Dy
        interp; Tm, Yb = Er. t4, t6 per ion from Table 5.7, interpolated where
        missing; t2 = 7 (no Ln3+ Cl value in the book; see above).
  S2- : only Er:ZnS (Table 5.2, 2.34 A) and Eu in oxysulfides (Table 5.8).
        Er:ZnS adopted for all ions except Eu. NO t_k for S2- in the book:
        t_k copied from the O2- Er entry (pure placeholder, flagged).
"""
from __future__ import annotations

STEVENS_TO_WYBOURNE = {2: 2.0, 4: 8.0, 6: 16.0}   # Bbar_k = factor * Abar_k
PDF_PAGE_OFFSET = 12                              # PDF page = book page + 12
REF = ("D. J. Newman & B. Ng (eds.), Crystal Field Handbook, CUP 2000, "
       "ch. 5 (Newman & Ng), doi:10.1017/CBO9780511524295")

_W, _S = "Wybourne Bbar_k", "Stevens Abar_k (x2/x8/x16 applied)"
_Z = dict(La=57, Ce=58, Pr=59, Nd=60, Pm=61, Sm=62, Eu=63, Gd=64, Tb=65,
          Dy=66, Ho=67, Er=68, Tm=69, Yb=70)
IONS = ["Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb"]


def bbar2_contact(B4, B6):
    """Eq. (5.30), p. 116: lower-bound estimate of the contact part of Bbar2."""
    return 5.0 / 11.0 * (2.0 * B4 - 7.0 / 13.0 * B6)


def _st(a, k):
    return None if a is None else round(a * STEVENS_TO_WYBOURNE[k], 1)


def _row(ion, ligand, host, B2, B4, B6, R0, t2=None, t4=None, t6=None,
         norm=_W, page="", notes=""):
    return dict(ion=ion, ligand=ligand, host=host, Bbar2=B2, Bbar4=B4,
                Bbar6=B6, R0=R0, t2=t2, t4=t4, t6=t6,
                normalisation_original=norm, page=page, notes=notes)


# ---------------------------------------------------------------------------
# 1. Everything tabulated for Ln3+ with O2-, F-, Cl-, S2- (all in Wybourne
#    Bbar_k, cm^-1; Stevens originals already multiplied by 2/8/16).
# ---------------------------------------------------------------------------
SPM_TABLE = [
    # ---- Table 5.2 (p. 91): Er3+, Wybourne -------------------------------
    _row("Er", "O", "YVO4", 725.0, 441.0, 378.0, 2.24, t2=7.0,
         page="p. 91, Table 5.2; Bbar2 & t2 p. 108",
         notes="Bbar2 = 725(25) zircon-structure value, t2 = 7+-1 (p. 108); "
               "Bbar4 = lower bound of 441-483 (p. 103); Table 5.11 gives "
               "R1 = 2.29 A for YVO4"),
    _row("Er", "O", "YGG", None, 762.0, 378.0, 2.34,
         page="p. 91, Table 5.2 (= Table 5.14 x8/x16)",
         notes="nearest-neighbour shell of garnet fit (sec. 5.4.5)"),
    _row("Er", "O", "YAG", None, 885.0, 403.0, 2.30,
         page="p. 91, Table 5.2 (= Table 5.14 x8/x16)",
         notes="nearest-neighbour shell of garnet fit"),
    _row("Er", "S", "ZnS", 400.0, 376.0, 432.0, 2.34, page="p. 91, Table 5.2",
         notes="uncertainties B2 200, B4 32, B6 96 [NN88]"),
    _row("Er", "O", "CaWO4", 820.0, 400.0, 320.0, 2.48, page="p. 91, Table 5.2",
         notes="uncertainties B2 80, B4 24, B6 96 [NN88]"),
    _row("Er", "F", "LaF3", 480.0, 592.0, 310.0, 2.42, page="p. 91, Table 5.2",
         notes="uncertainties B2 60, B4 40, B6 30 [NN88]"),
    _row("Er", "Cl", "LaCl3", 438.0, 264.0, 162.0, 2.42,
         page="p. 91, Table 5.2",
         notes="R = 2.42 as printed is suspect (La-Cl ~2.95 A, cf. p. 114)"),
]

# ---- Table 5.4 (p. 93): F- in cubic fluorites, Stevens Abar4/Abar6 --------
_T54_R = dict(CaF2=2.366, SrF2=2.512, BaF2=2.685, CdF2=2.333, PbF2=2.567)
_T54 = {  # ion: {host: (Abar4, Abar6)}
    "Eu": dict(SrF2=(79.6, 23.9), CdF2=(90.6, 28.1), PbF2=(79.9, 17.1)),
    "Dy": dict(CaF2=(87.8, 25.8), SrF2=(81.5, 23.0), BaF2=(76.5, 20.7),
               CdF2=(90.3, 26.6)),
    "Er": dict(CaF2=(76.6, 22.9), SrF2=(70.5, 19.9), BaF2=(64.4, 17.8),
               CdF2=(80.6, 24.3), PbF2=(66.9, 19.0)),
    "Gd": dict(CaF2=(86.8, 27.9), SrF2=(77.0, 23.3), BaF2=(69.7, 19.2)),
    "Ho2+": dict(CaF2=(83.0, 18.7), SrF2=(70.2, 15.5), BaF2=(59.1, 13.3)),
}
for _ion, _d in _T54.items():
    for _h, (_a4, _a6) in _d.items():
        _n = ("cubic site; R = nominal host spacing (not relaxed); "
              "Abar4=%.1f Abar6=%.1f" % (_a4, _a6))
        if _ion == "Ho2+":
            _n += "; DIVALENT Ho - not for Ho3+"
        if _ion == "Eu" and _h == "PbF2":
            _n += "; Abar6/Abar4 outside 0.28-0.32, handbook suspects error"
        SPM_TABLE.append(_row(_ion, "F", _h, None, _st(_a4, 4), _st(_a6, 6),
                              _T54_R[_h], norm=_S, page="p. 93, Table 5.4",
                              notes=_n))

SPM_TABLE.append(_row("Pr", "F", "LaF3", None, _st(80, 4), _st(32, 6), 2.44,
                      norm=_S, page="p. 93 (text, Yeung & Reid [YR89])",
                      notes="Abar4(2.44)=80, Abar6(2.44)=32"))

# ---- Table 5.6 (p. 95): Cl- in LaCl3, Wybourne ---------------------------
_T56 = dict(Pr=(42.2, 248, 319, 277), Nd=(41.6, 310, 335, 277),
            Pm=(42.0, 306, 378, 267), Eu=(41.7, 372, 282, 277),
            Tb=(41.8, 374, 284, 181), Ho=(41.8, 426, 279, 176),
            Er=(41.8, 438, 264, 162))
for _ion, (_th, _b2, _b4, _b6) in _T56.items():
    SPM_TABLE.append(_row(_ion, "Cl", "LaCl3", float(_b2), float(_b4),
                          float(_b6), 2.95, page="p. 95, Table 5.6 [NN89b]",
                          notes="theta=%.1f deg; R0 not printed in table, "
                                "2.95 A from p. 114" % _th))

# ---- Table 5.7 (p. 95): Cl-, Stevens, two hosts + t4, t6 [RR85] ----------
_T57 = dict(Pr=(82, 22, 40, 17, 16, 6), Nd=(70, 24, 40, 18, 13, 7),
            Sm=(58, 23, 32, 17, 14, 8), Tb=(58, 14, 34, 11, 12, 5),
            Ho=(60, 15, 33, 11, 13, 6), Er=(53, 15, 32, 10, 12, 9))
for _ion, (_a4i, _a6i, _a4ii, _a6ii, _t4, _t6) in _T57.items():
    SPM_TABLE.append(_row(_ion, "Cl", "Cs2NaLnCl6 / Cs2NaYCl6 (cubic)", None,
                          _st(_a4i, 4), _st(_a6i, 6), None, t4=float(_t4),
                          t6=float(_t6), norm=_S,
                          page="p. 95, Table 5.7 set (i) [RR85]",
                          notes="R not given in handbook"))
    SPM_TABLE.append(_row(_ion, "Cl", "anhydrous LnCl3 (C3h)", None,
                          _st(_a4ii, 4), _st(_a6ii, 6), 2.95, t4=float(_t4),
                          t6=float(_t6), norm=_S,
                          page="p. 95, Table 5.7 set (ii) [RR85]",
                          notes="R not printed; 2.95 A assumed (LaCl3, p. 114)"))

SPM_TABLE += [
    _row("Pr", "Cl", "LaCl3 (pressure)", None, _st(30, 4), _st(18, 6), 2.95,
         t4=6.0, t6=5.5, norm=_S, page="p. 114 [TGH93]",
         notes="same for Nd:LaCl3; Abar4=30(4), t4=6(2), Abar6=18(2), t6=5.5(2)"),
    _row("Nd", "Cl", "LaCl3 (pressure)", None, _st(30, 4), _st(18, 6), 2.95,
         t4=6.0, t6=5.5, norm=_S, page="p. 114 [TGH93]", notes="as Pr"),
    _row("Pr", "Cl", "LaCl3 (pressure, GdSH89)", None, _st(38, 4), _st(16, 6),
         2.95, t4=7.0, t6=7.5, norm=_S, page="p. 114 [GdSH89]",
         notes="Abar4=38(2), t4=7(2) without distortion correction; "
               "Abar6=16(2), t6=7.5(1.5), d6=-47(9) cm^-1/A"),
]

# ---- Table 5.8 (p. 96): Eu3+ with O2- / S2-, Wybourne ---------------------
for _h, _lig, _R, _b4, _b6, _ref in [
        ("La2O2S", "O", 2.42, 706, 302, "MN73"),
        ("La2O2S", "S", 3.04, 150, 82, "MN73"),
        ("Gd2O2S", "O", 2.27, 713, 315, "NS71"),
        ("Gd2O2S", "S", 2.95, 389, 87, "NS71"),
        ("YVO4", "O", 2.24, 400, 432, "NS71"),
        ("YGG", "O", 2.34, 475, 493, "NS69"),
        ("LaAlO3", "O", 2.68, 280, 262, "LL75")]:
    SPM_TABLE.append(_row("Eu", _lig, _h, None, float(_b4), float(_b6), _R,
                          page="p. 96, Table 5.8 [%s]" % _ref))

# ---- Table 5.14 (p. 112) + Table 5.13 (p. 110): garnets, Stevens ----------
_R_GARNET = dict(YAG=(2.303, 2.432), YGG=(2.338, 2.428))
_T514 = {  # label: (ion, host, A2(1), A2(2), A4(1), A4(2), A6(1), A6(2))
    "Nd:YAG": ("Nd", "YAG", 2901.7, 2820.4, 144.6, 74.1, 47.0, 33.6),
    "ErGG": ("Er", "ErGG (Er3Ga5O12)", -350.8, -371.4, 90.9, 54.6, 24.1, 20.5),
    "Er:YGG": ("Er", "YGG", -75.5, -100.0, 95.2, 57.3, 23.6, 21.1),
    "Er:YAG": ("Er", "YAG", 2167.5, 2096.3, 110.6, 67.0, 25.2, 18.6),
    "Dy:YGG": ("Dy", "YGG", -118.3, -157.1, 104.9, 74.1, 28.4, 21.5),
    "Dy:YAG": ("Dy", "YAG", 2474.7, 2395.3, 119.4, 70.0, 29.2, 22.6),
}


def _tk(b1, b2, r1, r2):
    """Eq. (5.6), p. 88."""
    import math
    return math.log(b1 / b2) / math.log(r2 / r1)


GARNET_T = {}
for _lab, (_ion, _h, _a21, _a22, _a41, _a42, _a61, _a62) in _T514.items():
    _key = "YAG" if "YAG" in _h else "YGG"
    _r1, _r2 = _R_GARNET[_key]
    _t4 = round(_tk(_a41, _a42, _r1, _r2), 2)
    _t6 = round(_tk(_a61, _a62, _r1, _r2), 2)
    GARNET_T[_lab] = (_t4, _t6)
    _rn = ("" if _h != "ErGG (Er3Ga5O12)" else
           "; ErGG distances not given, YGG R1/R2 assumed for t_k")
    for _i, (_a2, _a4, _a6, _R) in enumerate(
            [(_a21, _a41, _a61, _r1), (_a22, _a42, _a62, _r2)], start=1):
        SPM_TABLE.append(_row(
            _ion, "O", "%s shell i=%d" % (_h, _i), _st(_a2, 2), _st(_a4, 4),
            _st(_a6, 6), _R, t4=_t4, t6=_t6, norm=_S,
            page="p. 112, Table 5.14; R from Table 5.13 p. 110",
            notes="Bbar2 'unphysical' per handbook; t4/t6 from the two shells "
                  "via eq. (5.6) (Nd quoted 12.2/6.2 p. 112; others computed "
                  "here)" + _rn))

SPM_TABLE += [
    _row("Nd", "O", "YAG (Newman & Edgar, distortion-corrected)",
         _st(573, 2), None, None, 2.303, norm=_S, page="p. 112 [NE76]",
         notes="Abar2 = 573 cm^-1, nearest neighbours"),
    _row("Er", "O", "YAG (Newman & Edgar, distortion-corrected)",
         _st(451, 2), None, None, 2.303, norm=_S, page="p. 112 [NE76]",
         notes="Abar2 = 451 cm^-1, nearest neighbours"),
    _row("Er", "O", "YVO4 (sec. 5.4.1 inversion), shell 1",
         None, _st(58.8, 4), _st(23.7, 6), 2.29, norm=_S,
         page="p. 106, eq. (5.22); R Table 5.11 p. 101",
         notes="Abar4(R1)=58.8, Abar6(R1)=23.7"),
    _row("Er", "O", "YVO4 (sec. 5.4.1 inversion), shell 2",
         None, _st(18.7, 4), _st(22.3, 6), 2.43, norm=_S,
         page="p. 106, eq. (5.22)",
         notes="Abar4(R2)=18.7, Abar6(R2)=22.3; YPO4 inversion ill-conditioned "
               "(A4 = 185/-1854, rejected by handbook; A4 = 51+-2 preferred)"),
    _row("Ho/Tb", "O", "phosphate glass (semiclassical fit, ch. 9)",
         None, 272.0, 272.0, None, page="p. 181 (K. S. Chan, ch. 9)",
         notes="8-fold cubic average site, Bbar4 = Bbar6 = 272 cm^-1; smaller "
               "than Er:YVO4 as glass distances are larger"),
]


# ---------------------------------------------------------------------------
# 2. Per-ion model recommendations (Wybourne), see module docstring recipe.
# ---------------------------------------------------------------------------
def _scale(B, R_from, R_to, t):
    return B * (R_from / R_to) ** t


def _lin(z, z1, v1, z2, v2):
    return v1 + (v2 - v1) * (z - z1) / (z2 - z1)


def _entry(R0, B2, B4, B6, t2, t4, t6, src):
    r = lambda x: round(float(x), 1)
    return dict(R0=R0, B2=r(B2), B4=r(B4), B6=r(B6), t2=round(t2, 2),
                t4=round(t4, 2), t6=round(t6, 2), src=src)


def _build_O():
    R0 = 2.303
    t2 = 7.0
    base = {  # ion: (B2, B4, B6, t4, t6) at R0 = R1(YAG)
        "Nd": (1146.0, _st(144.6, 4), _st(47.0, 6)) + GARNET_T["Nd:YAG"],
        "Dy": (None, _st(119.4, 4), _st(29.2, 6)) + GARNET_T["Dy:YAG"],
        "Er": (902.0, _st(110.6, 4), _st(25.2, 6)) + GARNET_T["Er:YAG"],
    }
    # Dy Bbar2: linear in Z between NE76 Nd and Er values
    b2dy = _lin(66, 60, 1146.0, 68, 902.0)
    base["Dy"] = (b2dy,) + base["Dy"][1:]
    src_d = ("direct: Table 5.14 x8/x16 nearest O shell of {h} (p. 112), "
             "R0 = R1(YAG) Table 5.13; t4,t6 eq.(5.6) from the two shells; "
             "Bbar2 {b2}; t2 = 7 (p. 108)")
    out = {
        "Nd": _entry(R0, *base["Nd"][:3], t2, *base["Nd"][3:],
                     src_d.format(h="Nd:YAG", b2="= 2*573 [NE76] p. 112")),
        "Er": _entry(R0, *base["Er"][:3], t2, *base["Er"][3:],
                     src_d.format(h="Er:YAG", b2="= 2*451 [NE76] p. 112")),
        "Dy": _entry(R0, *base["Dy"][:3], t2, *base["Dy"][3:],
                     src_d.format(h="Dy:YAG",
                                  b2="INTERP (Z) between NE76 Nd and Er")),
    }

    def interp(z, a, b):
        za, zb = _Z[a], _Z[b]
        return [_lin(z, za, base[a][i], zb, base[b][i]) for i in range(5)]

    for ion, (a, b) in dict(Sm=("Nd", "Dy"), Gd=("Nd", "Dy"), Tb=("Nd", "Dy"),
                            Ho=("Dy", "Er")).items():
        v = interp(_Z[ion], a, b)
        out[ion] = _entry(R0, v[0], v[1], v[2], t2, v[3], v[4],
                          "INTERP linear in Z between YAG %s and %s "
                          "(Table 5.14); not a handbook value" % (a, b))
    out["Pr"] = dict(out["Nd"], src="ADOPTED Nd:YAG values (lighter Ln "
                     "'somewhat larger', p. 91; no Pr-O data in handbook)")
    for ion in ("Tm", "Yb"):
        out[ion] = dict(out["Er"], src="ADOPTED Er:YAG values (Er values are "
                        "a first approximation for >half-filled Ln, p. 91)")
    # Eu: direct Table 5.8 YGG, t's / B2 interpolated
    t4 = _lin(63, 60, base["Nd"][3], 66, base["Dy"][3])
    t6 = _lin(63, 60, base["Nd"][4], 66, base["Dy"][4])
    b2 = _scale(_lin(63, 60, 1146.0, 68, 902.0), R0, 2.34, t2)
    out["Eu"] = _entry(2.34, b2, 475.0, 493.0, t2, t4, t6,
                       "direct Bbar4/Bbar6: Table 5.8 Eu:YGG 2.34 A (p. 96) - "
                       "far below the Nd/Dy/Er garnet trend for Bbar4 "
                       "(Eu-O values span 280-713 across hosts); Bbar2, t4, t6 "
                       "INTERP (Z) from the YAG set, Bbar2 moved to 2.34 A")
    return out


def _build_F():
    tF = (3.4, 5.2, 5.1)        # ab initio Pr-F, p. 25
    R0 = 2.366                  # CaF2 nominal spacing
    # Cl- Bbar2 ion trend (Table 5.6), used to scale the single Er F value
    b2cl = {k: v[1] for k, v in _T56.items()}
    b2cl["Sm"] = _lin(62, 61, b2cl["Pm"], 63, b2cl["Eu"])
    b2cl["Gd"] = _lin(64, 63, b2cl["Eu"], 65, b2cl["Tb"])
    b2cl["Dy"] = _lin(66, 65, b2cl["Tb"], 67, b2cl["Ho"])
    b2cl["Tm"] = b2cl["Yb"] = b2cl["Er"]

    def b2F(ion, R):
        return _scale(480.0 * b2cl[ion] / b2cl["Er"], 2.42, R, tF[0])

    b2src = ("Bbar2 = Er:LaF3 480(60)@2.42 A (Table 5.2) x Cl- ion trend "
             "Bbar2_Cl(ion)/Bbar2_Cl(Er) (Table 5.6) - ASSUMPTION")
    tsrc = "t_k = ab initio Pr3+-F- values (p. 25), not empirical"
    out = {}
    direct = dict(Gd=("CaF2", 86.8, 27.9), Dy=("CaF2", 87.8, 25.8),
                  Er=("CaF2", 76.6, 22.9), Eu=("CdF2", 90.6, 28.1))
    for ion, (h, a4, a6) in direct.items():
        R = _T54_R[h]
        out[ion] = _entry(R, b2F(ion, R), _st(a4, 4), _st(a6, 6), *tF,
                          "direct Bbar4/Bbar6: Table 5.4 %s x8/x16 (p. 93), "
                          "R0 = nominal %s spacing; %s; %s"
                          % (ion, h, b2src, tsrc))
    out["Pr"] = _entry(2.44, b2F("Pr", 2.44), _st(80, 4), _st(32, 6), *tF,
                       "direct Bbar4/Bbar6: Pr:LaF3 Yeung & Reid (p. 93), "
                       "Bbar6/Bbar4 = 0.8 (outside fluorite 0.56-0.64); "
                       "%s; %s" % (b2src, tsrc))

    def at_R0(ion):
        e = out[ion]
        return (_scale(e["B4"], e["R0"], R0, tF[1]),
                _scale(e["B6"], e["R0"], R0, tF[2]))

    for ion, (a, b) in dict(Nd=("Pr", "Eu"), Sm=("Pr", "Eu"),
                            Tb=("Gd", "Dy"), Ho=("Dy", "Er")).items():
        va, vb = at_R0(a), at_R0(b)
        z = _Z[ion]
        B4 = _lin(z, _Z[a], va[0], _Z[b], vb[0])
        B6 = _lin(z, _Z[a], va[1], _Z[b], vb[1])
        out[ion] = _entry(R0, b2F(ion, R0), B4, B6, *tF,
                          "INTERP linear in Z between %s and %s after moving "
                          "both to R0 = 2.366 A with the power law; %s; %s"
                          % (a, b, b2src, tsrc))
    for ion in ("Tm", "Yb"):
        out[ion] = dict(out["Er"], src="ADOPTED Er:CaF2 values (p. 91 rule); "
                        + out["Er"]["src"])
    return out


def _build_Cl():
    R0 = 2.95
    t57 = {k: (float(v[4]), float(v[5])) for k, v in _T57.items()}
    t57["Eu"] = tuple(_lin(63, 62, t57["Sm"][i], 65, t57["Tb"][i])
                      for i in range(2))
    t57["Gd"] = tuple(_lin(64, 62, t57["Sm"][i], 65, t57["Tb"][i])
                      for i in range(2))
    t57["Dy"] = tuple(_lin(66, 65, t57["Tb"][i], 67, t57["Ho"][i])
                      for i in range(2))
    t2 = 7.0
    rsrc = ("R0 = 2.95 A (LaCl3, p. 114; Table 5.2's 2.42 A looks like a "
            "misprint); t2 = 7 (no Ln3+-Cl value; zircon O t2 = 7+-1 p. 108, "
            "An4+ Cl/Br t2 ~ 7 p. 96)")
    out = {}
    for ion in ("Pr", "Nd", "Eu", "Tb", "Ho", "Er"):
        _, b2, b4, b6 = _T56[ion]
        tsrc = ("t4,t6 Table 5.7 (p. 95)" if ion in _T57 else
                "t4,t6 INTERP (Z) from Table 5.7")
        out[ion] = _entry(R0, b2, b4, b6, t2, *t57[ion],
                          "direct Bbar_k: Table 5.6 LaCl3 (p. 95); %s; %s"
                          % (tsrc, rsrc))
    b2sm = _lin(62, 61, _T56["Pm"][1], 63, _T56["Eu"][1])
    out["Sm"] = _entry(R0, b2sm, _st(32, 4), _st(17, 6), t2, *t57["Sm"],
                       "direct Bbar4/Bbar6: Table 5.7 set (ii) x8/x16, t4,t6 "
                       "Table 5.7 (p. 95); Bbar2 INTERP Pm/Eu Table 5.6; " + rsrc)
    for ion, (a, b) in dict(Gd=("Eu", "Tb"), Dy=("Tb", "Ho")).items():
        z = _Z[ion]
        v = [_lin(z, _Z[a], _T56[a][i], _Z[b], _T56[b][i]) for i in (1, 2, 3)]
        out[ion] = _entry(R0, *v, t2, *t57[ion],
                          "INTERP linear in Z between Table 5.6 %s and %s; "
                          "t4,t6 INTERP from Table 5.7; %s" % (a, b, rsrc))
    for ion in ("Tm", "Yb"):
        out[ion] = dict(out["Er"], src="ADOPTED Er values (p. 91 rule); "
                        + out["Er"]["src"])
    return out


def _build_S(o_er):
    t2, t4, t6 = o_er["t2"], o_er["t4"], o_er["t6"]
    tsrc = ("t_k: NO S2- power laws in the handbook - copied from O2- Er:YAG "
            "entry as a PLACEHOLDER")
    er = _entry(2.34, 400.0, 376.0, 432.0, t2, t4, t6,
                "direct: Table 5.2 Er:ZnS (p. 91), uncertainties B2 200, "
                "B4 32, B6 96; " + tsrc)
    out = {}
    for ion in IONS:
        out[ion] = dict(er, src=("ADOPTED Er:ZnS values (only Ln3+-S2- set "
                                 "with Bbar2; p. 91 rule); " + tsrc))
    out["Er"] = er
    out["Eu"] = _entry(2.95, bbar2_contact(389.0, 87.0), 389.0, 87.0, t2, t4,
                       t6, "direct Bbar4/Bbar6: Table 5.8 Eu:Gd2O2S S2- 2.95 A "
                       "(p. 96; La2O2S gives 150/82 at 3.04 A); Bbar2 = eq. "
                       "(5.30) lower bound; " + tsrc)
    return out


_O = _build_O()
SPM_FOR_MODEL = {
    "O": {i: _O[i] for i in IONS},
    "F": {i: v for i, v in sorted(_build_F().items(), key=lambda x: _Z[x[0]])},
    "Cl": {i: v for i, v in sorted(_build_Cl().items(),
                                   key=lambda x: _Z[x[0]])},
    "S": _build_S(_O["Er"]),
}


def bbar(ligand, ion, k, R):
    """Wybourne Bbar_k(R) in cm^-1 from SPM_FOR_MODEL (eq. 5.5)."""
    e = SPM_FOR_MODEL[ligand][ion]
    return e["B%d" % k] * (e["R0"] / R) ** e["t%d" % k]

# Alternative to the direct (but trend-breaking) Eu-O entry: the same linear-
# in-Z interpolation within the YAG garnet set used for Sm/Gd/Tb.
SPM_ALTERNATIVES = {
    ("O", "Eu"): _entry(
        2.303,
        *[_lin(63, 60, a, 66, b) for a, b in zip(
            (_O["Nd"]["B2"], _O["Nd"]["B4"], _O["Nd"]["B6"]),
            (_O["Dy"]["B2"], _O["Dy"]["B4"], _O["Dy"]["B6"]))],
        7.0, _O["Eu"]["t4"], _O["Eu"]["t6"],
        "INTERP linear in Z between YAG Nd and Dy (Table 5.14); use instead "
        "of SPM_FOR_MODEL['O']['Eu'] for a smooth series"),
}
