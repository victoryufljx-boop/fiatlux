"""Yb3+ Stark levels at room temperature in eight glasses (Lorentzian decomposition of absorption and
emission).  Yang, Liu, Wang, Zhang, Hu & Zhang, Opt. Lett. 39, 1772 (2014), Fig. 2 (levels, cm^-1),
Table 1 (compositions, mol%), Table 2 (N_Yb, E_d, peak emission wavelength).  Numbers only.
2F7/2 levels a, b, c, d relative to a; 2F5/2 levels e, f, g absolute (above a)."""
YANG2014 = {
    "Ge": dict(family="germanate", proxy=False, comp="55GeO2-20PbO-15BaO-10ZnO-1.5Yb2O3",
               f72=[0, 316, 509, 901], f52=[10246, 10684, 11038], N=3.26e20, lam_p=1022),
    "Si": dict(family="aluminosilicate", proxy=True, comp="55SiO2-30PbO-15Al2O3-1.5Yb2O3 (lead aluminosilicate)",
               f72=[0, 265, 488, 874], f52=[10246, 10593, 10977], N=3.45e20, lam_p=1024),
    "B": dict(family=None, proxy=True, comp="58H3BO3-27BaO-15La2O3-1.5Yb2O3",
              f72=[0, 246, 452, 830], f52=[10246, 10537, 10917], N=2.72e20, lam_p=1018),
    "Bi": dict(family="germanate", proxy=True, comp="55Bi2O3-28H3BO3-10SiO2-7Ga2O3-1.5Yb2O3 (bismuthate)",
               f72=[0, 296, 509, 758], f52=[10246, 10493, 10846], N=2.35e20, lam_p=1024),
    "FP": dict(family=None, proxy=True, comp="15Al(PO3)3-10MgF2-10CaF2-17SrF2-35BaF2-8Ga2O3-5AlF3-3YbF3",
               f72=[0, 277, 511, 734], f52=[10267, 10504, 10846], N=4.63e20, lam_p=1004),
    "Te": dict(family="tellurite", proxy=False, comp="70TeO2-17WO3-8La2O3-5Na2O-1.5Yb2O3",
               f72=[0, 286, 509, 722], f52=[10246, 10504, 10834], N=2.76e20, lam_p=1006),
    "F": dict(family="ZBLAN", proxy=True, comp="54ZrF4-20BaF2-20NaF-3LaF3-3AlF3-3YbF3 (ZBNLA)",
              f72=[0, 257, 453, 624], f52=[10267, 10449, 10753], N=4.22e20, lam_p=1008),
    "P": dict(family="phosphate", proxy=False, comp="70P2O5-10H3BO3-10BaO-3Al2O3-3Nb2O5-4K2O-1.5Yb2O3",
              f72=[0, 207, 307, 596], f52=[10267, 10460, 10776], N=1.77e20, lam_p=1004),
}
SRC = "Yang et al., Opt. Lett. 39, 1772 (2014), Fig. 2 and Tables 1-2 (room temperature, Lorentzian fits)"
