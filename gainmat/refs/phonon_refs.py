"""Experimental zone-centre (Gamma) phonons, eps_inf and Born charges for the
phonon bake-off (BAKEOFF_PHONONS.md).  All values in cm^-1.

Rules (same as host_refs.py): every number was seen on the page given in
`url` during the bake-off build (Nov 2026).  Where the page is a secondary
source quoting a primary paper, both are named.  Primary papers that were
paywalled and NOT opened are named as "cites ...".  Lists per irrep are in
ascending order but may be INCOMPLETE (missing members were not found on a
readable page); the comparison code matches each experimental value to the
nearest unused computed mode of the same irrep, so gaps are allowed.

Keys
  raman   : {irrep: [freq, ...]}              Raman-active modes
  ir_to   : {irrep: [TO, ...]}                IR-active transverse modes
  ir_lo   : {irrep: [LO, ...]}                IR-active longitudinal modes
  eps_inf : scalar or (perp, par) w.r.t. the c axis
  born    : {element: Z*} (scalar, isotropic) -- may be DFT, see note
"""
from __future__ import annotations

PHONON_REFS = {
    "CaF2": dict(
        raman={"T2g": [322]},
        ir_to={"T1u": [257]},
        ir_lo={"T1u": [463]},
        eps_inf=2.045,
        born={"Ca": 2.370, "F": -1.185},
        source="Y. Wang, L.A. Zhang, S. Shang, Z.-K. Liu, L.-Q. Chen, 'Accurate "
               "calculations of phonon dispersion in CaF2 and CeO2', Phys. Rev. B "
               "88, 024304 (2013), Table I 'Experiment' column (cites Kaiser et al. "
               "Phys. Rev. 127, 1950 (1962) and Russell, Proc. Phys. Soc. 85, 194 "
               "(1965), not opened). Z*(Ca)=2.370 is that paper's PBE value (no "
               "experimental Z* given); Z*(F) = -Z*(Ca)/2 by charge neutrality.",
        url="https://personal.ems.psu.edu/~chen/publications/Wang%20et%20al._2013_"
            "Physical%20Review%20B_Accurate%20calculations%20of%20phonon%20"
            "dispersion%20in%20CaF2%20and%20CeO2.pdf",
        cross_check="J.R. Love et al., APS March Meeting 2022 abstract T00.279 "
                    "(FTIR ellipsometry, 300 K): TO 261, LO 477, eps_inf 1.98. "
                    "https://archive.aps.org/mar/2022/t00/279",
    ),
    "Al2O3": dict(
        # Eg list is incomplete: 429 (labelled Eg(2)) and 748 (Eg(5)) from
        # arXiv:2602.19617 (cites Thapa et al. Appl. Opt. 56, 8598 (2017));
        # 576.3 labelled Eg in PMC9693301.  A1g 417 from arXiv:2602.19617.
        raman={"A1g": [417], "Eg": [429, 576.3, 748]},
        ir_to={"A2u": [397.5, 582.41], "Eu": [384.9, 439.1, 569.0, 633.6]},
        ir_lo={"A2u": [510.9, 881.1], "Eu": [387.6, 481.7, 629.5, 906.6]},
        eps_inf=(3.077, 3.072),
        born={"Al": 2.92, "O": -1.94},
        source="IR: M. Stokey et al., Phys. Rev. Materials 6, 014601 (2022), "
               "arXiv:2108.12058, Tables XII-XIII, 'same data for alpha-Al2O3 "
               "determined previously by Schubert et al.' (M. Schubert, T.E. "
               "Tiwald, C.M. Herzinger, PRB 61, 8187 (2000), not opened). "
               "Raman: arXiv:2602.19617 and PMC9693301 (Micromachines 13, 1961). "
               "Born charges: scalar bulk values quoted in arXiv:1712.03915 "
               "(origin, DFT vs. other, not stated on the page).",
        url="https://arxiv.org/pdf/2108.12058 ; https://arxiv.org/pdf/2602.19617 ; "
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC9693301/ ; "
            "https://ar5iv.arxiv.org/html/1712.03915",
        cross_check="Zeidler, Posch & Mutschke via arXiv:1304.1717 (300 K): Eu "
                    "TO 384.0/439.2/569.4/634.4, LO 387.5/482.0/908.2/630.6.",
    ),
    "SiO2-alpha-quartz": dict(
        raman={"A1": [207, 355, 465, 1080]},
        ir_to={"E": [128, 264.13, 393.55, 449.78, 694.31, 795.36, 1064.57, 1160.03],
               "A2": [363.03, 494.99, 775.58, 1073.37]},
        ir_lo={"E": [265.36, 402.88, 507.89, 697.31, 808.72, 1226.34, 1157.02],
               "A2": [386.42, 552.82, 789.17, 1238.42]},
        eps_inf=(2.356, 2.383),
        born=None,
        source="IR (300 K): Zeidler, Posch & Mutschke as tabulated in "
               "arXiv:1304.1717 Tables 4-5. E 128 cm-1 and A1 207/355/465: "
               "arXiv:2602.19617 (cites Masso et al. PRB 1970); A1 1080 "
               "('around'): Zhong et al., Contrib. Mineral. Petrol. 176:89 (2021). "
               "eps_inf exp: arXiv:cond-mat/0308077 Table III. Born charges: not "
               "found on a readable page -> no NAC for quartz.",
        url="https://arxiv.org/pdf/1304.1717 ; https://arxiv.org/pdf/2602.19617 ; "
            "https://link.springer.com/content/pdf/10.1007/s00410-021-01845-x.pdf ; "
            "https://arxiv.org/pdf/cond-mat/0308077",
        cross_check="Winta, Wolf, Paarmann PRB 99, 144308 (2019), 1.5 K: A2 TO "
                    "360.7/497.9/773.7/1073.0 LO 384.8/553.6/789.9/1238.7.",
    ),
    "LiYF4": dict(
        raman={"Ag": [151, 265, 427], "Bg": [174, 246, 327, 379, 427],
               "Eg": [154, 199, 326, 373, 447]},
        ir_to={"Au": [198, 251, 339, 372], "Eu": [137, 294, 323, 418]},
        ir_lo={},
        eps_inf=None,
        born=None,
        source="Salauen et al., J. Phys.: Condens. Matter 9, 6941 (1997) (300 K "
               "Raman + IR reflectivity; not opened) as reproduced in Tables 4-5 of "
               "arXiv:cond-mat/0504157. CAVEAT: the PDF text extraction lost the "
               "column alignment; mapping inferred (see BAKEOFF_PHONONS.md). No LO, "
               "eps_inf or Born charges found -> no NAC for LiYF4.",
        url="https://arxiv.org/pdf/cond-mat/0504157",
    ),
}


def match_modes(computed, exp_by_irrep):
    """Match experimental frequencies to computed modes of the same irrep.

    computed : list of dict(freq_cm1, irrep, deg)
    returns list of (irrep, exp, calc or None)."""
    pairs = []
    for irrep, flist in exp_by_irrep.items():
        pool = [m["freq_cm1"] for m in computed if m["irrep"] == irrep]
        cand = sorted(((abs(c - e), i, j) for i, e in enumerate(flist)
                       for j, c in enumerate(pool)))
        used_e, used_c, got = set(), set(), {}
        for _, i, j in cand:
            if i in used_e or j in used_c:
                continue
            used_e.add(i); used_c.add(j); got[i] = pool[j]
        for i, e in enumerate(flist):
            pairs.append((irrep, e, got.get(i)))
    return pairs
