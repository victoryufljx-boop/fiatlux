"""Reference data for the lanthanide 4f^N tier (gainmat.ln4f).

Every entry carries a `src` key. `verified=True` means the number was checked against
the cited source during development (web access to the source text); `verified=False`
means the value is quoted from widely-used literature from memory and MUST be checked
before being used as a calibration target. Never treat unverified values as ground truth.
"""

CARNALL_1989 = ("W. T. Carnall, G. L. Goodman, K. Rajnak, R. S. Rana, 'A systematic analysis of "
                "the spectra of the lanthanides doped into single crystal LaF3', J. Chem. Phys. 90, "
                "3443 (1989); doi:10.1063/1.455853. Values read from Table 4 of the ANL report version "
                "(OSTI 6995771). Bracketed (fixed) values flagged with *_fixed.")
QLANTH_2025 = ("arXiv:2503.17377 'Refining spectroscopic calculations for trivalent lanthanide ions' "
               "(qlanth), Table III, LaF3 refit of Carnall's data (used as cross-check).")

# Free-ion (+ host-averaged) parameters, Ln3+:LaF3, cm^-1. Carnall 1989 Table 4.
# The OSTI scan is OCR'd; values were cross-checked against the qlanth refit (all
# within ~0.5% for F^k, zeta; Pr gamma OCR flagged as uncertain by the reader).
FREE_ION = {
    "Pr3+": dict(N=2, F2=68878, F4=50347, F6=32901, zeta=751.7, alpha=16.23, beta=-566.6, gamma=1371,
                 src=CARNALL_1989, verified=True, note="gamma OCR uncertain (qlanth refit 1364)"),
    "Nd3+": dict(N=3, F2=73018, F4=52789, F6=35757, zeta=885.3, alpha=21.34, beta=-593.0, gamma=1445,
                 src=CARNALL_1989, verified=True),
    "Ho3+": dict(N=10, F2=94564, F4=66397, F6=52022, zeta=2145, alpha=17.15, beta=-607.9, gamma=1800,
                 gamma_fixed=True, src=CARNALL_1989, verified=True),
    "Er3+": dict(N=11, F2=97483, F4=67904, F6=54010, zeta=2376, alpha=17.79, beta=-582.1, gamma=1800,
                 gamma_fixed=True, src=CARNALL_1989, verified=True),
    "Tm3+": dict(N=12, F2=100134, F4=69613, F6=55975, zeta=2636, alpha=17.26, beta=-624.5, gamma=1820,
                 gamma_fixed=True, src=CARNALL_1989, verified=True),
}

# Tabulated squared reduced matrix elements (Carnall, Fields, Rajnak, J. Chem. Phys. 49,
# 4424 (1968) aquo-ion tables; reproduced in Kaminskii and in essentially every Er JO paper).
U_TABLE = {
    ("Er3+", "4I13/2", "4I15/2"): dict(U2=0.0195, U4=0.1173, U6=1.4316, verified=True,
                                       src="Carnall et al. JCP 49, 4424 (1968); quoted in task and in "
                                           "many Er JO papers"),
    ("Er3+", "4I11/2", "4I15/2"): dict(U2=0.0282, U4=0.0003, U6=0.3953, verified=False,
                                       src="Carnall 1968 (from memory)"),
    ("Er3+", "4I11/2", "4I13/2"): dict(U2=0.0331, U4=0.1708, U6=1.0864, verified=False,
                                       src="Weber, Phys. Rev. 157, 262 (1967) (from memory)"),
    ("Er3+", "4S3/2", "4I15/2"): dict(U2=0.0, U4=0.0, U6=0.2211, verified=False,
                                      src="Carnall 1968 (from memory)"),
    ("Er3+", "2H11/2", "4I15/2"): dict(U2=0.7125, U4=0.4123, U6=0.0925, verified=False,
                                       src="Carnall 1968 (from memory); hypersensitive"),
    ("Nd3+", "4F3/2", "4I9/2"): dict(U2=0.0, U4=0.2293, U6=0.0548, verified=False,
                                     src="Krupke, IEEE JQE 7, 153 (1971) (from memory)"),
    ("Nd3+", "4F3/2", "4I11/2"): dict(U2=0.0, U4=0.1423, U6=0.4083, verified=False,
                                      src="Krupke 1971 (from memory)"),
    ("Nd3+", "4F3/2", "4I13/2"): dict(U2=0.0, U4=0.0, U6=0.2140, verified=False,
                                      src="Krupke 1971 (from memory)"),
}

# Judd-Ofelt parameters, 1e-20 cm^2.
JUDD_OFELT = {
    "Nd:YAG": dict(O2=0.2, O4=2.7, O6=5.0, n=1.82, src="W. F. Krupke, IEEE J. Quantum Electron. 7, 153 "
                   "(1971) (values as commonly quoted)", verified=False,
                   note="NOT confirmed: UAB lecture notes (V. Fedorov, 'Judd Ofelt theory', "
                        "https://people.cas.uab.edu/~vfedorov/Notes/Judd%20Ofelt%20theory.pdf) give "
                        "O2=0.37, O4=2.29, O6=5.97 (unsourced). Krupke 1971 not read."),
    "Er:silica": dict(O2=4.23, O4=1.04, O6=0.61, n=1.45, src="Er-doped silicate/silica fiber, typical "
                      "literature values (from memory; spread across Al/Ge-silica is ~30%)", verified=False),
    "Er:ZBLAN": dict(O2=2.54, O4=1.39, O6=0.97, n=1.50, src="Er:ZBLAN, typical literature values "
                     "(from memory; e.g. Wetenkamp et al. J. Non-Cryst. Solids 140, 35 (1992))", verified=False),
    "Er:LiYF4": dict(O2=1.92, O4=0.61, O6=1.57, n=1.45, src="Er:YLF (from memory; orientation-averaged)",
                     verified=False),
}

# Experimental J-multiplet barycentres (cm^-1), Ln3+:LaF3 / Dieke diagram level
# positions, rounded. UNVERIFIED (from memory, ~+-100 cm^-1); intended as a coarse check.
EXP_LEVELS = {
    "Er3+": {"4I13/2": 6500, "4I11/2": 10150, "4I9/2": 12400, "4F9/2": 15250, "4S3/2": 18400,
             "2H11/2": 19150, "4F7/2": 20450, "4F5/2": 22150, "4F3/2": 22500, "2G9/2": 24550,
             "4G11/2": 26350},
    "Nd3+": {"4I11/2": 2000, "4I13/2": 3950, "4I15/2": 5950, "4F3/2": 11550, "4F5/2": 12550,
             "2H9/2": 12650, "4S3/2": 13500, "4F7/2": 13450, "4F9/2": 14750, "2H11/2": 15950,
             "4G5/2": 17300, "4G9/2": 19550, "2K13/2": 19650},
    "src": "Dieke diagram / Carnall 1989 LaF3 levels, rounded (from memory)",
    "verified": False,
}

# Measured lifetimes (s)
LIFETIMES = {
    ("Nd:YAG", "4F3/2"): dict(tau=230e-6, src="Nd:YAG fluorescence lifetime 230 us: RP Photonics "
                              "Encyclopedia, 'YAG lasers', Table 2 (https://www.rp-photonics.com/yag_lasers.html); "
                              "also quoted in Koechner, Solid-State Laser Engineering (not re-read)", verified=True),
    ("Er:silica", "4I13/2"): dict(tau=10e-3, src="RP Photonics Encyclopedia, 'Erbium-doped gain media': "
                                  "'upper-state lifetime of 4I13/2 is of the order of 8-10 ms' "
                                  "(https://www.rp-photonics.com/erbium_doped_gain_media.html); Lehigh OPG lecture 26 "
                                  "'~10 ms'. Desurvire (1994) not re-read.", verified=True),
    ("Er:silica", "4I11/2"): dict(tau=7e-6, src="few us (strongly quenched), typical", verified=False),
    ("Er:ZBLAN", "4I13/2"): dict(tau=9.0e-3, src="X. Zhu, N. Peyghambarian, Adv. OptoElectron. 2010, "
                                 "501956: 'lower laser level (4I13/2, tau1 = 9 ms)'", verified=True),
    ("Er:ZBLAN", "4I11/2"): dict(tau=6.9e-3, src="X. Zhu, N. Peyghambarian, Adv. OptoElectron. 2010, "
                                 "501956: 'upper laser level (4I11/2, tau2 = 6.9 ms)'", verified=True),
}
BRANCHING = {
    ("Nd:YAG", "4F3/2"): dict(values={"4I9/2": 0.32, "4I11/2": 0.53, "4I13/2": 0.14, "4I15/2": 0.005},
                              src="Krupke 1971 calculated branching (from memory; 0.25-0.32/0.53-0.60/0.12-0.14)",
                              verified=False),
}

# Energy-gap law W = B exp(-alpha dE): B in s^-1, alpha in cm, hw_max in cm^-1.
GAP_LAW = {
    "silicate":  dict(B=1.4e12, alpha=4.7e-3, hw=1100, src="Layne, Lowdermilk, Weber, PRB 16, 10 (1977)"),
    "phosphate": dict(B=5.4e12, alpha=4.7e-3, hw=1200, src="Layne, Lowdermilk, Weber, PRB 16, 10 (1977)"),
    "germanate": dict(B=3.4e10, alpha=4.9e-3, hw=900,  src="Layne, Lowdermilk, Weber, PRB 16, 10 (1977)"),
    "tellurite": dict(B=6.3e10, alpha=4.7e-3, hw=700,  src="Reisfeld & Jorgensen (1987) compilation"),
    "fluoride":  dict(B=1.88e10, alpha=5.77e-3, hw=500, src="Shinn, Sibley, Drexhage, Brown, PRB 27, 6635 "
                      "(1983), fluorozirconate (ZBLA)"),
    "chalcogenide": dict(B=1.0e6, alpha=2.9e-3, hw=350, src="GLS/chalcogenide, e.g. Schweizer et al. "
                         "(order of magnitude)"),
    "verified": False,
    "note": "Numbers recalled from the cited papers; check before calibration.",
}
VAN_DIJK = dict(beta_el=1e7, alpha=3.5e-3,
                src="J. M. F. van Dijk, M. F. H. Schuurmans, J. Chem. Phys. 78, 5317 (1983): "
                    "W = beta_el exp[-alpha (dE - 2 hw_max)]", verified=False)

HOSTS = {
    "YAG":    dict(family="garnet", n=1.82, hw=850),
    "silica": dict(family="silicate", n=1.45, hw=1100),
    "ZBLAN":  dict(family="fluoride", n=1.50, hw=600),  # Zhu & Peyghambarian 2010, Table 1 (verified)
    "LiYF4":  dict(family="fluoride", n=1.45, hw=560),
}


# ---------------------------------------------------------------------------------------
# Complete Carnall 1989 LaF3 parameter sets (Table 4), all operators, cm^-1.
# Source of the machine-readable numbers: qlanth repository (zia-lab/qlanth commit 14563fd,
# data/carnall.py "table 1 from Carnall et al. 1989" and data/carnallParams.m), cross-checked
# against the OSTI 6995771 ANL-report text of Table 4 (F^k, zeta, alpha, beta, gamma, T^i, M0,
# P2 agree for Nd/Er/Ho/Tm). M2=0.56 M0, M4=0.31 M0, P4=0.5 P2, P6=0.1 P2 (Table 4 footnotes).
FREE_ION_CARNALL = {
    "Pr3+": dict(N=2, F2=68878, F4=50347, F6=32901, zeta=751.7, alpha=16.23, beta=-566.6, gamma=1371,
                 T2=0, T3=0, T4=0, T6=0, T7=0, T8=0, M0=2.08, P2=-88.6, src=CARNALL_1989, verified=True),
    "Nd3+": dict(N=3, F2=73018, F4=52789, F6=35757, zeta=885.3, alpha=21.34, beta=-593.0, gamma=1445,
                 T2=298, T3=35, T4=59, T6=-285, T7=332, T8=305, M0=2.11, P2=192, src=CARNALL_1989, verified=True),
    "Ho3+": dict(N=10, F2=94564, F4=66397, F6=52022, zeta=2145, alpha=17.15, beta=-607.9, gamma=1800,
                 T2=400, T3=37, T4=107, T6=-264, T7=316, T8=336, M0=2.54, P2=605, src=CARNALL_1989, verified=True),
    "Er3+": dict(N=11, F2=97483, F4=67904, F6=54010, zeta=2376, alpha=17.79, beta=-582.1, gamma=1800,
                 T2=400, T3=43, T4=73, T6=-271, T7=308, T8=299, M0=3.86, P2=594, src=CARNALL_1989, verified=True),
    "Tm3+": dict(N=12, F2=100134, F4=69613, F6=55975, zeta=2636, alpha=17.26, beta=-624.5, gamma=1820,
                 T2=400, T3=0, T4=0, T6=0, T7=0, T8=0, M0=3.81, P2=695, src=CARNALL_1989, verified=True,
                 note="Tm T3..T8 not recoverable from OSTI OCR; qlanth file gives 0"),
}

# Experimental multiplet barycentres, Ln3+:LaF3 (cm^-1, relative to lowest Stark level):
# unweighted mean of the observed Kramers-doublet energies of each multiplet in
# Carnall, Goodman, Rajnak, Rana, J. Chem. Phys. 90, 3443 (1989) / ANL-88-8 (OSTI 6995771),
# Appendix II (Nd3+:LaF3) and Appendix VI (Er3+:LaF3), "Experimental and computed energy
# level structure". The OSTI scan text available here stops before the appendices, so the
# Stark levels were taken from the machine-readable transcription of those appendices in
# qlanth (zia-lab/qlanth commit 14563fd, data/Carnall.m, keys "appendix:Nd:RawTable",
# "appendix:Er:RawTable"; columns Exp/Calc/Exp-Calc). Only multiplets whose Stark levels in
# the appendix are contiguous under one label with the correct count (J+1/2) are kept;
# multiplets with interleaved labels (Nd 4F5/2/2H9/2, 2H11/2, 4G5/2-4G11/2, 2K15/2;
# Er 2H11/2, 4G11/2) are omitted. Levels with missing Exp entries are averaged over the
# observed ones (n_obs < n_states). Barycentres compared after subtracting the ground
# multiplet barycentre on both sides.
EXP_BARY_CARNALL1989 = {
    "Nd3+": {
        "4I9/2": 195.4,  # n_obs=5/5
        "4I11/2": 2097.3,  # n_obs=6/6
        "4I13/2": 4087.7,  # n_obs=7/7
        "4I15/2": 6164.6,  # n_obs=8/8
        "4F3/2": 11613.0,  # n_obs=2/2
        "4S3/2": 13673.5,  # n_obs=2/2
        "4F7/2": 13632.5,  # n_obs=4/4
        "4F9/2": 14894.4,  # n_obs=5/5
        "2G9/2": 21202.6,  # n_obs=5/5
        "2D3/2": 21345.5,  # n_obs=2/2
        "2P1/2": 23473.0,  # n_obs=1/1
        "2D5/2": 24034.7,  # n_obs=3/3
        "2P3/2": 26402.0,  # n_obs=2/2
        "4D3/2": 28357.5,  # n_obs=2/2
        "4D5/2": 28567.3,  # n_obs=3/3
        "4D1/2": 28962.0,  # n_obs=1/1
        "2I11/2": 29587.4,  # n_obs=5/6
    },
    "Er3+": {
        "4I15/2": 218.6,  # n_obs=8/8
        "4I13/2": 6700.6,  # n_obs=7/7
        "4I11/2": 10339.8,  # n_obs=6/6
        "4I9/2": 12596.6,  # n_obs=5/5
        "4F9/2": 15453.4,  # n_obs=5/5
        "4S3/2": 18572.5,  # n_obs=2/2
        "4F7/2": 20719.8,  # n_obs=4/4
        "4F5/2": 22383.7,  # n_obs=3/3
        "4F3/2": 22717.5,  # n_obs=2/2
        "2G9/2": 24747.6,  # n_obs=5/5
        "4G9/2": 27631.0,  # n_obs=5/5
        "2K15/2": 27912.0,  # n_obs=6/8
        "2G7/2": 28252.7,  # n_obs=3/4
    },
    "src": "Carnall et al. J. Chem. Phys. 90, 3443 (1989), ANL report OSTI 6995771, Appendix II (Nd) "
           "and Appendix VI (Er), via qlanth data/Carnall.m transcription (commit 14563fd)",
    "verified": True,
}
