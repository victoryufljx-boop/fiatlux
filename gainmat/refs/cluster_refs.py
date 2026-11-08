"""Experimental reference data for the embedded-cluster tier (gainmat.cluster).

Every entry carries a citation and a `verified` flag:
  verified=True  : number read from the cited source (or a paper quoting it, named in `via`)
                   during development with web access.
  verified=False : textbook/commonly quoted value; check the original before relying on it.
"""

# ------------------------------------------------------------------ structures
STRUCTURES = {
    "LiYF4": dict(
        spacegroup=88, setting=1,        # I4_1/a, origin choice 1 (Li 4a, Y 4b)
        a=5.1668, c=10.733,
        sites=[("Li", (0.0, 0.0, 0.0), +1), ("Y", (0.0, 0.0, 0.5), +3),
               ("F", (0.2820, 0.1642, 0.0815), -1)],
        cite="E. Garcia, R. R. Ryan, Acta Cryst. C 49, 2053 (1993) (scheelite LiYF4)",
        verified=False,
        note="Coordinates as commonly quoted; reproduce Y-F = 2.246 x4, 2.293 x4 A, "
             "consistent with published Y-F 2.24/2.29 A.",
    ),
    "Al2O3": dict(
        spacegroup=167, setting=1,       # R-3c, hexagonal axes
        a=4.7602, c=12.9933,
        sites=[("Al", (0.0, 0.0, 0.35216), +3), ("O", (0.30624, 0.0, 0.25), -2)],
        cite="A. Kirfel, K. Eichhorn, Acta Cryst. A 46, 271 (1990) (corundum, 295 K)",
        verified=False,
        note="Gives Al-O = 1.855 x3, 1.971 x3 A (literature 1.86/1.97 A).",
    ),
}

# ------------------------------------------------------------------ Ce3+:LiYF4
CE_LIYF4 = dict(
    # five 4f -> 5d crystal-field components (low temperature excitation spectra)
    abs_5d_cm1=[33433, 41101, 48564, 50499, 52790],
    abs_cite="L. van Pieterson, M. F. Reid, R. T. Wegh, S. Soverna, A. Meijerink, "
             "Phys. Rev. B 65, 045113 (2002); also J. Chem. Phys. 115, 9382 (2001)",
    abs_via="L. Hu, M. F. Reid, C. K. Duan, S. Xia, M. Yin, arXiv:1011.1709, Table 1 "
            "(experimental column, footnote a -> van Pieterson et al.)",
    abs_verified=True,
    # 5d -> 4f emission doublet (to 2F5/2 and 2F7/2), room-temperature laser emission
    em_nm=[306.0, 325.0],
    em_cite="D. J. Ehrlich, P. F. Moulton, R. M. Osgood Jr., Opt. Lett. 4, 184 (1979) "
            "('Ultraviolet solid-state Ce:YLF laser at 325 nm')",
    em_verified=False,
    # derived: lowest 5d absorption minus 5d->2F5/2 emission maximum
    stokes_cm1=33433 - 1e7 / 306.0,
    stokes_note="derived from the two rows above (~750 cm^-1); not an independent measurement",
)

# Ce3+ spin-orbit parameters from fits of LiYF4:Ce3+ 4f and 5d levels
CE_ZETA = dict(
    zeta_4f_cm1=615.0, zeta_5d_cm1=1082.0,
    cite="L. van Pieterson et al., J. Chem. Phys. 115, 9382 (2001) / PRB 65, 045113 (2002)",
    via="arXiv:1011.1709 Table 2 (experimental column)",
    verified=True,
)

# ------------------------------------------------------------------ Ti3+:Al2O3
TI_AL2O3 = dict(
    abs_nm=[490.0, 560.0],   # main 2T2g->2Eg band and Jahn-Teller shoulder
    abs_cite="P. F. Moulton, J. Opt. Soc. Am. B 3, 125 (1986)",
    abs_verified=False,
    abs_note="Band maximum ~490 nm, shoulder ~550-560 nm (pi/sigma dependent). "
             "Dynamic Jahn-Teller splitting of 2Eg ~ 2000-3000 cm^-1. "
             "Moulton 1986 full text NOT read (abstract gives no peak values). "
             "Secondary: RP Photonics 'Titanium-sapphire lasers' ('absorption peak at "
             "~490 nm'); vendor datasheets (Red Optronics, Crylink) give 488 nm.",
    em_nm=780.0,
    em_cite="P. F. Moulton, J. Opt. Soc. Am. B 3, 125 (1986)",
    em_verified=False,
    em_note="Broad emission peak quoted 780-800 nm depending on polarisation/temperature. "
            "Moulton 1986 full text NOT read. Secondary: Crylink datasheet ~780 nm, "
            "Red Optronics datasheet 795 nm; RP Photonics gives max gain ~800 nm.",
    stokes_cm1=1e7 / 490.0 - 1e7 / 780.0,
    stokes_note="derived (band-max difference), not an independent measurement",
)

TI_ZETA = dict(zeta_3d_cm1=154.0,
               cite="Free-ion Ti3+ value, e.g. A. Abragam, B. Bleaney, "
                    "Electron Paramagnetic Resonance of Transition Ions (1970), Table 7.6",
               verified=False)

def summary():
    return dict(CE_LIYF4=CE_LIYF4, CE_ZETA=CE_ZETA, TI_AL2O3=TI_AL2O3, TI_ZETA=TI_ZETA)
