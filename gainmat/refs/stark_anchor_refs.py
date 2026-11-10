"""Measured Stark levels used by the measured-level mode (spectra.ensemble_spectra(level_anchor=...)).
Numbers only, with their sources.  Levels are cm^-1 above the multiplet's lowest level, one entry per
Kramers doublet for Er/Yb (None = not resolved); E0 = lowest level of the multiplet above the lowest
ground level.  Each anchor comes from ONE source and is tested on data from OTHER sources
(scripts/validate_level_anchor.py)."""
ANCHORS = {
    "Yb3+": dict(
        levels={"2F7/2": dict(levels=[0, 425, 530, 865]), "2F5/2": dict(levels=[0, 630], E0=1e7 / 977.0)},
        families=("silica", "aluminosilicate", "germanosilicate"),
        src="Newell, Peterson, Gavrielides & Sharma, Opt. Commun. 273, 256 (2007): 77 K emission lines "
            "977, 1020, 1032, 1069 nm; 2F5/2 gap 630 cm^-1 (Yb silica fibre, composition not stated)"),
    "Er3+": dict(
        levels={"4I15/2": dict(levels=[0, 30, 55, 129, 234, None, None, 420]),
                "4I13/2": dict(levels=[0, None, None, None, None, None, 330])},
        families=("silica", "aluminosilicate", "germanosilicate"),
        src="Desurvire & Simpson, Opt. Lett. 15, 547 (1990) (4.2 K FLN, silicate glasses): 4I15/2 lowest "
            "five levels and total splitting 395-500 (420 used); 4I13/2 total 305-357 (330 used)"),
    "Tm3+": dict(
        levels={"3H6": dict(levels=[0, 250, 430])},
        families=(), auto=False,     # partial (3 of 13 levels): tested, makes spectra worse, not used
        src="Lincoln et al., J. Lumin. 50, 297 (1991): 9 K FLN in Ge-silica, ground 3H6 components "
            "~0, ~250, ~430 cm^-1 (more above ~500, not resolved)"),
}
