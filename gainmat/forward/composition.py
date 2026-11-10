"""Composition front end: a glass composition (mol%) -> the model's glass family.

    from gainmat.forward.composition import family_for, parse_components
    family_for({"SiO2": 92, "Al2O3": 6, "Yb2O3": 2})      -> ("aluminosilicate", "rule: SiO2>=85, Al2O3 largest dopant")
    family_for("70TeO2-15K2O-15Nb2O5")                     -> ("tellurite", ...)

The rules use the dominant network former, as glass scientists classify these hosts; the model's
family then supplies the host parameters (phonon energy, coordination, disorder, refractive index,
Omegas).  Returns (None, reason) for glasses the model has no family for (borates, multicomponent
silicates with < 85 mol% SiO2 are mapped as 'proxy' where a neighbour family exists).
"""
from __future__ import annotations

import re

_SUB = str.maketrans("₀₁₂₃₄₅₆₇₈₉₋", "0123456789-")
RE_OXIDES = {"Er2O3", "Yb2O3", "Tm2O3", "Ho2O3", "Nd2O3", "Pr2O3", "Dy2O3", "Sm2O3", "Eu2O3", "Tb2O3", "Pr6O11",
             "Tb4O7", "ErF3", "YbF3", "TmF3", "HoF3", "NdF3", "PrF3", "DyF3", "ErCl3", "NdCl3"}


def parse_components(comp):
    """dict {name: mol%} from a dict, a LOMS 'A=x mol% | B=y mol%' string or a '70TeO2-15K2O' formula."""
    if isinstance(comp, dict):
        return {str(k).translate(_SUB): float(v) for k, v in comp.items()}
    s = str(comp).translate(_SUB)
    out = {}
    if "=" in s:
        for part in s.split("|"):
            m = re.match(r"\s*(.+?)\s*=\s*([\d.]+)", part)
            if m:
                out[m.group(1).strip()] = float(m.group(2))
        return out
    for m in re.finditer(r"([\d.]+)\s*([A-Z][A-Za-z0-9]*)", s):
        out[m.group(2)] = out.get(m.group(2), 0.0) + float(m.group(1))
    return out


def family_for(comp):
    c = {k: v for k, v in parse_components(comp).items() if k not in RE_OXIDES}
    tot = sum(c.values()) or 1.0
    x = {k: 100.0 * v / tot for k, v in c.items()}
    g = lambda *names: sum(x.get(n, 0.0) for n in names)
    if any(re.search(r"(S|Se)\d*$", k) and not k.endswith("O4") for k in x) and not any(k.endswith(("O2", "O3", "O5")) and x[k] > 30 for k in x):
        return "chalcogenide", "rule: sulfide/selenide former"
    fl = sum(v for k, v in x.items() if re.search(r"F\d*$", k) and not k.startswith(("Fe",)))
    if fl >= 50:
        if g("ZrF4", "HfF4") >= 30:
            return "ZBLAN", "rule: fluorozirconate (ZrF4 >= 30)"
        if g("InF3") >= 20:
            return "fluoroindate", "rule: InF3 >= 20"
        return "ZBLAN", "proxy: non-zirconate fluoride glass (fluorides >= 50 mol%)"
    if g("TeO2") >= 40:
        return "tellurite", "rule: TeO2 >= 40"
    if g("P2O5") >= 40 and g("SiO2") < 40:
        return "phosphate", "rule: P2O5 >= 40"
    if g("GeO2") >= 40 and g("GeO2") > g("SiO2"):
        return "germanate", "rule: GeO2 >= 40"
    if g("SiO2") >= 85:
        dop = {"aluminosilicate": g("Al2O3"), "germanosilicate": g("GeO2"), "phosphosilicate": g("P2O5")}
        best = max(dop, key=dop.get)
        if dop[best] < 0.5:
            return "silica", "rule: SiO2 >= 85, no former dopant >= 0.5 mol%"
        return best, f"rule: SiO2 >= 85, {best.replace('silicate', '')} dopant largest"
    if g("SiO2") >= 40:
        return "aluminosilicate", "proxy: multicomponent silicate (40 <= SiO2 < 85)"
    if g("GeO2") >= 25:
        return "germanate", "proxy: heavy-metal germanate (GeO2 >= 25)"
    return None, "no matching family (e.g. borate, bismuthate, aluminate)"
