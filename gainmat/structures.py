"""Crystal-structure registry for the host tier.

get_structure(key, primitive=True) -> ase.Atoms

Every structure is built from space group + experimental lattice constants +
Wyckoff positions with ase.spacegroup.crystal, then (optionally) reduced to the
primitive cell with spglib.  Each entry cites its crystallographic source.
Each entry also stores the experimental density and the expected short
cation-anion bond lengths so that `verify_structure` can catch setting / origin
mistakes (the usual failure mode when typing Wyckoff positions).

NOTE on origin choices: ase's `setting=1/2` follows the International Tables
origin choices.  The entries below were checked numerically (bond lengths) on
build; see tests/test_host.py::test_structures.
"""
from __future__ import annotations

import numpy as np

# --------------------------------------------------------------------------
# key -> spec.  'bonds' = (cation, anion, expected nearest distance in A)
# --------------------------------------------------------------------------
STRUCTURES = {
    # alpha-quartz, P3_221 (No. 154).
    # Source: L. Levien, C.T. Prewitt, D.J. Weidner, "Structure and elastic
    # properties of quartz at pressure", Am. Mineral. 65, 920-930 (1980):
    # a = 4.9160, c = 5.4054 A; Si x = 0.4697; O (0.4135, 0.2669, 0.1191).
    # The z offsets below put the same coordinates into ase's setting of
    # P3_221 (z_Si = 2/3, z_O = 0.1191 + 2/3); verified to give 4 Si-O bonds
    # of 1.605/1.614 A (Levien et al.: 1.605(1)/1.614(1) A).
    "SiO2-alpha-quartz": dict(
        symbols=["Si", "O"],
        basis=[(0.4697, 0.0, 2.0 / 3.0), (0.4135, 0.2669, 0.1191 + 2.0 / 3.0)],
        spacegroup=154, setting=1,
        cellpar=[4.9160, 4.9160, 5.4054, 90, 90, 120],
        formula="SiO2", Z=3, density=2.648,
        bonds=[("Si", "O", 1.61)],
    ),
    # Fluorite CaF2, Fm-3m (No. 225), Ca 4a (0,0,0), F 8c (1/4,1/4,1/4).
    # a = 5.4626 A at room temperature: D.N. Batchelder, R.O. Simmons,
    # J. Chem. Phys. 41, 2324 (1964).  Density 3.18 g/cm3 (CRC Handbook).
    "CaF2": dict(
        symbols=["Ca", "F"],
        basis=[(0, 0, 0), (0.25, 0.25, 0.25)],
        spacegroup=225, setting=1,
        cellpar=[5.4626, 5.4626, 5.4626, 90, 90, 90],
        formula="CaF2", Z=4, density=3.18,
        bonds=[("Ca", "F", 2.365)],
    ),
    # LiYF4 scheelite, I4_1/a (No. 88), origin choice 1 in ase:
    # Li 4a, Y 4b, F 16f.
    # Source: E. Garcia, R.R. Ryan, "Structure of the laser host material
    # LiYF4", Acta Cryst. C49, 2053-2054 (1993).  a = 5.1668, c = 10.733 A,
    # F (0.2821, 0.1642, 0.0815).
    # CAVEAT: the coordinates were taken from memory of the Garcia & Ryan
    # values and checked only through the resulting bond lengths
    # (Li-F 1.90 A x4, Y-F 2.246/2.293 A, matching the published scheelite
    # LiYF4 geometry) and the density (3.99 g/cm3); the original CIF (COD
    # entry) could not be downloaded from this sandbox for a line-by-line check.
    # Lattice constants a = 5.1695(1), c = 10.7309(1) A were also reported
    # by Ma et al., J. Mater. Chem. C (2023) SI, Table 1 -- consistent.
    "LiYF4": dict(
        symbols=["Li", "Y", "F"],
        basis=[(0, 0, 0), (0, 0, 0.5), (0.2821, 0.1642, 0.0815)],
        spacegroup=88, setting=1,
        cellpar=[5.1668, 5.1668, 10.733, 90, 90, 90],
        formula="LiYF4", Z=4, density=3.99,
        bonds=[("Li", "F", 1.90), ("Y", "F", 2.25)],
    ),
    # alpha-Al2O3 corundum, R-3c (No. 167), hexagonal axes:
    # Al 12c (0,0,z) z = 0.35216, O 18e (x,0,1/4) x = 0.30624,
    # a = 4.7589, c = 12.991 A.
    # Source: J. Lewis, D. Schwarzenbach, H.D. Flack, Acta Cryst. A38,
    # 733-739 (1982) (coordinates); lattice constants as tabulated in the
    # ICSD/AMCSD corundum entries.  Density 3.987 g/cm3.
    "Al2O3": dict(
        symbols=["Al", "O"],
        basis=[(0, 0, 0.35216), (0.30624, 0, 0.25)],
        spacegroup=167, setting=1,
        cellpar=[4.7589, 4.7589, 12.991, 90, 90, 120],
        formula="Al2O3", Z=6, density=3.987,
        bonds=[("Al", "O", 1.855)],
    ),
    # Y3Al5O12 garnet, Ia-3d (No. 230):
    # Y 24c (1/8,0,1/4), Al(oct) 16a (0,0,0), Al(tet) 24d (3/8,0,1/4),
    # O 96h (-0.0306, 0.0512, 0.1500), a = 12.000 A.
    # Source: F. Euler, J.A. Bruce, "Oxygen coordinates of compounds with
    # garnet structure", Acta Cryst. 19, 971-978 (1965).  Density 4.55 g/cm3.
    # Conventional cell 160 atoms, primitive (bcc) cell 80 atoms.
    "Y3Al5O12": dict(
        symbols=["Y", "Al", "Al", "O"],
        basis=[(0.125, 0, 0.25), (0, 0, 0), (0.375, 0, 0.25),
               (-0.0306, 0.0512, 0.1500)],
        spacegroup=230, setting=1,
        cellpar=[12.000, 12.000, 12.000, 90, 90, 90],
        formula="Al5O12Y3", Z=8, density=4.55,
        bonds=[("Al", "O", 1.76), ("Y", "O", 2.30)],
    ),
}

ALIASES = {
    "quartz": "SiO2-alpha-quartz", "SiO2": "SiO2-alpha-quartz",
    "alpha-SiO2": "SiO2-alpha-quartz", "a-SiO2": "SiO2-alpha-quartz",
    "YLF": "LiYF4", "sapphire": "Al2O3", "corundum": "Al2O3",
    "YAG": "Y3Al5O12", "fluorite": "CaF2",
}


def list_structures():
    return list(STRUCTURES)


def _resolve(key: str) -> str:
    if key in STRUCTURES:
        return key
    if key in ALIASES:
        return ALIASES[key]
    raise KeyError(f"unknown structure key {key!r}; known: {list(STRUCTURES)}")


def get_structure(key: str, primitive: bool = True):
    """Return an ase.Atoms for `key` (conventional cell if primitive=False)."""
    from ase.spacegroup import crystal
    spec = STRUCTURES[_resolve(key)]
    atoms = crystal(spec["symbols"], spec["basis"], spacegroup=spec["spacegroup"],
                    setting=spec["setting"], cellpar=spec["cellpar"],
                    primitive_cell=False, size=(1, 1, 1))
    if primitive:
        atoms = to_primitive(atoms)
    atoms.info["structure_key"] = _resolve(key)
    return atoms


def to_primitive(atoms, symprec=1e-3):
    """Standardised primitive cell via spglib (keeps ase.Atoms type)."""
    import spglib
    from ase import Atoms
    cell = (atoms.cell[:], atoms.get_scaled_positions(), atoms.numbers)
    lat, pos, nums = spglib.standardize_cell(cell, to_primitive=True,
                                             no_idealize=False, symprec=symprec)
    return Atoms(numbers=nums, cell=lat, scaled_positions=pos, pbc=True)


def density_gcm3(atoms) -> float:
    return atoms.get_masses().sum() * 1.66053906660 / atoms.get_volume()


def nearest_distance(atoms, a: str, b: str, rc: float = 4.0) -> float:
    from ase.neighborlist import neighbor_list
    i, j, d = neighbor_list("ijd", atoms, rc)
    s = np.array(atoms.get_chemical_symbols())
    m = (s[i] == a) & (s[j] == b)
    return float(d[m].min())


def verify_structure(key: str, tol_density: float = 0.03, tol_bond: float = 0.06):
    """Check formula, density and shortest cation-anion distances.

    Returns a dict; raises AssertionError on failure."""
    k = _resolve(key)
    spec = STRUCTURES[k]
    conv = get_structure(k, primitive=False)
    prim = get_structure(k, primitive=True)
    from ase.formula import Formula
    red = Formula(conv.get_chemical_formula()).reduce()[0]
    want = Formula(spec["formula"]).reduce()[0]
    assert red == want, f"{k}: formula {red} != {want}"
    assert len(conv) % len(prim) == 0
    rho = density_gcm3(conv)
    assert abs(rho - spec["density"]) / spec["density"] < tol_density, \
        f"{k}: density {rho:.3f} vs exp {spec['density']}"
    bonds = {}
    for a, b, d0 in spec["bonds"]:
        d = nearest_distance(conv, a, b)
        assert abs(d - d0) < tol_bond, f"{k}: {a}-{b} {d:.3f} vs {d0}"
        bonds[f"{a}-{b}"] = round(d, 4)
    # no unphysically close pairs of any kind
    from ase.neighborlist import neighbor_list
    dmin = float(neighbor_list("d", conv, 3.0).min())
    assert dmin > 1.3, f"{k}: atoms too close ({dmin:.2f} A)"
    return dict(key=k, formula=str(red), n_conv=len(conv), n_prim=len(prim),
                density=round(rho, 4), density_exp=spec["density"],
                bonds=bonds, dmin=round(dmin, 4))


if __name__ == "__main__":
    for key in STRUCTURES:
        print(verify_structure(key))
