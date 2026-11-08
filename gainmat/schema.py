"""Input/output records shared by every tier. Keep tiers decoupled: each tier
fills only its own fields and returns a dict that can be merged into Result."""
from __future__ import annotations
from dataclasses import dataclass, field, asdict
from typing import Optional


# ---------------------------------------------------------------- inputs
@dataclass
class Host:
    name: str                       # e.g. "LiYF4"
    composition: str                # formula or mole-fraction string
    kind: str = "crystal"           # "crystal" | "glass"
    structure_key: Optional[str] = None   # key into gainmat.structures registry
    quench: Optional[dict] = None   # glasses: {"T_melt":..,"rate_K_ps":..,"n_atoms":..}


@dataclass
class Dopant:
    ion: str                        # "Er3+"
    site: Optional[str] = None      # substituted cation, "Y"
    concentration_at_pct: float = 1.0


@dataclass
class Candidate:
    host: Host
    dopant: Optional[Dopant] = None
    target_nm: Optional[float] = None


# ---------------------------------------------------------------- outputs
@dataclass
class HostProps:
    gap_pbe_eV: Optional[float] = None
    gap_gllbsc_eV: Optional[float] = None
    uv_edge_nm: Optional[float] = None
    eps_inf: Optional[float] = None
    n_inf: Optional[float] = None
    phonon_max_cm1: Optional[float] = None      # highest optical phonon
    phonon_eff_cm1: Optional[float] = None      # DOS-weighted cutoff used in gap law
    method: dict = field(default_factory=dict)


@dataclass
class Level:
    label: str                      # "4I13/2"
    energy_cm1: float
    J: float


@dataclass
class Transition:
    upper: str
    lower: str
    wavelength_nm: float
    A_ed_s1: float = 0.0
    A_md_s1: float = 0.0
    branching: float = 0.0
    W_nr_s1: float = 0.0


@dataclass
class DopantProps:
    levels: list = field(default_factory=list)        # [Level]
    transitions: list = field(default_factory=list)   # [Transition]
    tau_rad_s: dict = field(default_factory=dict)     # level -> s
    quantum_eff: dict = field(default_factory=dict)   # level -> 0..1
    judd_ofelt: Optional[dict] = None                 # {"O2":..,"O4":..,"O6":..} 1e-20 cm2
    abs_peaks_nm: list = field(default_factory=list)  # 5d / 3d centers
    em_peak_nm: Optional[float] = None
    stokes_cm1: Optional[float] = None
    method: dict = field(default_factory=dict)


@dataclass
class Result:
    candidate: Candidate
    host: HostProps = field(default_factory=HostProps)
    dopant: DopantProps = field(default_factory=DopantProps)

    def to_dict(self):
        return asdict(self)


CM1_PER_EV = 8065.544
def ev_to_nm(e):  return 1239.842 / e
def cm1_to_nm(k): return 1e7 / k
