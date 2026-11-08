# gainmat — inputs, outputs, and how each output is computed

Goal: predict whether a (host, dopant) pair can amplify light at a target wavelength,
from first principles, so an ML model can later be trained on the simulated outputs
and calibrated against a small experimental set.

## 1. Inputs (what defines a candidate)

| Field | Example | Notes |
|---|---|---|
| `host.composition` | `LiYF4`, `0.53ZrF4-0.20BaF2-0.04LaF3-0.03AlF3-0.20NaF` | Formula or oxide/fluoride mole fractions for glasses |
| `host.kind` | `crystal` / `glass` | Glasses need an MD melt-quench to get a structure ensemble |
| `host.structure` | space group + lattice + Wyckoff sites, or CIF | Crystals only |
| `host.quench` | T_melt, cooling rate, N atoms, seeds | Glasses only |
| `dopant.ion` | `Er3+`, `Nd3+`, `Ce3+`, `Ti3+` | Element + oxidation state |
| `dopant.site` | `Y`, `Ca`, `Al` | Substituted cation |
| `dopant.concentration` | 1 at.% | Enters cross-relaxation / quenching (later tier) |
| `codopants` | `Yb3+` sensitizer | Optional; energy transfer is a later tier |
| `target` | wavelength band, pump wavelength | Used to rank, not to compute |

## 2. Outputs

### Host (dopant-independent)
| Output | Why it matters | Method (tier) |
|---|---|---|
| Band gap E_g → UV edge λ_UV = 1240/E_g nm | Sets shortest usable wavelength; 5d levels must sit below the conduction band | GPAW GLLB-SC (cheap, good for wide-gap ionics); PBE reported for reference |
| ε∞ → refractive index n | Local-field factor in radiative rates; waveguide design | GPAW RPA (scissor-corrected) |
| Max phonon energy ħω_max, phonon DOS | Multiphonon quenching → sets longest usable wavelength | MACE foundation MLIP + phonopy; spot-checked with DFT finite differences |
| Glass vibrational DOS cutoff | Same, for amorphous hosts | MACE melt-quench MD on GPU |

### Dopant in host
| Output | Why it matters | Method |
|---|---|---|
| 4f^N energy levels (2S+1)L_J and Stark splits | Which transitions exist | Determinant CI: Slater F^k, ζ_so, α,β,γ + crystal-field B^k_q |
| Reduced matrix elements U^(2,4,6) | Intensities | From CI eigenvectors |
| Judd–Ofelt Ω_2,4,6 | Host-specific intensity | Literature (validation) → later ab initio |
| Radiative rates A, branching β, τ_rad | Gain strength | Judd–Ofelt + magnetic dipole |
| Nonradiative rate W_nr, quantum efficiency η | Whether the level actually emits | Energy-gap law using computed ħω_max |
| Emission cross section σ_e(λ) | Gain coefficient | Füchtbauer–Ladenburg (later: McCumber) |
| 4f→5d / 3d→3d absorption + emission, Stokes shift | Broadband / UV gain centers | Point-charge-embedded cluster CASSCF/NEVPT2 |
| 5d position vs conduction band | Photoionization / ESA loss in UV | Combine cluster 5d energy with host E_g |

### Derived amplifier metrics
gain window (λ_min, λ_max), σ_e·τ, η, ground-state-absorption overlap.

## 3. Validation set

Hosts: α-SiO2 (silica proxy), CaF2, LiYF4, Al2O3, Y3Al5O12, ZBLAN glass, SiO2 glass.
Centers: Nd3+:YAG (1064 nm), Er3+:silica (1.55 µm), Er3+:ZBLAN (2.7 µm), Ce3+:LiYF4 (UV 5d), Ti3+:Al2O3 (broadband NIR).

Each computed number is compared with an experimental value whose source is cited in
`gainmat/reference_data.py`.

## 4. Where ML goes later
Tier-1 host outputs are cheap → screen thousands of compositions.
Tier-2 cluster/CI outputs are expensive → active learning picks which to run.
Multi-fidelity model: PBE → GLLB-SC/MLIP → CASSCF/NEVPT2 → experiment.
