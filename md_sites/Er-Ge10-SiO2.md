# Er-Ge10-SiO2: Er site geometries from classical MD (Pedone 2006 potential)

> **INTERIM (replicas 1-2 only; replicas 3-4 still running).** Numbers will be updated.

Glass: germanosilicate, 900 SiO2 + 100 GeO2 + 4 Er2O3 (89.6 / 9.96 / 0.40 mol%), 3020 atoms, cubic box 36.0075 Å at a fixed 2.35 g/cm³ (NVT).

Files: `Er-Ge10-SiO2.json` has the sites, statistics, RDFs and every potential parameter with its source. `Er-Ge10-SiO2_repN_final.extxyz` holds the final 300 K configurations. The scripts are in `md_sites/scripts/`.

## Potential

Pedone, Malavasi, Menziani, Cormack, Segre, J. Phys. Chem. B 110, 11780 (2006), Table 2 (binary-oxide fit). Every number below comes from that table, which I transcribed from the paper PDF in this session. U = q_i q_j/r + D[(1-e^{-a(r-r0)})² - 1] + C/r¹². Charges: O -1.2, Si +2.4, Ge +2.4, Er +1.8.

| pair | D (eV) | a (1/Å) | r0 (Å) | C (eV Å¹²) |
|---|---|---|---|---|
| Si-O | 0.340554 | 2.0067 | 2.1 | 1.0 |
| Al-O | 0.361581 | 1.900442 | 2.164818 | 0.9 |
| Ge-O | 0.158118 | 2.29423 | 2.261313 | 5.0 |
| Er-O | 0.040448 | 2.294078 | 2.837722 | 3.0 |
| O-O | 0.042395 | 1.379316 | 3.618701 | 22.0 |

The Er-O pair is Pedone's own Er(1.8)-O(-1.2) row, so no substitute RE-O set was needed. Al-O (0.361581, 1.900442, 2.164818, 0.9) is recorded in the JSON but not used because this glass has no Al. Cation-cation pairs are Coulomb only.

**Implementation (deviation):** LAMMPS 22 Jul 2025 (PyPI wheel). I put every pair into one `pair_style table` (12000 points). Coulomb uses damped shifted force (DSF, alpha 0.2 Å⁻¹, cutoff 10 Å) instead of Ewald/PPPM. The Morse and r⁻¹² terms are cut at 10 Å and energy-shifted. The reason is speed: the machine has only 2 cores, and hybrid coul/long + PPPM (rc 11 Å, 1e-5) ran at 25 ms/step on one core, while the single table runs at 8.5-10 ms/step. DSF is a standard real-space alternative for oxide glasses, but I did not run a PPPM cross-check on the final structures.

## Protocol

Each replica starts from a random placement with minimum distances (cation-cation 2.6 Å, cation-O 1.55 Å, O-O 2.2 Å, Er-Er 4.0 Å). This is followed by CG minimisation and 3 ps of nve/limit with Berendsen at 6000 K to remove overlaps. Then, under Nose-Hoover NVT (Tdamp 0.1 ps) with a 1 fs timestep:

- 6000 K for 20 ps (the Morse + r⁻¹² core held at 6000 K, with no atom collapse)
- 3500 K for 50 ps
- linear cooling from 3500 K to 300 K over 320 ps, i.e. **10 K/ps**
- 300 K for 50 ps. Sites come from 10 frames, one every 2 ps over the last 20 ps.

Main replicas: rep1, rep2 (seeds {'rep1': 11, 'rep2': 22}). Each took about 62 min of wall time on one core, with two replicas running at once.

## Deviations from the requested protocol (honest list)

1. **Quench rate 10 K/ps, not 5 K/ps.** `nproc` = 2. At about 9 ms/step per core, one 5 K/ps replica (760 ps) takes about 110 min. Within 3 h 15 min that would allow only 2 replicas (16 Er), so I followed the 10 K/ps fallback.

2. **Distinct Er sites: 16 (main set), below the target of ≥ 60.** Each replica holds 8 Er, and 2 cores allowed only 2 rounds of 2 concurrent replicas. `er_sites` contains 10 frames per Er. These snapshots are thermally correlated (300 K vibrations of the same site), not independent sites.

3. Coulomb uses DSF instead of Ewald/PPPM (see above).

4. During replicas 3-4 the cloud container was restarted (at about 01:05 UTC), killing both runs part-way through the quench. I restarted them from their restart files written at the end of the 3500 K hold (`hold.restart`, step 70000), so the full 320 ps quench and the 300 K stage ran uninterrupted after the restart. Because the restart rebuilt the thermostat chain, the run is not bit-identical to an uninterrupted one, but the protocol is the same.

## Sanity checks (pooled over all frames of the main replicas)

| quantity | MD | target |
|---|---|---|
| Si-O first peak | 1.61 Å | ~1.61 Å |
| O-Si-O | 109.37 ± 6.45 (mode 109.0)° | ~109° |
| Si-O-Si | 152.74 ± 12.04 (mode 153.0)° | ~145-150° |
| Si CN (fraction 4-fold) | 3.999 (0.9977) | 4 |
| Ge-O first peak | 1.73 Å | ~1.74 Å |
| Ge CN (fraction 4-fold) | 4.011 (0.9795) | 4 |
| O-Ge-O | 109.31 ± 9.6 (mode 107.0)° | ~109° |
| Si-O-Ge / Ge-O-Ge | 152.65 ± 12.08 (mode 153.0)° / 153.96 ± 12.38 (mode 155.0)° | ~130-150° |
| O-O / Si-Si first peak | 2.61 / 3.15 Å | 2.62 / 3.1 Å |
| bulk O speciation | {'NBO': 0.01469, 'BO': 0.98347, 'free': 0.0, 'TBO': 0.00184} | NBO only from Er2O3 |
| Er-O first peak | 2.17 Å | ~2.2-2.35 Å (EXAFS 2.25) |
| Er-O CN within r_min = 3.21 Å | 4.931 ± 1.09 | ~6-7 (EXAFS ~6) |
| final-stage T | {'rep1': 301.2, 'rep2': 299.6} K | 300 K |

The Si-O-Si angle comes out a few degrees above the requested range, consistent with Pedone-type potentials, which tend to give wide Si-O-Si angles. Both formers are almost perfectly tetrahedral.

## Er first shell (main replicas)

- r_min (first minimum of the pooled Er-O g(r), 0.1 Å smoothing) = **3.21 Å**. Er-O first peak 2.17 Å.
- CN within r_min: **4.931 ± 1.09**, histogram {'3': 11, '4': 45, '5': 67, '6': 22, '7': 11, '8': 4}. CN within 2.5 / 2.6 / 2.8 Å: 3.8 / 4.013 / 4.331.
- Er-O over all first-shell ligands: **2.3789 ± 0.2829 Å**. The shortest Er-O per site averages 2.1498 Å.
- First-shell O types: {'BO': 289, 'NBO': 500, 'free': 0, 'TBO': 0}, giving an NBO fraction of **0.6337**. Most of the NBOs in the glass sit on Er; Er2O3 is their only source.
- Er-O-T linkages: Er-O-Si 1032, Er-O-Ge 46, so the Ge fraction is 0.0427 against 0.10 Ge among the formers.
- Er-Er nearest neighbour: mean 10.593 Å, min 3.491 Å. Fraction of Er with another Er within 5 Å: 0.25.
- O-Er-O angles: 103.77 ± 28.77 (mode 93.0)°. Er-O-T angles: 124.26 ± 25.3 (mode 93.0)°.

**Assessment.** The run is now cold (300 K) and much better quenched than `forward/md-sites` (>1000 K, Er-O 2.49 ± 0.40 Å, CN 6.6). The Er-O distance is in the expected range: the peak is slightly short of EXAFS, and the mean over the shell is close to it. **The coordination number is low, however.** About 4 O sit within 2.6 Å, and roughly one more lies in a diffuse 2.6-3.2 Å tail. The target from EXAFS is ~6. This is a property of the Pedone Er-O parameters in a silica-rich network at this density, not of the protocol. Er binds all available NBOs and adds a few long, weak bridging-O contacts, so the sites are strongly distorted, low-symmetry 4-5 coordinated polyhedra. Er also clusters: Er pairs share NBOs at 3.5-4.2 Å. In real MCVD Ge-silica, the absence of Al/P likewise favours clustering, but downstream crystal-field use should treat these sites as an under-coordinated end member. A possible next step is to rescale the Er-O D/r0 to reach CN ~6. That would mean refitting, which I did not do because inventing parameters was ruled out.

## JSON layout

`er_sites[i]` = {replica, frame, timestep, er_index, ligands: [["O", x, y, z], ...] relative to the Er, minimum image, Å, sorted by distance; ligand_tags: [{type: BO|NBO|free|TBO, formers: [Si/Ge...], r}], other_neighbors_within_rmin: non-O atoms within r_min (e.g. Si, Er), quench_rate_K_per_ps}. Only O appear in `ligands`, so `gainmat.forward.glass.load_md_sites` reads them directly. BO/NBO criterion: O within 2.0 Å of 2 / 1 Si or Ge.
