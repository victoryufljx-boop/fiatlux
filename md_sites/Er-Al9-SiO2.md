# Er-Al9-SiO2: classical-MD Er sites (Pedone 2006 potential, slow quench to 300 K)

Glass: Er-doped alumino-silicate, 9.0 mol% Al2O3 (Simpson thesis fibre Al-7).
Run 2026-11-09 23:40 to 2026-11-10 ~02:30 UTC on **2 CPU cores**. Data: `Er-Al9-SiO2.json`.
Final configurations: `Er-Al9-SiO2_r{1,2}_final.extxyz`. Scripts, LAMMPS logs and the pooled RDF: `md2_Er-Al9-SiO2/`.

## What was run

- **Composition.** 870 SiO2 + 86 Al2O3 + 4 Er2O3, which is 3060 atoms (Si 870, Al 172, Er 8, O 2010). That gives 8.958 mol% Al2O3 and 0.417 mol% Er2O3. The cell is charge-neutral.
- **Box.** Cubic, L = 35.719 Å, NVT at a fixed density of 2.28 g/cm3.
- **Potential.** Pedone et al., J. Phys. Chem. B 110, 11780 (2006): partial charges + Morse + C/r^12, using LAMMPS `pair_style pedone` (7.5 Å cutoff) overlaid with `coul/long` (10 Å real-space cutoff, PPPM accuracy 1e-4). Cation-cation pairs are Coulomb only. LAMMPS version: 22 Jul 2025 (PyPI wheel).

| pair | D (eV) | a (1/Å) | r0 (Å) | C (eV Å^12) | source |
|---|---|---|---|---|---|
| O-O | 0.042395 | 1.379316 | 3.618701 | 22.0 | Welch et al., Front. Mater. 6:121 (2019) Tab. 1; Chen et al., arXiv:2605.11604 Tab. S1; Gong & Olivetti, arXiv:2209.02046 Tab. 2 |
| Si-O | 0.340554 | 2.006700 | 2.100000 | 1.0 | same three |
| Al-O | 0.361581 | 1.900442 | 2.164818 | 0.9 | arXiv:2605.11604 Tab. S1; arXiv:2209.02046 Tab. 2 (rounded) |
| Er-O | 0.040448 | 2.294078 | 2.837722 | 3.0 | **only** amorphouspy `pmmcs_potential.py` (GitHub glasagent/amorphouspy @37e6d2e) |

Charges: O −1.2, Si +2.4, Al +1.8 (from the two arXiv tables) and Er +1.8 (from amorphouspy).

**Provenance caveat.** I could not open the 2006 paper itself in this session: the ACS page and the free PDF mirror were both blocked. Si-O, Al-O and O-O are cross-checked in 2–3 independent published reproductions. Er-O comes from a single source, a code repository that says it implements the Pedone (PMMCS) table. Its Si, Al and O entries match the published values digit for digit, but **someone with journal access should check the Er-O row against Table 1 of the paper.** No parameter was invented. Ge-O was not used because there is no Ge in this glass.

### Protocol (each replica)

| Step | Details |
|---|---|
| 1. Placement and minimisation | Random placement (minimum distance 1.7 Å), then CG minimisation |
| 2. High-temperature stage | 6000 K for 20 ps. This was stable: the Morse/r^-12 core held |
| 3. Melt | 3500 K for 50 ps |
| 4. Quench | Linear 3500 → 300 K at **15 K/ps** (213.3 ps) |
| 5. Anneal | 300 K for **30 ps** |

Throughout: timestep 1 fs, Nose-Hoover NVT with Tdamp = 0.1 ps. There are **2 independent replicas** (different placement and velocity seeds), so **16 distinct Er ions**. Sampling: 10 frames, one every 2 ps over the last 20 ps of the 300 K stage, giving 160 site snapshots.

Timings (s per replica, 1 core each, run concurrently):

| Stage | r1 (s) | r2 (s) |
|---|---|---|
| 6000 K, 20 ps | 306 | 304 |
| 3500 K, 50 ps | 722 | 724 |
| Quench, 213 ps | 2965 | 2975 |
| 300 K, 30 ps | 399 | 398 |

Speed was about 65–80 steps/s per replica.

## Deviations from the requested protocol

1. **Quench rate 15 K/ps instead of 5 K/ps, 300 K stage 30 ps instead of 50 ps.** The sandbox has only 2 cores, and 3060 atoms with PPPM ran at about 70 steps/s per replica. At that speed 5 K/ps was only possible for a single replica, so I planned 10 K/ps for two replicas. That first attempt (10 K/ps) was **killed by a sandbox restart at 370 ps of the quench (T ≈ 490 K)**. I had not written restart files during the quench, so it was lost. I restarted both replicas from their saved end-of-3500 K states and used 15 K/ps with a 30 ps anneal to finish inside the budget. Restart files were written every 10 ps after that.
2. **Only 2 replicas, so 16 distinct Er instead of ≥ 60.** Each Er appears in 10 thermal snapshots, but snapshots of the same ion are correlated. Treat the effective sample as 16 sites with thermal jitter.
3. **PPPM accuracy 1e-4** instead of a tighter value, and a 10 Å real-space cutoff, both for speed.

## Sanity checks (pooled over the 20 sampled frames)

| quantity | result | target |
|---|---|---|
| Si-O first peak | 1.61 Å | ~1.61 |
| Si coordination | 99.6 % 4-fold | 4 |
| O-Si-O | 109.4 ± 6.6° | ~109 |
| Si-O-Si | 154.0 ± 12.1° (mode ~157°) | 145–150 (a little wide/open) |
| Si-O-Al | 150.2 ± 13.7° | — |
| Al-O first peak | 1.71 Å | ~1.75 |
| Al coordination (cut 2.35 Å = g(r) min) | 85 % 4, 12 % 3, 3 % 5 | mostly 4 |
| O-Al-O | 109.8 ± 12.1° | — |
| O-O peak | 2.63 Å | ~2.62 |
| Bulk O species | 93.1 % BO, 0.8 % NBO, 6.2 % triclusters (mostly Al-rich) | — |
| **Er-O RDF peak** | **2.19 Å** | ~2.2–2.35 |
| **Er-O mean distance (within r_min)** | **2.34 ± 0.19 Å** (median 2.29) | EXAFS ~2.25 |
| **r_min (first minimum of Er-O g(r))** | **2.93 Å** | 2.9–3.2 |
| **Er CN** | **5.2 ± 0.9** (hist: 4:40, 5:59, 6:46, 7:15 of 160) | ~6–7 (EXAFS ~6) |

Other results:

- **NBO fraction in the Er first shell: 36 %.** Bulk NBO is only 0.8 %, so Er gathers NBOs around itself strongly.
- **Er sits near Al.** 33 % of Er-O-former linkages are Er-O-Al, while Al-O bonds are only 16 % of all O-former bonds in the bulk. That is a ~2× enrichment of Al in the second shell.
- **Er ligand tags** (counts over 160 snapshots): NBO-Si 291, BO-Al-Si 262, BO-Si-Si 188, BO-Al-Al 62, TBO (Al triclusters) 23, NBO-Al 10.
- **Er-Er nearest-neighbour distance:** 11.0 Å mean, 6.8 Å minimum. There is no Er clustering (none closer than 5 Å) at this low concentration in a 3060-atom box.

## Assessment and caveats

- These are fully quenched 300 K glasses. Much better than the earlier `forward/md-sites` attempt: Er-O is 2.34 ± 0.19 Å here versus 2.49 ± 0.40 Å there, and production now ends at 300 K instead of > 1000 K.
- **The Er CN (5.2) is lower than EXAFS (~6).** The Er-O well in this parameter set is very shallow (D = 0.04 eV), so coordination is set largely by the Coulomb term. Er-O lengths and the O-Er-O angle spread (103 ± 29°) are realistic, but low-CN sites (CN 4: 25 % of snapshots) may be over-represented. If the downstream crystal-field fit wants CN ≈ 6, consider weighting toward the CN 6–7 sites, or re-running with a verified Er-O row.
- The pressure at 300 K is strongly negative (−3.6 to −3.8 GPa at the end of the 300 K stage). At this composition the Pedone model prefers a density above the experimental 2.28 g/cm3 used here (fixed NVT, as requested). This, plus the fast quench, likely explains Si-O-Si being slightly too open.
- 12 % 3-fold Al plus O triclusters are expected for Al2O3-SiO2 with no network modifier: the only charge compensator for AlO4⁻ is the small amount of Er³⁺.
- The fastest affordable quench rate here (15 K/ps) is still ~10¹³× faster than in a real fibre draw. All classical MD shares this limitation.

## JSON layout

- `er_sites[]` has fields `{replica, frame, timestep, er_index, ligands: [["O", x, y, z], ...], ligand_tags, ligand_dist}`.
  - Ligands are the O atoms within r_min = 2.93 Å of the Er, given relative to the Er (minimum image, Å) and sorted by distance.
  - Tags: `BO-X-Y` = bridging O (Si within 2.0 Å, Al within 2.1 Å); `NBO-X` = one network former; `TBO-...` = oxygen bonded to 3 formers.
- `summary` holds the CN histogram, R mean/std, r_min, NBO fraction, Er-O-Al/Si linkages, Er-Er nearest neighbours, and the sanity checks including 5° angle histograms.
- `potential` gives the source for every number.
- `protocol` gives the stages, seeds and timings.
