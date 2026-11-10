"""Write md_sites/Er-SiO2.md from md_sites/Er-SiO2.json (numbers) + fixed narrative."""
import json, os
here = os.path.dirname(os.path.abspath(__file__))
J = json.load(open(f"{here}/../Er-SiO2.json"))
S, P = J["summary"], J["protocol"]
san = S["sanity"]
nrep = P["replicas"]
t = P["replica_timing"]
rows = "\n".join(f"| {r} | {v['seed']} | {v['wall_s_total']/60:.0f} | {v['T_mean_300K_stage']} | {'yes' if v.get('resumed_after_container_restart') else 'no'} |" for r, v in t.items())
cnv = ", ".join(f"{k} A: {v['CN_mean']:.2f} +- {v['CN_std']:.2f}" for k, v in S["CN_vs_cutoff"].items())
ch = ", ".join(f"CN{k}: {v}" for k, v in S["CN_hist"].items())
ee = S["er_er_nn"]
chk = J["dsf_vs_ewald_check"]
md = f"""# Er-SiO2: classical-MD Er sites (Pedone 2006 potential, slow quench to 300 K)

Branch `md2/Er-SiO2`. Data: `Er-SiO2.json`. Final structures: `Er-SiO2_<rep>_final.extxyz`. Scripts: `tools/`.

## What was run

* **Glass:** 1000 SiO2 + 4 Er2O3 (8 Er, 0.40 mol% Er2O3), 3020 atoms, cubic box 35.961 A, fixed density 2.20 g/cm3 (NVT).
* **Potential:** Pedone et al., J. Phys. Chem. B 110, 11780 (2006), Table 2. Partial charges O -1.2, Si +2.4, Er +1.8. Coulomb + Morse + C/r^12.
  * O-O: D 0.042395 eV, a 1.379316 1/A, r0 3.618701 A, C 22.0 eV A^12.
  * Si-O: D 0.340554, a 2.006700, r0 2.100000, C 1.0.
  * Er-O: D 0.040448, a 2.294078, r0 2.837722, C 3.0.
  * All values come from Table 2 of the paper (PDF copy at lptmc.jussieu.fr; see JSON). Er-O is in the original Pedone table, so no other parameter source was needed. Al-O and Ge-O values from the same table are listed in the JSON for the other glasses but are not used here.
* **Engine:** LAMMPS (PyPI wheel 2025.7.22.4.0), one serial process per replica, replicas run concurrently on the 2 available cores.
* **Protocol, per replica:**
  1. Random placement with a minimum-distance check, then CG minimisation.
  2. 6000 K for 20 ps, then 3500 K for 50 ps.
  3. Linear cool 3500 -> 300 K at **10 K/ps** (320 ps).
  4. 300 K for 50 ps.
  * Nose-Hoover NVT, Tdamp 0.1 ps, timestep 1 fs.
  * Sites come from the last 20 ps: 10 frames, 2 ps apart.
* **Replicas:** {nrep} independent random starts, {len(J['er_sites'])} site-frames, **{nrep*8} distinct Er ions**.

| replica | velocity seed | wall time (min) | mean T in 300 K stage (K) | resumed after container restart |
|---|---|---|---|---|
{rows}

## Deviations from the requested protocol (please read)

1. **Quench rate 10 K/ps, not 5 K/ps.** The machine had only 2 cores, and the speed I measured was about 120 steps/s per core. At 5 K/ps one replica would take about 1.8 h, so I used the 10 K/ps fallback.
2. **Coulomb is DSF, not Ewald/PPPM.** I used damped shifted force (alpha 0.2 1/A, cutoff 10 A), tabulated together with Morse + C/r^12 into one `pair_style table`. This was 3.5x faster than PPPM 1e-5 with a 10 A cutoff (35 -> 120 steps/s). The short-range cutoff is 10 A, against 15 A in the paper; the Morse tails are below 2e-5 eV at 10 A.
   * **Check on replica rep3:**
     * Single point: instantaneous forces differ noticeably from Ewald. The rms force difference is {chk['single_point']['rms_force_difference_eV_A']} eV/A, against rms forces of {chk['single_point']['rms_force_Ewald_eV_A']} eV/A.
     * Pressure: {chk['single_point']['pressure_bar_DSF']} bar (DSF) vs {chk['single_point']['pressure_bar_Ewald']} bar (Ewald).
     * PPPM relaxation: I re-equilibrated the final structure for 19 ps at 300 K with PPPM. The Er shell did not change. At a 3.0 A cutoff: CN {chk['Er_shell_at_fixed_cutoff_3.0A']['DSF']['CN_mean']} -> {chk['Er_shell_at_fixed_cutoff_3.0A']['PPPM_relaxed']['CN_mean']}, R {chk['Er_shell_at_fixed_cutoff_3.0A']['DSF']['R_mean']} -> {chk['Er_shell_at_fixed_cutoff_3.0A']['PPPM_relaxed']['R_mean']} A.
3. **Fewer replicas than hoped: {nrep} replicas, {nrep*8} distinct Er (target was >= 60).**
   * The cloud container restarted at about 01:25 UTC and killed the running jobs. Replica rep2 was resumed from its end-of-quench restart file. rep1 and rep4 had to restart their quench from the end of the 3500 K stage, which cost about 40 min of compute.
   * A planned fifth replica was dropped to stay inside the 3 h 15 min budget.
   * The 10 frames per Er are correlated snapshots of the same ion, not independent sites.
4. **Large residual tension.** At the fixed density of 2.20 g/cm3, the 300 K structure is under about -3 GPa (DSF and Ewald agree). The Pedone potential prefers a denser silica, so this follows from the requested NVT density, not from the Coulomb treatment.

## Sanity checks (pooled over replicas, last 20 ps at 300 K)

| quantity | result | target |
|---|---|---|
| Si-O first peak | {san['SiO_rdf_peak']:.2f} A | ~1.61 A |
| Si coordination | {san['Si_CN_mean']:.3f} ({100*san['Si_CN4_fraction']:.1f}% four-fold) | 4 |
| O-Si-O angle | mean {san['OSiO_mean']:.1f}, std {san['OSiO_std']:.1f} deg | ~109.5 deg |
| Si-O-Si angle | mean {san['SiOSi_mean']:.1f}, std {san['SiOSi_std']:.1f} deg, mode {san['SiOSi_mode']:.0f} deg | ~145-150 deg |
| Er-O RDF first peak | {S['ErO_rdf_peak']:.2f} A | EXAFS ~2.25 A |
| Er-O first minimum r_min | {S['r_min']:.2f} A | 2.9-3.2 A expected |
| Er CN (O within r_min) | {S['CN_mean']:.2f} +- {S['CN_std']:.2f} ({ch}) | ~6 (EXAFS) |
| mean Er-O within r_min | {S['R_mean']:.3f} +- {S['R_std']:.3f} A | ~2.25 A |
| NBO fraction in the first shell | {S['nbo_fraction']:.2f} | |
| O-Er-O mean angle | {S['O-Er-O']['mean']:.1f} deg | |

* **Er CN depends on the cutoff.** {cnv}.
* **Er-Er nearest neighbours:** mean {ee['mean']:.2f} A, minimum {ee['min']:.2f} A. In {100*ee['fraction_below_5A']:.0f}% of Er-frames the nearest Er is closer than 5 A, i.e. Er-O-Er pairs, a clustering signature. Per-Er values for the last frame are in the JSON.

### Assessment

* **Silica network: fine.** It is fully tetrahedral and O-Si-O is correct. The Si-O-Si angle is a few degrees above the experimental ~147 deg, a known trait of rigid-ion potentials at fast quench rates.
* **Er shell: under-coordinated.** This is the main discrepancy with EXAFS.
  * The mean Er-O distance within r_min ({S['R_mean']:.2f} A) is close to the EXAFS 2.25 A. The most probable distance (the RDF peak at {S['ErO_rdf_peak']:.2f} A) is somewhat short.
  * There is a sharp peak of about 3 O, mostly non-bridging, then a tail out to about 2.7 A. g(r) then stays low but non-zero to about 3.5 A, so the "first minimum" is shallow and depends on the cutoff (see the CN-vs-cutoff line).
  * Within r_min, Er has about 4-5 O instead of ~6, and many shells are strongly distorted.
  * This is what the Pedone Er-O parameters give in pure silica at this quench rate. The Coulomb treatment is not the cause (see the PPPM check above).
  * The downstream crystal-field step should treat these as **low-CN, strongly distorted Er sites**, probably more distorted than in real fused silica.
* **Extra neighbours stored.** Every site record also lists the O between r_min and 3.3 A (`outer_O_rmin_to_3.3A`), so a larger cutoff can be applied downstream without rerunning anything.

## JSON layout

* **Top-level keys:** `glass`, `composition`, `potential` (every parameter with its source), `protocol`, `dsf_vs_ewald_check`, `er_sites`, `summary`.
* **`er_sites[i]` fields:**
  * `replica`, `frame` (0-9), `timestep`, `er_index` (LAMMPS atom id).
  * `ligands`: `[element, x, y, z]` relative to the Er, in A, minimum image, sorted by distance.
  * `ligand_tags`: one entry per ligand. `tag` is BO, NBO, free_O or O_3c (Si within 2.0 A). It also gives `formers` (the bonded Si ids) and `atom_id`.
  * `outer_O_rmin_to_3.3A`: the extra O described above.
* **`summary` fields:** RDFs, angle histograms, CN histogram, Er-Er nearest-neighbour statistics.
"""
open(f"{here}/../Er-SiO2.md", "w").write(md)
print("ok")
