import json,sys
an=json.load(open(sys.argv[1])); timing=json.load(open(sys.argv[2])) if len(sys.argv)>2 else {}
src_frontiers="R. S. Welch, C. J. Wilkinson, J. C. Mauro, C. B. Bragatto, Front. Mater. 6:121 (2019), Table 1 ('from Pedone et al., 2006'), https://www.frontiersin.org/journals/materials/articles/10.3389/fmats.2019.00121/full"
src_gong="Y. Gong, E. Olivetti, arXiv:2209.02046, Table 2 ('force field parameters used for the Pedone force field'), https://arxiv.org/pdf/2209.02046"
src_chen="Chen, He, Gong, arXiv:2605.11604, Supplementary Table S1 (N-A-S glass, 'Based on the empirical interatomic potential model of Pedone et al.'), https://arxiv.org/pdf/2605.11604"
src_amorph="amorphouspy (glasagent/amorphouspy, BAM), amorphouspy/src/amorphouspy/lammps/potentials/pmmcs_potential.py, dict pmmcs_potential_params, commit 37e6d2eb7206fe6de05e8cd0fbdc136d7cf611d8, https://github.com/glasagent/amorphouspy/blob/37e6d2eb7206fe6de05e8cd0fbdc136d7cf611d8/amorphouspy/src/amorphouspy/lammps/potentials/pmmcs_potential.py (docstring: 'PMMCS (Pedone) parameters'; docs cite Pedone et al. 2006)"
pot=dict(
 name="Pedone, Malavasi, Menziani, Cormack, Segre, J. Phys. Chem. B 110, 11780 (2006) (PMMCS), rigid-ion partial-charge",
 functional_form="U_ij(r) = q_i q_j e^2/(4 pi eps0 r) + D_ij[(1-exp(-a_ij(r-r0_ij)))^2 - 1] + C_ij/r^12  (LAMMPS pair_style pedone + coul/long)",
 primary_paper_note="The 2006 paper itself (DOI 10.1021/jp0611018) could not be read in this session (ACS page and the only free PDF mirror were blocked). Every number below was taken from a secondary source that reproduces the Pedone table; Si-O, Al-O, O-O are cross-checked in 2-3 independent sources; Er-O comes from a single source (code), whose Si/Al/O entries agree digit-for-digit with the cross-checked literature values.",
 charges=dict(O=dict(q=-1.2,source=[src_gong,src_chen]),Si=dict(q=2.4,source=[src_gong,src_chen]),Al=dict(q=1.8,source=[src_gong,src_chen]),Er=dict(q=1.8,source=[src_amorph])),
 pairs={
  "O-O":dict(D_eV=0.042395,a_invA=1.379316,r0_A=3.618701,C_eVA12=22.0,source=[src_frontiers,src_chen,src_gong+" (rounded: 1.37932, 3.6187)",src_amorph]),
  "Si-O":dict(D_eV=0.340554,a_invA=2.006700,r0_A=2.100000,C_eVA12=1.0,source=[src_frontiers,src_chen,src_gong,src_amorph]),
  "Al-O":dict(D_eV=0.361581,a_invA=1.900442,r0_A=2.164818,C_eVA12=0.9,source=[src_chen,src_gong+" (rounded: 1.90044, 2.16482)",src_amorph]),
  "Er-O":dict(D_eV=0.040448,a_invA=2.294078,r0_A=2.837722,C_eVA12=3.0,source=[src_amorph],note="single (code) source; not verified against Table 1 of the 2006 paper — please check with journal access"),
 },
 not_used=dict(**{"Ge-O":"no Ge in this glass; amorphouspy lists (0.158118, 2.294230, 2.261313, C=5.0) - not used, not verified"}),
 cation_cation="Coulomb only (as in the Pedone model)",
 cutoffs=dict(short_range_A=7.5,coulomb_real_space_A=10.0,kspace="PPPM, relative force accuracy 1e-4",note="7.5 A short-range cutoff as in Gong & Olivetti Table 2"),
 code="LAMMPS 22 Jul 2025 (PyPI wheel lammps-2025.7.22.4.0)")
comp=dict(SiO2=870,Al2O3=86,Er2O3=4,mol_percent=dict(SiO2=round(870/960*100,3),Al2O3=round(86/960*100,3),Er2O3=round(4/960*100,3)),
          atoms=dict(Si=870,Al=172,Er=8,O=2010),n_atoms=3060)
prot=dict(ensemble="NVT, fixed experimental density 2.28 g/cm3, cubic box L=35.719 A",n_atoms=3060,density_g_cm3=2.28,
 initial="random placement, min distance 1.7 A, per-replica seed; CG minimisation (<=2000 iter)",
 thermostat="Nose-Hoover (LAMMPS fix nvt), Tdamp 0.1 ps",timestep_fs=1.0,
 stages=[dict(T_K=6000,ps=20),dict(T_K=3500,ps=50),dict(T_K="3500->300 linear",ps=213.333,rate_K_per_ps=15),dict(T_K=300,ps=30)],
 quench_rate_K_per_ps=15,quench_rate_note="Planned 10 K/ps (5 K/ps unaffordable: only 2 cores, ~70 steps/s per replica). The first 10 K/ps quench of both replicas was killed at 370 ps (T~490 K) by a sandbox restart with no mid-quench restart file, so both replicas were re-quenched from their saved end-of-3500-K states at 15 K/ps with a 30 ps (not 50 ps) 300 K stage to fit the budget.",
 replicas=2,replica_seeds={"1":dict(placement_seed=1,velocity_seed=7930),"2":dict(placement_seed=2,velocity_seed=15849)},
 sampling="10 frames, every 2 ps over the last 20 ps (t=12..30 ps of the 300 K stage)",timing=timing)
summ=an['summary']
out=dict(glass="Er-Al9-SiO2",description="Er-doped aluminosilicate, 9.0 mol% Al2O3 (Simpson thesis fibre Al-7); classical MD melt-quench",
 composition=comp,potential=pot,protocol=prot,er_sites=an['er_sites'],summary=summ,
 ligand_tag_legend="BO-X-Y: bridging O bonded to formers X,Y (Si within 2.0 A, Al within 2.1 A); NBO-X: non-bridging, one former; free_O: no former; TBO: >=3 formers (tricluster). ligands are O within r_min of the Er, minimum image, A, sorted by distance")
json.dump(out,open(sys.argv[3],'w'),separators=(',',':'))
