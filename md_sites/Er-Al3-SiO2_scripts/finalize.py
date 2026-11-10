import json, glob, os, numpy as np, subprocess
from ase.io import read, write
A=json.load(open('analysis.json')); S=A['summary']
reps=S['replicas']
timing={}
for r in reps:
    t0=int(open(f'rep{r}.start').read()); t1=int(open(f'rep{r}.end').read()); timing[r]=t1-t0
# per-replica final extxyz
os.makedirs('/home/claude/fiatlux/md_sites',exist_ok=True)
for r in reps:
    out=f'/home/claude/fiatlux/md_sites/Er-Al3-SiO2_rep{r}_300K.extxyz'
    if not os.path.exists(out):
        a=read(f'rep{r}_300K.data',format='lammps-data',style='charge',Z_of_type={1:14,2:13,3:68,4:8})
        write(out,a,format='extxyz')
pot={'name':'Pedone et al. (PMMCS) partial-charge rigid-ion potential',
 'reference':'A. Pedone, G. Malavasi, M. C. Menziani, A. N. Cormack, U. Segre, J. Phys. Chem. B 110, 11780-11795 (2006), doi:10.1021/jp0611018',
 'functional_form':'U(r)=q_i q_j e^2/r + D_ij[(1-exp(-a_ij(r-r0)))^2-1] + C_ij/r^12 (LAMMPS pair_style pedone)',
 'charges':{'O':-1.2,'Si':2.4,'Al':1.8,'Er':1.8},
 'pairs':{'Si-O':{'D_eV':0.340554,'a_invA':2.006700,'r0_A':2.100000,'C_eVA12':1.0},
          'Al-O':{'D_eV':0.361581,'a_invA':1.900442,'r0_A':2.164818,'C_eVA12':0.9},
          'Er-O':{'D_eV':0.040448,'a_invA':2.294078,'r0_A':2.837722,'C_eVA12':3.0},
          'O-O':{'D_eV':0.042395,'a_invA':1.379316,'r0_A':3.618701,'C_eVA12':22.0},
          'Ge-O':'not used (no Ge in this glass); amorphouspy lists D=0.158118 a=2.294230 r0=2.261313 C=5.0'},
 'cation_cation':'Coulomb only',
 'parameter_source':('All numbers (charges, D, a, r0, C for Si-O, Al-O, Er-O, O-O) were copied from the PMMCS parameter dictionary '
   'in amorphouspy v0.8.1 (PyPI wheel amorphouspy-0.8.1-py3-none-any.whl, file amorphouspy/lammps/potentials/pmmcs_potential.py, '
   'BAM / A. Atila), which states they are the Pedone 2006 (PMMCS) parameters. The original J. Phys. Chem. B table could not be '
   'fetched in this run (ACS paywall/blocked). Si-O, Al-O and O-O agree with the values widely reproduced in the literature; '
   'the Er-O row could NOT be independently cross-checked against the original paper or a second reproduction in this run.'),
 'electrostatics':'coul/dsf (damped shifted force), alpha=0.25 1/A, cutoff 8.0 A (amorphouspy PMMCS default); Morse+C/r^12 cutoff 5.5 A (Pedone 2006)',
}
proto={'ensemble':'NVT (fixed box), Nose-Hoover chain, Tdamp 0.1 ps','density_g_cm3':2.22,'n_atoms':3088,
 'counts':{'Si':966,'Al':68,'Er':8,'O':2046},'box_A':36.1271,'timestep_fs':1.0,
 'stages':[{'stage':'random placement (min dist O-O 2.2, cation-O 1.5, cation-cation 2.8 A) + CG minimisation','time_ps':0},
   {'stage':'Langevin 3000->6000 K with nve/limit 0.05 A/step (contact relaxation)','time_ps':3},
   {'stage':'NVT 6000 K','time_ps':20},{'stage':'NVT ramp 6000->3500 K','time_ps':2},{'stage':'NVT 3500 K','time_ps':50},
   {'stage':'NVT linear quench 3500->300 K','time_ps':320,'rate_K_per_ps':10.0},{'stage':'NVT 300 K','time_ps':50}],
 'quench_rate_K_per_ps':10.0,
 'quench_rate_note':'10 K/ps used instead of 5 K/ps: only 2 CPU cores; measured ~7 ms/step per single-core replica (2 replicas run concurrently)',
 'sampling':'10 frames every 2 ps over the last 20 ps of the 300 K stage',
 'replicas':len(reps),'replica_ids':reps,'replica_seeds':{r:{'placement':1000+r,'velocity':5000+7*r} for r in reps},
 'wall_time_s_per_replica':timing,'software':'LAMMPS 22 Jul 2025 (PyPI lammps 2025.7.22.4.0), serial per replica'}
glass={'name':'Er-Al3-SiO2','description':'alumino-silicate, 3.4 mol% Al2O3 (Simpson thesis fibre Al-2), Er as single RE species (stands in for Er/Yb/Tm)'}
comp={'formula_units':{'SiO2':966,'Al2O3':34,'Er2O3':4},'mol_percent':{'SiO2':round(966/1004*100,3),'Al2O3':round(34/1004*100,3),'Er2O3':round(4/1004*100,3)}}
out={'glass':glass,'composition':comp,'potential':pot,'protocol':proto,'summary':S,'er_sites':A['er_sites']}
json.dump(out,open('/home/claude/fiatlux/md_sites/Er-Al3-SiO2.json','w'),separators=(',',':'))
print('ok',os.path.getsize('/home/claude/fiatlux/md_sites/Er-Al3-SiO2.json'))
