import numpy as np, time
from pyscf import lib; lib.num_threads(1)
from gainmat.cluster import embedding as emb, relax
from gainmat.cluster.casci_lf import run_scf
t=time.time()
cl = emb.make_cluster("LiYF4","Y","Ce",3,r_qm=2.6,r_cap=4.0,r_pc=12.0,cap_mode="aimp")
print(cl.summary()['notes'], cl.total_charge(), len(cl.aimp_symbols))
mol,mf = run_scf(cl, {"Ce":"sbkjc","F":"def2-svp"},{"Ce":"sbkjc"})
print('AIMP ROHF', mf.e_tot, mf.converged, time.time()-t)
eg,ee,inf = relax.scf_energies(cl, {"Ce":"sbkjc","F":"def2-svp"},{"Ce":"sbkjc"})
print('dSCF', inf, time.time()-t)
