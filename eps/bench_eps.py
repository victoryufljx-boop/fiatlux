"""RPA eps_inf (LFE) with/without scissor. Usage: python bench_eps.py MAT KMESH NBFACTOR ECUTRESP SCISSOR1 [SCISSOR2 ...]"""
import os, sys, time, json, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "/home/claude/fiatlux")
import numpy as np
from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac
from gpaw.response.df import DielectricFunction

def struct(m):
    if m == "LiF": return bulk("LiF", "rocksalt", a=4.0263)
    if m == "MgO": return bulk("MgO", "rocksalt", a=4.212)
    from gainmat.structures import get_structure
    return get_structure(m)

mat = sys.argv[1]; k = int(sys.argv[2]); nbf = float(sys.argv[3]); ecr = float(sys.argv[4])
sciss = [float(x) for x in sys.argv[5:]]
kmesh = (k, k, k) if mat not in ("Al2O3", "SiO2-alpha-quartz") else (k, k, k)
d = "/home/claude/work/gap/eps"
atoms = struct(mat)
t0 = time.time()
calc = GPAW(mode=PW(400), xc="PBE", kpts={"size": kmesh, "gamma": True},
            occupations=FermiDirac(0.01), txt=f"{d}/{mat}_gs.txt")
atoms.calc = calc; atoms.get_potential_energy()
nocc = int(round(calc.get_number_of_electrons())) // 2
nb = int(nbf * nocc)
calc = calc.fixed_density(nbands=nb + 8, convergence={"bands": nb}, txt=f"{d}/{mat}_nscf.txt")
gpw = f"{d}/{mat}_k{k}.gpw"; calc.write(gpw, mode="all")
tgs = time.time() - t0
out = dict(mat=mat, kpts=list(kmesh), nbands=nb, nocc=nocc, ecut_resp=ecr, t_gs=tgs)
dirs = ["x"] if mat in ("LiF", "MgO", "CaF2") else ["x", "z"]
for s in [0.0] + sciss:
    for dr in dirs:
        t1 = time.time()
        df = DielectricFunction(gpw, frequencies=[0.0], ecut=ecr, nbands=nb, intraband=False,
                                hilbert=False, eta=1e-4, eshift=(s or None), txt=f"{d}/{mat}_df.txt")
        e_nlf, e_lf = df.get_macroscopic_dielectric_constant(direction=dr)
        out[f"eps_nlf_s{s:.2f}_{dr}"] = float(np.real(e_nlf))
        out[f"eps_lf_s{s:.2f}_{dr}"] = float(np.real(e_lf))
        out[f"t_df_s{s:.2f}_{dr}"] = time.time() - t1
        print(mat, s, dr, e_nlf, e_lf, flush=True)
out["t_total"] = time.time() - t0
print(json.dumps(out))
with open(f"{d}/results_eps.jsonl", "a") as f: f.write(json.dumps(out) + "\n")
