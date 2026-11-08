"""Gap bake-off: PBE, GLLB-SC, TB09 (mBJ), HSE06@PBE (nsc). Usage: python bench_gap.py MAT METHOD"""
import os, sys, time, json, warnings
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
warnings.filterwarnings("ignore")
import numpy as np
from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac, Mixer

def struct(m):
    if m == "LiF": return bulk("LiF", "rocksalt", a=4.0263)
    if m == "MgO": return bulk("MgO", "rocksalt", a=4.212)
    from gainmat.structures import get_structure
    return get_structure(m)

KP = {"LiF": (6,6,6), "MgO": (6,6,6), "CaF2": (6,6,6), "Al2O3": (4,4,4), "SiO2-alpha-quartz": (4,4,4)}
ECUT = float(os.environ.get("ECUT", 400))
mat, meth = sys.argv[1], sys.argv[2]
atoms = struct(mat)
kp = {"size": KP[mat], "gamma": True}
out = dict(mat=mat, method=meth, ecut=ECUT, kpts=KP[mat], natoms=len(atoms))
txt = f"/home/claude/work/gap/bench/{mat}_{meth}.txt"
t0 = time.time()
if meth in ("PBE", "HSE06"):
    calc = GPAW(mode=PW(ECUT), xc="PBE", kpts=kp, occupations=FermiDirac(0.01),
                nbands=-8 if meth=="PBE" else None, convergence={"bands": "occupied"}, txt=txt)
    atoms.calc = calc; atoms.get_potential_energy()
    h, l = calc.get_homo_lumo(); out["gap_pbe"] = float(l - h); out["t_scf"] = time.time()-t0
    if meth == "HSE06":
        from gpaw.hybrids.eigenvalues import non_self_consistent_eigenvalues as nsc
        eig = calc.get_eigenvalues  # find VBM/CBM k
        nk = len(calc.get_ibz_k_points()); nel = int(round(calc.get_number_of_electrons()))
        nv = nel // 2
        e = np.array([calc.get_eigenvalues(kpt=k) for k in range(nk)])
        kv = int(np.argmax(e[:, nv-1])); kc = int(np.argmin(e[:, nv]))
        ks = sorted({kv, kc, 0})
        t1 = time.time()
        edft, vdft, vhyb = nsc(calc, "HSE06", n1=nv-1, n2=nv+1, kpt_indices=ks)
        ehyb = edft - vdft + vhyb
        vb = ehyb[0, :, 0]; cb = ehyb[0, :, 1]
        out["kpts_used"] = ks; out["kv"] = kv; out["kc"] = kc
        out["gap_hse"] = float(cb.min() - vb.max())
        out["gap_hse_direct_G"] = float(cb[ks.index(0)] - vb[ks.index(0)])
        out["t_hse_nsc"] = time.time() - t1
elif meth == "GLLBSC":
    calc = GPAW(mode=PW(ECUT), xc="GLLBSC", kpts=kp, occupations=FermiDirac(0.01),
                nbands=-8, convergence={"bands": "occupied"}, txt=txt)
    atoms.calc = calc; atoms.get_potential_energy()
    h, l = calc.get_homo_lumo()
    r = calc.hamiltonian.xc.response
    ks_gap, dxc = r.calculate_discontinuity(r.calculate_discontinuity_potential(h, l))
    out.update(gap_ks=float(ks_gap), dxc=float(dxc), gap_gllbsc=float(ks_gap + dxc))
elif meth == "TB09":
    calc = GPAW(mode=PW(ECUT), xc={"name": "TB09", "stencil": 1}, kpts=kp,
                occupations=FermiDirac(0.01), mixer=Mixer(0.5, 5, 50.0), nbands=-8,
                convergence={"bands": "occupied"}, txt=txt)
    atoms.calc = calc; atoms.get_potential_energy()
    h, l = calc.get_homo_lumo(); out["gap_tb09"] = float(l - h)
    try: out["c_tb09"] = float(calc.hamiltonian.xc.c)
    except Exception: pass
out["t_total"] = time.time() - t0
print(json.dumps(out))
with open("/home/claude/work/gap/bench/results.jsonl", "a") as f: f.write(json.dumps(out) + "\n")
