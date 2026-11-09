"""Time a few SCF iterations of YAG (80 atoms) PBE and GLLB-SC at the 'full' settings
(ecut 500, Gamma-centred 2x2x2) to extrapolate the full run time.
Usage: GPAW_MPI_BACKEND=cgpaw mpiexec -n N python time_yag.py XC NITER"""
import os, sys, time, json, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from gpaw import GPAW, PW, FermiDirac
from gpaw.mpi import world
from gainmat.structures import get_structure
xc, nit = sys.argv[1], int(sys.argv[2])
a = get_structure("Y3Al5O12")
d = os.path.dirname(os.path.abspath(__file__))
calc = GPAW(mode=PW(500), xc=xc, kpts={"size": (2, 2, 2), "gamma": True},
            occupations=FermiDirac(0.01), convergence={"bands": "occupied"},
            maxiter=nit, txt=f"{d}/yag_{xc}_np{world.size}.txt")
a.calc = calc
t0 = time.time()
try:
    a.get_potential_energy()
except Exception as e:
    pass
t = time.time() - t0
nk = len(calc.get_ibz_k_points()); nb = calc.get_number_of_bands()
if world.rank == 0:
    out = dict(xc=xc, niter=nit, nproc=world.size, t_s=t, nibzk=nk, nbands=nb, natoms=len(a))
    print(json.dumps(out))
    open(f"{d}/yag_timing.jsonl", "a").write(json.dumps(out) + "\n")
