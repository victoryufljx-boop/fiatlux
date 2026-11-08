"""Calibrate: CaF2 'full'-level gap (PBE+GLLB-SC, ecut 500, path band edges) serial vs MPI launcher."""
import json, os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
os.environ["OMP_NUM_THREADS"] = "1"
from gainmat.host.run_host import LEVELS
from gainmat.structures import get_structure
from gainmat.host.mpi_run import run_gpaw_mpi
d = os.path.dirname(os.path.abspath(__file__))
mode = sys.argv[1]
if mode == "serial":
    from gainmat.host.gap import compute_gap
    t0 = time.time(); g = compute_gap(get_structure("CaF2"), **LEVELS["full"]["gap"]); t = time.time() - t0
    out = dict(mode="serial", nproc=1, t_s=t, gap=g)
else:
    n = int(sys.argv[2])
    t0 = time.time(); r = run_gpaw_mpi("CaF2", level="full", nproc=n, workdir=f"{d}/caf2_mpi"); t = time.time() - t0
    out = dict(mode="mpi", nproc=n, t_s=t, gap=r["gap"])
print(json.dumps(out, default=float))
open(f"{d}/calib_caf2.jsonl", "a").write(json.dumps(out, default=float) + "\n")
