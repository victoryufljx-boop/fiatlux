#!/bin/bash
# Pass 2: TB09 for the remaining hosts + GLLB-SC convergence, all under MPI.
cd "$(dirname "$0")"
NP=${NP:-$(nproc)}
MPI="mpiexec -n $NP --allow-run-as-root --oversubscribe"
export OMP_NUM_THREADS=1 GPAW_MPI_BACKEND=cgpaw
for m in CaF2 Al2O3 SiO2-alpha-quartz; do $MPI python bench_gap.py $m TB09 > log_${m}_TB09 2>&1; done
echo TB09DONE > doneB
