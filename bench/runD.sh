#!/bin/bash
# Pass 2c: eps_inf (RPA+LFE, +/- scissor to GLLB-SC gap) for CaF2/Al2O3/quartz, LiF nbands check; then quartz TB09 retry.
cd "$(dirname "$0")"
NP=${NP:-$(nproc)}
MPI="mpiexec -n $NP --allow-run-as-root --oversubscribe"
export OMP_NUM_THREADS=1 GPAW_MPI_BACKEND=cgpaw
./runC.sh
cd ../eps
# scissor = GLLBSC - PBE gap (first pass, ecut 400)
$MPI python bench_eps.py CaF2 6 6 60 5.69 > log_CaF2 2>&1
$MPI python bench_eps.py LiF 8 8 60 5.78 > log_LiF_nb8 2>&1
$MPI python bench_eps.py Al2O3 4 4 50 3.42 > log_Al2O3 2>&1
$MPI python bench_eps.py SiO2-alpha-quartz 4 4 50 3.95 > log_quartz 2>&1
cd ../bench
MIXBETA=0.05 $MPI python bench_gap.py SiO2-alpha-quartz TB09 > log_SiO2-alpha-quartz_TB09 2>&1
echo DONE > doneD
