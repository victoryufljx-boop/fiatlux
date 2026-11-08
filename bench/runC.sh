#!/bin/bash
# Pass 2b: GLLB-SC convergence (ecut, k-mesh) on LiF and CaF2, under MPI.
cd "$(dirname "$0")"
NP=${NP:-$(nproc)}
MPI="mpiexec -n $NP --allow-run-as-root --oversubscribe"
export OMP_NUM_THREADS=1 GPAW_MPI_BACKEND=cgpaw
for m in LiF CaF2; do
  for e in 400 500 600 700; do ECUT=$e KMESH=6,6,6 $MPI python bench_gap.py $m GLLBSC > log_conv_${m}_e$e 2>&1; done
  for k in 4 8 10; do ECUT=500 KMESH=$k,$k,$k $MPI python bench_gap.py $m GLLBSC > log_conv_${m}_k$k 2>&1; done
done
echo CONVDONE > doneC
