cd /home/claude/work/gap/bench
for m in CaF2 Al2O3 SiO2-alpha-quartz; do for meth in PBE GLLBSC HSE06; do
 OMP_NUM_THREADS=1 python bench_gap.py $m $meth > log_${m}_${meth} 2>&1; done; done
for m in CaF2 Al2O3 SiO2-alpha-quartz; do OMP_NUM_THREADS=1 python bench_gap.py $m TB09 > log_${m}_TB09 2>&1; done
echo DONE > doneA
