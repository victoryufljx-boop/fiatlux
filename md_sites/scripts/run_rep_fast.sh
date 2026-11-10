#!/bin/bash
# usage: run_rep.sh REP SEED
cd /home/claude/md
d=rep$1; mkdir -p $d
python3 build.py $2 $d > $d/build.out
python3 dsftab.py $d/all.table
cd $d; cp ../in.glass .
date +%s > t_start
lmp -in in.glass -var seed $2 -var tmelt 6000 -var nmelt 20000 -var nhold 20000 -var ncool 160000 -var nfinal 30000 -log log.lammps -screen none
date +%s > t_end
