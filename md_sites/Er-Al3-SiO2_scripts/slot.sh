#!/bin/bash
# usage: slot.sh rep1 rep2 ...
cd /home/claude/md
export OMP_NUM_THREADS=1
for r in "$@"; do
  python gen.py $((1000+r)) rep$r.data > rep$r.gen
  date -u +%s > rep$r.start
  lmp -in md.in -var data rep$r.data -var seed $((5000+7*r)) -var tag rep$r -var nquench 320000 -log rep$r.log -screen none
  date -u +%s > rep$r.end
done
