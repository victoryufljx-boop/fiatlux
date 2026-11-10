#!/bin/bash
# wait for rep $1 to finish then start rep $2 seed $3
while [ ! -f rep$1/t_end ]; do sleep 5; done
./run_rep.sh $2 $3
