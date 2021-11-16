#!/bin/bash -l
source /users_home/opa/ww3_cst/setup.sh

area=
cd /users_home/opa/ww3_cst/models/$area/tasks
module list
bsub -P 0410 -e /work/opa/ww3_cst/runs/$area/`date +%Y%m%d`/main.err -o /work/opa/ww3_cst/runs/$area/`date +%Y%m%d`/main.out python /users_home/opa/ww3_cst/models/$area/tasks/runme.py  -s `date +%Y%m%d` -d 5

