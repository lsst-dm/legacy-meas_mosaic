#!/bin/sh
#PBS -l ncpus=1
#PBS -l nodes=1
#PBS -q default
#
cd $PBS_O_WORKDIR
#

python ./run_stack.py --rerun=yasuda --instrument=suprimecam --program=ACTJ0022M0036 --filter=W-S-R+ All
