#!/bin/sh
#PBS -l ncpus=1
#PBS -l nodes=1:hsca-06
#PBS -q default
#
cd $PBS_O_WORKDIR
#
setup -r /home/yasuda/work/hscMosaic

python ./run_stack.py All
