#!/bin/sh
#PBS -l ncpus=10
#PBS -l nodes=1:hsca-06:ppn=10
#PBS -q default
#
OMP_NUM_THREADS=10; export OMP_NUM_THREADS
if [ -n "$PBS_O_WORKDIR" ]
then
  cd $PBS_O_WORKDIR
fi
#
setup -r /home/yasuda/temp/hscMosaic

python ./run_mosaic.py
