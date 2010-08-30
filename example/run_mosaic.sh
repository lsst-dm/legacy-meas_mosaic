#!/bin/sh
#PBS -l ncpus=8
#PBS -l nodes=1:hsca-06:ppn=8
#PBS -q default
#
OMP_NUM_THREADS=6; export OMP_NUM_THREADS
if [ -n "$PBS_O_WORKDIR" ]
then
  cd $PBS_O_WORKDIR
fi
#
setup -r /home/yasuda/temp/hscMosaic

python ./run_mosaic.py
