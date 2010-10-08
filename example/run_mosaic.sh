#!/bin/sh
#PBS -l ncpus=10
#PBS -l nodes=1:ppn=10
#PBS -q default
#
OMP_NUM_THREADS=10; export OMP_NUM_THREADS
if [ -n "$PBS_O_WORKDIR" ]
then
  cd $PBS_O_WORKDIR
fi
#
setup -r /home/yasuda/temp/hscMosaic

python ./run_mosaic.py COSMOS_0
#python ./run_mosaic.py SXDS
#python ./run_mosaic.py CFHQS
#python ./run_mosaic.py CFHTLS-D3
#python ./run_mosaic.py GALPLANE
#python ./run_mosaic.py LOCKMANHOLE
