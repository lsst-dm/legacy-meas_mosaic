#!/bin/sh
#PBS -l ncpus=10
#PBS -l nodes=1:ppn=10
#PBS -q default
#PBS -j oe
#
OMP_NUM_THREADS=10; export OMP_NUM_THREADS
if [ -n "$PBS_O_WORKDIR" ]
then
  cd $PBS_O_WORKDIR
fi
#
source /data/ana/products-yasuda/loadLSST.sh
setup -t price-DC2 -r /home/yasuda/temp/hscMosaic
#setup pipette -t DC2.2 -t DC2

python $HSCMOSAIC_DIR/example/run_mosaic.py --rerun=yasuda-test2 --frameid=23 --outputDir="/data/yasuda/DC2/yasuda-test2"
