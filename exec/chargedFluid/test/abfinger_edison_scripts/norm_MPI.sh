#!/bin/bash

# request an interative job and run this script
# salloc --qos=debug --time=30 --nodes=1

# norm_MPI.ex is obtained from amrex/Tools/C_util/Convergence/DiffSameDomainRefined.cpp (3d,MPI)
# need to comment:
# //(*error[iLevel])[mfi].minus(data2Coarse  , 0, iComp, 1);

RUNNAME=TEST
RES=${RUNNAME}.norm

for ((i=0;i<=5000;i+=20))
do
  printf -v filename "plt%08d" $i
  fullname=$RUNNAME/$filename
  echo "*** $fullname ***"

  #./norm.ex infile1= $fullname reffile= $fullname norm=1 | tail -1 | tee -a $RES
  srun -n 24 ./norm_MPI.ex infile1= $fullname reffile= $fullname norm=1 | tail -1 | tee -a $RES
done
