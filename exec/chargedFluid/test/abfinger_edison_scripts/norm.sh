#!/bin/bash -l

# norm_MPI.ex is obtained from amrex/Tools/C_util/Convergence/DiffSameDomainRefined.cpp (3d,MPI)
# need to comment:
# //(*error[iLevel])[mfi].minus(data2Coarse  , 0, iComp, 1);

NPROC=$1
RUNNAME=$2

RES=${RUNNAME}.norm

for ((i=0;i<=5000;i+=20))
do
  printf -v filename "plt%08d" $i
  fullname=$RUNNAME/$filename
  echo "*** $fullname ***"

  srun -n $NPROC ./norm_MPI.ex infile1= $fullname reffile= $fullname norm=1 | tail -1 | tee -a $RES
done
