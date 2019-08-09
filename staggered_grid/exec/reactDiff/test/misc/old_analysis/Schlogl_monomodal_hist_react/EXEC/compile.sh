#!/bin/bash

DIR1=`pwd`
DIR2=../../../../../../src_reactDiff/
DIR3=../../../
EXEC0=main.Linux.gfortran.mpi.exe
INPUTS=../inputs_Schlogl_hist_react

for version in {prev,H,H0,H1,H2}
do
  EXEC=${EXEC0}_$version
  SRCFILE=compute_reaction_rates_${version}.f90

  # remove previous exec and copy src code
  if [ -f $EXEC ]
  then
    rm $EXEC
  fi
  cp $SRCFILE $DIR2/compute_reaction_rates.f90

  # compile and copy exec
  cd $DIR3
  if [ -f $EXEC0 ]
  then
    rm $EXEC0 
  fi
  make
  if [ -f $EXEC0 ]
  then
    cp $EXEC0 $DIR1/$EXEC
  else
    echo "ERROR: no exec obtained ($SRCFILE)"
    exit
  fi

  # run exec 
  cd $DIR1
  ./$EXEC $INPUTS --seed 1 > scr_out_$version
  rm fort.9 fort.10 Hist1D.means.inst.dat
done
