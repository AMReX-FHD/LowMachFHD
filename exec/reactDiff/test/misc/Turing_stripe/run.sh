#!/bin/bash

RUNNAME=TEST

EXEC=../../main.Linux.gfortran.debug.mpi.exe
INPUTS=inputs_Turing_stripe_3d

if [ -d $RUNNAME ]
then
  echo "directory $RUNNAME already exisits"
else
  mkdir $RUNNAME
fi

cd $RUNNAME
mpiexec -np 4 ../$EXEC ../$INPUTS
