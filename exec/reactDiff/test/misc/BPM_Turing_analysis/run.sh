#!/bin/bash

RUNNAME=TEST
RUNNO=
INPUTSFILE=inputs-256x256-RDME

EXEC=../../main.Linux.gfortran.debug.mpi.exe
NPROC=4

if [ "$RUNNO" != "" ]
then
  RUNNAME=${RUNNAME}_RUN$RUNNO
fi

if [ -d $RUNNAME ]
then
  echo "directory $RUNNAME already exists..."
  sleep 3s
else
  mkdir $RUNNAME
fi

cd $RUNNAME

if [ "$RUNNO" = "" ]
then
  mpiexec -np $NPROC ../$EXEC ../$INPUTSFILE | tee screen_out
else
  echo "mpiexec -np $NPROC ../$EXEC ../$INPUTSFILE > screen_out &"
  mpiexec -np $NPROC ../$EXEC ../$INPUTSFILE > screen_out &
fi
