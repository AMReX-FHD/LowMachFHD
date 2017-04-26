#!/bin/bash

EXEC=../../main.Linux.gfortran.mpi.exe
INPUTS=../inputs_paper_3spec_equil_fluct_2d
RUNNAME=alg0_short
NPROC=4

# algorithm type
OPT1="--algorithm_type 0"               # inertial trapezoidal
#OPT1="--algorithm_type 2"              # overdamped
#OPT1="--algorithm_type 5"              # inertial midpoint

# smaller time average
OPT2="--max_step 12000  --print_int 100  --n_steps_skip 2000"
#OPT2="--max_step 120000 --print_int 1000 --n_steps_skip 20000"

OPTS="$OPT1 $OPT2"

#######
# RUN #
#######

if [ ! -f $EXEC ]
then
  echo "ERROR: $EXEC does not exist"
  exit
fi

if [ ! -f $INPUTS ]
then
  echo "ERROR: $INPUTS does not exist"
  exit
fi

if [ -d $RUNNAME ]
then
  echo "ERROR: $RUNNAME already exists"
  exit
fi

mkdir $RUNNAME
cp $INPUTS $RUNNAME
cp $0 $RUNNAME

cd $RUNNAME

mpiexec -n $NPROC ../$EXEC ../$INPUTS $OPTS | tee scr_out
