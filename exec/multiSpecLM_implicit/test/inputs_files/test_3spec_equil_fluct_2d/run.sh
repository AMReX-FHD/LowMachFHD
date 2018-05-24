#!/bin/bash

EXEC=../../main.Linux.gfortran.mpi.exe
INPUTS=inputs_3spec_equil_fluct_2d
RUNNAME=TEST
NPROC=4

# algorithm type
OPT1="--algorithm_type 6"               # boussinesq
#OPT1="--algorithm_type 0"              # inertial trapezoidal
#OPT1="--algorithm_type 2"              # overdamped
#OPT1="--algorithm_type 5"              # inertial midpoint

# smaller time average
#OPT2="--max_step 10  --print_int 10  --n_steps_skip 0 --hydro_grid_int 0 --gmres_verbose 1"
OPT2="--max_step 12000  --print_int 100  --n_steps_skip 2000"
#OPT2="--max_step 120000 --print_int 1000 --n_steps_skip 20000"

# turn off gmres solver (velocity becomes exactly zero)
#OPT3="--variance_coef_mom 0. --initial_variance 0. --gmres_abs_tol 1.e-9 --mg_abs_tol 1.e-9"

OPTS="$OPT1 $OPT2 $OPT3"

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
