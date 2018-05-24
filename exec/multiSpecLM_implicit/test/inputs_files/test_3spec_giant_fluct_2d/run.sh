#!/bin/bash

EXEC=../../main.Linux.gfortran.mpi.exe
INPUTS=inputs_3spec_giant_fluct_2d
RUNNAME=TEST
NPROC=4

# visc_coef
OPT1="--visc_coef 1.e3"

# algorithm type
OPT2="--algorithm_type 6"               # Boussinesq 
#OPT2="--algorithm_type 0"              # inertial trapezoidal
#OPT2="--algorithm_type 2"              # overdamped
#OPT2="--algorithm_type 5"              # inertial midpoint

# smaller system
OPT3="--prob_hi_x 64. --prob_hi_y 32. --n_cells_x 64 --n_cells_y 32"

# max_grid_size
OPT4="--max_grid_size_x 32 --max_grid_size_y 16"

# smaller time average
OPT5="--max_step 10000  --print_int 100 --n_steps_skip 2500  --stats_int 1000"
#OPT5="--max_step 100000 --print_int 1000 --n_steps_skip 25000 --stats_int 10000"

OPTS="$OPT1 $OPT2 $OPT3 $OPT4 $OPT5"

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
