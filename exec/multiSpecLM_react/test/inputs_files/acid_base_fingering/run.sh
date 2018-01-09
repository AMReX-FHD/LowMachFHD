#!/bin/bash

EXEC=../../main.Linux.gfortran.mpi.exe
INPUTS=../inputs_asymm_convect_fingers
RUNNAME=TEST
NPROC=8

# domain size
#OPT1="--prob_hi_x 0.8 --n_cells_x 128 --max_grid_size_x 64 --prob_hi_y 1.6 --n_cells_y 256 --max_grid_size_y 64"
OPT1="--prob_hi_x 0.4 --n_cells_x 64  --max_grid_size_x 32 --prob_hi_y 0.8 --n_cells_y 128 --max_grid_size_y 32"

# time step
OPT2="--fixed_dt 0.05 --max_step 400 --print_int 10 --plot_int 10"
#OPT2="--fixed_dt 0.1  --max_step 200 --print_int 5  --plot_int 5"
#OPT2="--fixed_dt 0.15 --max_step 150 --print_int 3  --plot_int 3"

# reaction
OPT3="--nreactions 1 --use_Poisson_rng -1"     # det reactions
#OPT3="--nreactions 1 --use_Poisson_rng 0"      # cle reactions
#OPT3="--nreactions 1 --use_Poisson_rng 1"      # poisson reactions
#OPT3="--nreactions 0"                          # no reaction

# rate multiplier
OPT4="--rate_multiplier 1."

OPTS="$OPT1 $OPT2 $OPT3 $OPT4"

#######
# RUN #
#######

# check executable and inputs file

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

# single run

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
