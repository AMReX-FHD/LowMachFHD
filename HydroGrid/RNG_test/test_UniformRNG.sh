#!/bin/bash

EXEC=./main.Linux.gfortran.debug.mpi.exe

# You can change these options
RNG_ENG_TYPE=1                         # RNG_ENG_TYPE = 1 (HydroGrid), 2 (c++)
SEED=0                                 # SEED: random seed (0=randomly chosen)
PREC=1                                 # PREC: 1 (UniformRNG_sp), 2 (UniformRNG_dp)

HIST_N_SAMPLES=10000000                # how many samples for histogram
HIST_N_BINS=1000                       # number of bins
HIST_FILE=res.hist                     # output filename

# OPTS
OPT0="--rng_eng_type $RNG_ENG_TYPE --seed $SEED"
OPT1="--dist_type 1 --prec $PREC"
OPT2="--hist_n_samples $HIST_N_SAMPLES --hist_n_bins $HIST_N_BINS --hist_left_end 0. --hist_right_end 1. --hist_is_center_val F --hist_file $HIST_FILE"
OPTS="$OPT0 $OPT1 $OPT2"

if [ ! -f $EXEC ]
then
  echo "ERROR: $EXEC does not exist"
  exit 
fi

echo "$EXEC $OPTS"
$EXEC $OPTS
