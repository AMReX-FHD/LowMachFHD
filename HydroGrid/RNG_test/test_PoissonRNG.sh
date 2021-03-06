#!/bin/bash

EXEC=./main.Linux.gfortran.debug.mpi.exe

# You can change these options
RNG_ENG_TYPE=1                         # RNG_ENG_TYPE = 1 (HydroGrid), 2 (c++)
SEED=0                                 # SEED: random seed (0=randomly chosen)
PREC=1                                 # PREC: 1 (PoissonRNG_sp), 2 (PoissonRNG_dp)

PARAM_MEAN=10.                         # mean of Poisson distribution

HIST_N_SAMPLES=10000000                # how many samples for histogram
HIST_LEFT_END=0                        # left-end value of histogram
HIST_RIGHT_END=30                      # right-end value of histogram
HIST_FILE=res.hist                     # output filename
((HIST_N_BINS=HIST_RIGHT_END-HIST_LEFT_END+1))  # number of bins

# OPTS
OPT0="--rng_eng_type $RNG_ENG_TYPE --seed $SEED"
OPT1="--dist_type 3 --prec $PREC --param_mean $PARAM_MEAN"
OPT2="--hist_n_samples $HIST_N_SAMPLES --hist_n_bins $HIST_N_BINS --hist_left_end $HIST_LEFT_END --hist_right_end $HIST_RIGHT_END --hist_is_center_val T --hist_file $HIST_FILE"
OPTS="$OPT0 $OPT1 $OPT2"

if [ ! -f $EXEC ]
then
  echo "ERROR: $EXEC does not exist"
  exit 
fi

echo "$EXEC $OPTS"
$EXEC $OPTS
