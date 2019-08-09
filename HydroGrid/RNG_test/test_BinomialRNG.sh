#!/bin/bash

EXEC=./main.Linux.gfortran.debug.mpi.exe

# You can change these options
RNG_ENG_TYPE=1                         # RNG_ENG_TYPE = 1 (HydroGrid), 2 (c++)
SEED=0                                 # SEED: random seed (0=randomly chosen)
PREC=1                                 # PREC: 1 (BinomialRNG_sp), 2 (BinomialRNG_dp)

PARAM_N_TRIALS=100                     # n_trials of Binomial distribution
PARAM_SUCCESS_PROB=0.2                 # success_prob of Binomial distribution

HIST_N_SAMPLES=10000000                # how many samples for histogram
HIST_LEFT_END=0                        # left-end value of histogram
HIST_RIGHT_END=100                     # right-end value of histogram
HIST_FILE=res.hist                     # output filename
((HIST_N_BINS=HIST_RIGHT_END-HIST_LEFT_END+1))  # number of bins

# OPTS
OPT0="--rng_eng_type $RNG_ENG_TYPE --seed $SEED"
OPT1="--dist_type 4 --prec $PREC --param_n_trials $PARAM_N_TRIALS --param_success_prob $PARAM_SUCCESS_PROB"
OPT2="--hist_n_samples $HIST_N_SAMPLES --hist_n_bins $HIST_N_BINS --hist_left_end $HIST_LEFT_END --hist_right_end $HIST_RIGHT_END --hist_is_center_val T --hist_file $HIST_FILE"
OPTS="$OPT0 $OPT1 $OPT2"

if [ ! -f $EXEC ]
then
  echo "ERROR: $EXEC does not exist"
  exit 
fi

echo "$EXEC $OPTS"
$EXEC $OPTS
