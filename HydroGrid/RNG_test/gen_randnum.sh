#!/bin/bash

EXEC=./main.Linux.gfortran.debug.mpi.exe

# you can generate a number of random numbers and save them in a file
# for a quick histogram test, you can use gnuplot (see gnuplot_cmd_quick_histogram file)

# OPT0 
RNG_ENG_TYPE=1                         # RNG_ENG_TYPE = 1 (HydroGrid), 2 (c++)
SEED=0                                 # SEED: random seed (0=randomly chosen)
OPT0="--rng_eng_type $RNG_ENG_TYPE --seed $SEED"

# OPT1 
DIST_TYPE=1                            # DIST_TYPE = 1(uniform),2(std normal),3(Poisson),4(binomial)
PREC=1                                 # PREC: 1 (single precision), 2 (double precision)
OPT1="--dist_type $DIST_TYPE --prec $PREC"

# OPT2 (for Poisson, DIST_TYPE=3)
#PARAM_MEAN=10.                        # mean of Poisson distribution
#OPT2="--param_mean $PARAM_MEAN"

# OPT2 (for binomial, DIST_TYPE=4)
#PARAM_N_TRIALS=100                    # n_trials of Binomial distribution
#PARAM_SUCCESS_PROB=0.2                # success_prob of Binomial distribution
#OPT2="--param_n_trials $PARAM_N_TRIALS --param_success_prob $PARAM_SUCCESS_PROB"

# OPT3 
GEN_RANDNUM_INT=1000000                # how many random numbers to be generated
GEN_RANDNUM_FILE=res.gen_randnum       # output filename 
OPT3="--gen_randnum_int $GEN_RANDNUM_INT --gen_randnum_file $GEN_RANDNUM_FILE"

# OPTS
OPTS="$OPT0 $OPT1 $OPT2 $OPT3"

if [ ! -f $EXEC ]
then
  echo "ERROR: $EXEC does not exist"
  exit 
fi

echo "$EXEC $OPTS"
$EXEC $OPTS
