#!/bin/bash

EXEC=../../main.Linux.gfortran.mpi.exe
INPUTS=../inputs_dimerization_hist_react_only
RUNNAME=TEST
NRUN=16
NPROC=1

###########
# OPTIONS #
###########

# number of atoms per cell
#OPT1="--prob_hi_z 20. --rate_const_1 1. --rate_const_2 1.290564"       # N0=20
OPT1="--prob_hi_z 40. --rate_const_1 1. --rate_const_2 1.311528"        # N0=40
#OPT1="--prob_hi_z 60. --rate_const_1 1. --rate_const_2 1.318704"       # N0=60

# tau leaping or CLE
OPT2="--use_Poisson_rng 1"
#OPT2="--use_Poisson_rng 0"

# rate multiplier
OPT3="--rate_multiplier 0.021"

# time step size
OPT4="--fixed_dt 1.e-2"

# task (histogram)
OPT5="--max_step 100000  --plot_int 0   --print_int 1000 --hydro_grid_int 100 --n_steps_skip 20000 --histogram_unit 10"

#OPTMISC1="--max_grid_size_x 16 --max_grid_size_y 16"
#OPTMISC2=

OPTS="$OPT1 $OPT2 $OPT3 $OPT4 $OPT5 $OPTMISC1 $OPTMISC2"

#######
# RUN #
#######

# check executable / inputs file / run directory

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

# multiple run

for ((i=1;i<=$NRUN;i++))
do
  RUNDIR=${RUNNAME}_RUN$i

  if [ -d $RUNDIR ]
  then
    echo "ERROR: $RUNDIR already exists"
    exit
  fi

  echo $RUNDIR
  mkdir $RUNDIR
  cp $INPUTS $RUNDIR
  cp $0 $RUNDIR

  cd $RUNDIR
  echo "mpiexec -n $NPROC ../$EXEC ../$INPUTS $OPTS > scr_out &"
  mpiexec -n $NPROC ../$EXEC ../$INPUTS $OPTS > scr_out &

  sleep 2s
  cd ..
done

echo
echo "try: tail -f $RUNDIR/scr_out"
