#!/bin/bash

EXEC=../../main.Linux.gfortran.mpi.exe
INPUTS=../inputs_dimerization_hist
RUNNAME=TEST
NRUN=16
NPROC=1

###########
# OPTIONS #
###########

# number of atoms per cell
#OPT1="--rho0 20. --rate_const_1 1. --rate_const_2 1.290564"    # N0=20
OPT1="--rho0 40. --rate_const_1 1. --rate_const_2 1.311528"     # N0=40
#OPT1="--rho0 60. --rate_const_1 1. --rate_const_2 1.318704"    # N0=60

# reaction on/off
#OPT2="--rate_multiplier 0."            # reaction off
#OPT2="--rate_multiplier 0.421875"      # N0=20, d = sqrt(10)*dx, k1 = 27*N0/1280
OPT2="--rate_multiplier 0.84375"        # N0=40
#OPT2="--rate_multiplier 1.265625"      # N0=60

# advection type
OPT3="--advection_type 0"

# temperature
OPT4="--T_init_1 1.e4 --T_init_2 1.e4"

# time step size
OPT5="--fixed_dt 1.e-2 --max_step 100000 --plot_int 0 --print_int 1000 --hydro_grid_int 100 --n_steps_skip 20000 --histogram_unit 10"

# Stratonovich or Ito
OPT6="--midpoint_stoch_mass_flux_type 1"
#OPT6="--midpoint_stoch_mass_flux_type 2"

#OPTMISC1="--max_grid_size_x 16 --max_grid_size_y 16"
#OPTMISC2="--use_Poisson_rng 0"

OPTS="$OPT1 $OPT2 $OPT3 $OPT4 $OPT5 $OPT6 $OPTMISC1 $OPTMISC2"

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
