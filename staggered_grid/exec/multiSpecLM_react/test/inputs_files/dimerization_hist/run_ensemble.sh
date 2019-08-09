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
OPT1="--rho0 40. --rate_const_1 0.7755 --rate_const_2 1."               # N0=40, K=0.7755

# reaction on/off
OPT2="--rate_multiplier 1.125 --include_discrete_LMA_correction T"      # N0=40, d=sqrt(10)*dx, k2=9*N0/320
#OPT2="--nreactions 0"                                                  # reaction off

# advection on/off
OPT3="--variance_coef_mom 1. --initial_variance_mom 1. --advection_type 0"
#OPT3="--variance_coef_mom 0. --initial_variance_mom 0. --gmres_abs_tol 1.e-10 --mg_abs_tol 1.e-10"

# temperature
OPT4="--T_init_1 1.e3 --T_init_2 1.e3"

# time step size
OPT5="--fixed_dt 1.e-2 --max_step 100000 --plot_int 0 --print_int 1000 --hydro_grid_int 10 --n_steps_skip 10000 --histogram_unit 10"

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
