#!/bin/bash

EXEC=../../main.Linux.gfortran.mpi.exe
INPUTS=../inputs_paper_equilibrium_1d_hist_schemes
RUNNAME=ExMidTau_dt0.1
NRUN=16

###########
# OPTIONS #
###########

# which scheme
OPT1="--temporal_integrator -2 --use_Poisson_rng 1"   # ExMidTau
#OPT1="--temporal_integrator -2 --use_Poisson_rng 2"   # ExMidSSA
#OPT1="--temporal_integrator -4 --use_Poisson_rng 1"   # ImMidTau
#OPT1="--temporal_integrator -4 --use_Poisson_rng 2"   # ImMidSSA
#OPT1="--temporal_integrator  1 --use_Poisson_rng 2 --diffusion_type 3 --reaction_type 0" # SSA/2+MN+SSA/2

# timestep size
#OPT2="--fixed_dt 5.    --max_step 10200   --print_int 200    --hydro_grid_int 1   --n_steps_skip 200"
#OPT2="--fixed_dt 2.    --max_step 10500   --print_int 500    --hydro_grid_int 1   --n_steps_skip 500"
#OPT2="--fixed_dt 1.    --max_step 11000   --print_int 1000   --hydro_grid_int 1   --n_steps_skip 1000"
#OPT2="--fixed_dt 0.5   --max_step 22000   --print_int 2000   --hydro_grid_int 2   --n_steps_skip 2000"
#OPT2="--fixed_dt 0.4   --max_step 27500   --print_int 2500   --hydro_grid_int 5   --n_steps_skip 2500"
#OPT2="--fixed_dt 0.25  --max_step 44000   --print_int 4000   --hydro_grid_int 4   --n_steps_skip 4000"
#OPT2="--fixed_dt 0.2   --max_step 55000   --print_int 5000   --hydro_grid_int 5   --n_steps_skip 5000"
OPT2="--fixed_dt 0.1   --max_step 110000  --print_int 10000  --hydro_grid_int 10  --n_steps_skip 10000"
#OPT2="--fixed_dt 0.05  --max_step 220000  --print_int 20000  --hydro_grid_int 20  --n_steps_skip 20000"
#OPT2="--fixed_dt 0.025 --max_step 440000  --print_int 40000  --hydro_grid_int 40  --n_steps_skip 40000"
#OPT2="--fixed_dt 0.02  --max_step 550000  --print_int 50000  --hydro_grid_int 50  --n_steps_skip 50000"
#OPT2="--fixed_dt 0.01  --max_step 1100000 --print_int 100000 --hydro_grid_int 100 --n_steps_skip 100000"

# midpoint_stoch_flux_type: 1, 2, or 3
OPT3="--midpoint_stoch_flux_type 3"

# which random number generator 
OPTMISC1="--use_bl_rng F --seed 0"
#OPTMISC1="--use_bl_rng T --seed_diffusion 0 --seed_reaction 0 --seed_init 0"

OPTS="$OPT1 $OPT2 $OPT3 $OPTMISC1"

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

if [ $NRUN == 1 ]
then
  if [ -d $RUNNAME ]
  then
    echo "ERROR: $RUNNAME already exists"
    exit
  fi

  mkdir $RUNNAME
  cp $INPUTS $RUNNAME
  cp $0 $RUNNAME

  cd $RUNNAME
  ../$EXEC ../$INPUTS $OPTS | tee scr_out

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
  echo "mpiexec -n 1 ../$EXEC ../$INPUTS $OPTS > scr_out &"
  mpiexec -n 1 ../$EXEC ../$INPUTS $OPTS > scr_out &

  sleep 2s
  cd ..
done

tail -f $RUNDIR/scr_out

