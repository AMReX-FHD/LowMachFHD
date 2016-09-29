#!/bin/bash

###########
# OPTIONS #
###########

RUNNAME=TEST
NRUN=4

# diff-only or reaction-diffusion
OPT1="--nreactions 0"
#OPT1="--nreactions 4 --rate_multiplier 0.1"

# system paras: cross_section, dx, Ncell 
CS=10.
DX=1.
NCELL=64
DV=`echo print $CS*$DX | python`
LX=`echo print $NCELL*$DX | python`
OPT2="--prob_hi_x $LX --prob_hi_y $DX --n_cells_x $NCELL --n_cells_y 1 --cross_section $CS"

# numerical scheme
#OPT3="--temporal_integrator 0  --diffusion_type 3   --reaction_type 0  --use_Poisson_rng 2"
#OPT3="--temporal_integrator 2  --diffusion_type 3   --reaction_type 0  --use_Poisson_rng 2"
OPT3="--temporal_integrator -1 --use_Poisson_rng 1  --avg_type 1"
#OPT3="--temporal_integrator -2 --use_Poisson_rng 1  --avg_type 1  --midpoint_stoch_flux_type 3"
#OPT3="--temporal_integrator -4 --use_Poisson_rng 1  --avg_type 1  --midpoint_stoch_flux_type 3"

# initial configuration
OPT4="--initial_variance -1. --integer_populations T"

# time step
#OPT5="--fixed_dt 5.    --max_step 10200    --print_int 200     --hydro_grid_int 1    --n_steps_skip 200"
#OPT5="--fixed_dt 2.    --max_step 10500    --print_int 500     --hydro_grid_int 1    --n_steps_skip 500"
#OPT5="--fixed_dt 1.    --max_step 11000    --print_int 1000    --hydro_grid_int 1    --n_steps_skip 1000"
#OPT5="--fixed_dt 0.5   --max_step 22000    --print_int 2000    --hydro_grid_int 2    --n_steps_skip 2000"
#OPT5="--fixed_dt 0.4   --max_step 27500    --print_int 2500    --hydro_grid_int 5    --n_steps_skip 2500"
#OPT5="--fixed_dt 0.25  --max_step 44000    --print_int 4000    --hydro_grid_int 4    --n_steps_skip 4000"
#OPT5="--fixed_dt 0.2   --max_step 55000    --print_int 5000    --hydro_grid_int 5    --n_steps_skip 5000"
#OPT5="--fixed_dt 0.1   --max_step 110000   --print_int 10000   --hydro_grid_int 10   --n_steps_skip 10000"
#OPT5="--fixed_dt 0.05  --max_step 220000   --print_int 20000   --hydro_grid_int 20   --n_steps_skip 20000"
#OPT5="--fixed_dt 0.025 --max_step 440000   --print_int 40000   --hydro_grid_int 40   --n_steps_skip 40000"
#OPT5="--fixed_dt 0.02  --max_step 550000   --print_int 50000   --hydro_grid_int 50   --n_steps_skip 50000"
OPT5="--fixed_dt 0.01  --max_step 1100000  --print_int 100000  --hydro_grid_int 100  --n_steps_skip 100000"
#OPT5="--fixed_dt 0.005 --max_step 2200000  --print_int 200000  --hydro_grid_int 200  --n_steps_skip 200000"
#OPT5="--fixed_dt 0.002 --max_step 5500000  --print_int 500000  --hydro_grid_int 500  --n_steps_skip 500000"
#OPT5="--fixed_dt 0.001 --max_step 11000000 --print_int 1000000 --hydro_grid_int 1000 --n_steps_skip 1000000"
#OPT5="--fixed_dt 0.0001 --max_step 110000000 --print_int 10000000 --hydro_grid_int 10000 --n_steps_skip 10000000"

# rng
OPT6="--use_bl_rng F --seed 0"
#OPT6="--use_bl_rng T --seed_diffusion 0 --seed_reaction 0 --seed_init 0"

#######
# RUN #
#######

EXEC=../../main.Linux.gfortran.mpi.exe
INPUTS=inputs_Schlogl_hist_1d
OPTS="$OPT1 $OPT2 $OPT3 $OPT4 $OPT5 $OPT6"

for ((i=1;i<=$NRUN;i++))
do
  RUNDIR=${RUNNAME}_RUN$i
  echo $RUNDIR 

  mkdir $RUNDIR
  cp $INPUTS $RUNDIR
  cp run.sh $RUNDIR
  
  cd $RUNDIR
  echo "mpiexec -n 1 ../$EXEC ../$INPUTS $OPTS > scr_out &"
  mpiexec -n 1 ../$EXEC ../$INPUTS $OPTS > scr_out &

  sleep 3s
  cd ..
done
