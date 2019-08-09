#!/bin/bash

###########
# OPTIONS #
###########

EXEC=main.Linux.gfortran.mpi.exe_prev
INPUTS=inputs_Schlogl_hist_react

RUNNAME=prev_ImMidTau_dt1
NRUN=16

# how fast diffusion and reaction
OPT1="--D_Fick_1 1. --rate_multiplier 0.005"

# system geometry
CS=5.
DX=1.
NCELL=4
DV=`echo print $CS*$DX | python`
LX=`echo print $NCELL*$DX | python`
OPT2="--prob_hi_x $LX --prob_hi_y $DX --n_cells_x $NCELL --n_cells_y 1 --cross_section $CS"

# numerical scheme
#OPT3="--temporal_integrator -2 --use_Poisson_rng 1  --avg_type 1  --midpoint_stoch_flux_type 3"   # ExMidTau 
OPT3="--temporal_integrator -4 --use_Poisson_rng 1  --avg_type 1  --midpoint_stoch_flux_type 3"   # ImMidTau 

# initial configuration
OPT4="--initial_variance -1. --integer_populations T"
#OPT4="--initial_variance -1. --integer_populations F"

# time step
#OPT5="--fixed_dt 10.   --max_step 101000    --print_int 100     --hydro_grid_int 1    --n_steps_write_avg 1    --n_steps_skip 1000"
#OPT5="--fixed_dt 5.    --max_step 102000    --print_int 200     --hydro_grid_int 1    --n_steps_write_avg 1    --n_steps_skip 2000"
#OPT5="--fixed_dt 2.    --max_step 105000    --print_int 500     --hydro_grid_int 1    --n_steps_write_avg 1    --n_steps_skip 5000"
OPT5="--fixed_dt 1.    --max_step 110000    --print_int 1000    --hydro_grid_int 1    --n_steps_write_avg 1    --n_steps_skip 10000"
#OPT5="--fixed_dt 0.5   --max_step 220000    --print_int 2000    --hydro_grid_int 2    --n_steps_write_avg 2    --n_steps_skip 20000"
#OPT5="--fixed_dt 0.2   --max_step 550000    --print_int 5000    --hydro_grid_int 5    --n_steps_write_avg 5    --n_steps_skip 50000"
#OPT5="--fixed_dt 0.1   --max_step 1100000   --print_int 10000   --hydro_grid_int 10   --n_steps_write_avg 10   --n_steps_skip 100000"
#OPT5="--fixed_dt 0.05  --max_step 2200000   --print_int 20000   --hydro_grid_int 20   --n_steps_write_avg 20   --n_steps_skip 200000"
#OPT5="--fixed_dt 0.02  --max_step 5500000   --print_int 50000   --hydro_grid_int 50   --n_steps_write_avg 50   --n_steps_skip 500000"
#OPT5="--fixed_dt 0.01  --max_step 11000000  --print_int 100000  --hydro_grid_int 100  --n_steps_write_avg 100  --n_steps_skip 1000000"

# rng
OPT6="--use_bl_rng F --seed 0"
#OPT6="--use_bl_rng T --seed_diffusion 0 --seed_reaction 0 --seed_init_mass 0"

#######
# RUN #
#######

OPTS="$OPT1 $OPT2 $OPT3 $OPT4 $OPT5 $OPT6"

for ((i=1;i<=$NRUN;i++))
do
  RUNDIR=${RUNNAME}_RUN$i

  if [ -d $RUNDIR ]
  then
    echo "ERROR: $RUNDIR already exists..."
    exit
  fi

  echo $RUNDIR 
  mkdir $RUNDIR
  cp $INPUTS $RUNDIR
  cp run.sh $RUNDIR
  
  cd $RUNDIR
  echo "mpiexec -n 1 ../$EXEC ../$INPUTS $OPTS > scr_out &"
  mpiexec -n 1 ../$EXEC ../$INPUTS $OPTS > scr_out &

  sleep 1s
  cd ..
done

echo "try: tail -f ${RUNNAME}_RUN1/scr_out"
