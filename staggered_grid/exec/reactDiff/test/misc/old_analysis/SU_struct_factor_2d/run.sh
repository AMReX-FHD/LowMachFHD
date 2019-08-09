#!/bin/bash

EXEC=../../main.Linux.gfortran.mpi.exe
INPUTS=../../inputs_SU_2d

RUNNAME=TEST

# D_Fick
OPT0="--D_Fick_1 1.e-6 --D_Fick_2 1.e-6"

# rate_multiplier
OPT1="--rate_multiplier 1."

# grid 
OPT2="--prob_hi_x 6.4e-5 --prob_hi_y 6.4e-5 --n_cells_x 64 --n_cells_y 64 --max_grid_size_x 32 --max_grid_size_y 32"
#OPT2="--prob_hi_x 6.4e-5 --prob_hi_y 6.4e-5 --n_cells_x 128 --n_cells_y 128 --max_grid_size_x 64 --max_grid_size_y 64"

# cross_section
OPT3="--cross_section 1.e-5"
#OPT3="--cross_section 4.e-5"

# timestep 
OPT4="--max_step 10000 --plot_int 0 --print_int 100 --hydro_grid_int 1 --n_steps_skip 1000 --fixed_dt 1.e-7"
#OPT4="--max_step 40000 --plot_int 0 --print_int 400 --hydro_grid_int 4 --n_steps_skip 4000 --fixed_dt 0.25e-7"

# scheme
OPT5="--temporal_integrator -2 --use_Poisson_rng 1"
#OPT2="--temporal_integrator 1 --diffusion_type 3 --reaction_type 0 --use_Poisson_rng 2"

#####

if [ ! -f $EXEC ]
then
  echo "ERROR: $EXEC does not exist"
  exit
fi

RUNNAME=RUN_$RUNNAME
if [ -d $RUNNAME ]
then
  echo "ERROR: $RUNNAME exists"
  exit
else
  mkdir $RUNNAME
  cp $0 $RUNNAME
  cp $INPUTS $RUNNAME
  cd $RUNNAME
fi

OPTS="$OPT0 $OPT1 $OPT2 $OPT3 $OPT4 $OPT5"

echo "mpiexec -n 4 ../$EXEC ../$INPUTS $OPTS | tee scr_out"
mpiexec -n 4 ../$EXEC ../$INPUTS $OPTS | tee scr_out

#####

grep n_avg scr_out | awk '{print $4, $5}' > res.n_avg

if [ -f SU.S_k.pair=1.Re.dat ]
then
  sed -i 's/0.00000000/#0.00000000/g' SU.S_k.pair=*.Re.dat
fi
