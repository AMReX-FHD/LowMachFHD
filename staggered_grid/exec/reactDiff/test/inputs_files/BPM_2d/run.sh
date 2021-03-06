#!/bin/bash

EXEC=../../main.Linux.gfortran.mpi.exe
INPUTS=../inputs_paper_BPM_2d
RUNNAME=cs10_Nc64_ExMidTau_dt0.5
RUNNO=$1
NPROC=4

###########
# OPTIONS #
###########

# cross section
OPT1="--prob_hi_z 10."

# resolution
OPT2="--n_cells_x 64 --n_cells_y 64 --max_grid_size_x 32 --max_grid_size_y 32"
#OPT2="--n_cells_x 256 --n_cells_y 256 --max_grid_size_x 32 --max_grid_size_y 32"

# which scheme
OPT3="--temporal_integrator -2 --use_Poisson_rng 1"      # ExMidTau
#OPT3="--temporal_integrator -2 --use_Poisson_rng 2"     # ExMidSSA
#OPT3="--temporal_integrator -4 --use_Poisson_rng 1"     # ImMidTau
#OPT3="--temporal_integrator -4 --use_Poisson_rng 2"     # ImMidSSA
#OPT3="--temporal_integrator  1 --use_Poisson_rng 2 --diffusion_type 3 --reaction_type 0"    # SSA/2+MN+SSA/2

# timestep
OPT4="--fixed_dt 0.5 --max_step  20000 --print_int  20 --n_steps_write_avg 2 --plot_int 200"

# which random number generator 
OPTMISC1="--use_bl_rng F --seed 0"
#OPTMISC1="--use_bl_rng T --seed_diffusion 0 --seed_reaction 0 --seed_init_mass 0"

OPTS="$OPT1 $OPT2 $OPT3 $OPT4 $OPTMISC1"

#######
# RUN #
#######

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

if [ "$RUNNO" = "plot" ]
then
  RUNNAME=${RUNNAME}_plot
  RUNNO=
elif [ "$RUNNO" != "" ]
then
  RUNNAME=${RUNNAME}_RUN$RUNNO
fi

if [ -d $RUNNAME ]
then
  echo "ERROR: $RUNNAME already exists"
  exit
fi

mkdir $RUNNAME
cp $INPUTS $RUNNAME
cp $0 $RUNNAME

cd $RUNNAME

if [ "$RUNNO" = "" ]
then
  echo "mpiexec -n $NPROC ../$EXEC ../$INPUTS $OPTS"
  mpiexec -n $NPROC ../$EXEC ../$INPUTS $OPTS | tee scr_out
else
  echo "mpiexec -n $NPROC ../$EXEC ../$INPUTS $OPTS"
  mpiexec -n $NPROC ../$EXEC ../$INPUTS $OPTS > scr_out &
fi

echo
echo "try: tail -f $RUNNAME/scr_out"
