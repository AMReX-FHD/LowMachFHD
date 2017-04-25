#!/bin/bash

EXEC=../../main.Linux.gfortran.mpi.exe
INPUTS=../inputs_paper_Schlogl_1d_structure_factor
RUNNAME=ImMidTau_b5_a0.5

###########
# OPTIONS #
###########

# which scheme
#OPT1="--temporal_integrator -2 --use_Poisson_rng 1"     # ExMidTau
#OPT1="--temporal_integrator -2 --use_Poisson_rng 2"     # ExMidSSA
OPT1="--temporal_integrator -4 --use_Poisson_rng 1"      # ImMidTau
#OPT1="--temporal_integrator -4 --use_Poisson_rng 2"     # ImMidSSA
#OPT1="--temporal_integrator  1 --use_Poisson_rng 2 --diffusion_type 3 --reaction_type 0"    # SSA/2+MN+SSA/2
#OPT1="--temporal_integrator  2 --use_Poisson_rng 1 --diffusion_type 3 --reaction_type 1"    # MN/2+Tau2nd+MN/2
#OPT1="--temporal_integrator  2 --use_Poisson_rng 2 --diffusion_type 3 --reaction_type 0"    # MN/2+SSA+MN/2
#OPT1="--temporal_integrator  2 --use_Poisson_rng 1 --diffusion_type 1 --reaction_type 1"    # CN/2+Tau2nd+CN/2
#OPT1="--temporal_integrator  2 --use_Poisson_rng 2 --diffusion_type 1 --reaction_type 0"    # CN/2+SSA+CN/2

# alpha and beta
#DFICK=0.25              # beta=0.25
DFICK=5.                 # beta=5.
#RM=0.0142857            # alpha=0.025
RM=0.285714              # alpha=0.5
OPT2="--D_Fick_1 $DFICK --rate_multiplier $RM"

# time average 
OPT3="--max_step 60000  --n_steps_skip 10000 --hydro_grid_int 1 --print_int 1000"
#OPT3="--max_step 510000 --n_steps_skip 10000 --hydro_grid_int 1 --print_int 1000"

# which random number generator 
OPTMISC1="--use_bl_rng F --seed 0"
#OPTMISC1="--use_bl_rng T --seed_diffusion 0 --seed_reaction 0 --seed_init 0"

OPTS="$OPT1 $OPT2 $OPT3 $OPTMISC1"

#######
# RUN #
#######

PYSCR=Sk_xi.py

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

python ../$PYSCR $RM $DFICK
