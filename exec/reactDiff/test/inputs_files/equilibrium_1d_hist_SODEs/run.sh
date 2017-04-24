#!/bin/bash

EXEC=../../main.Linux.gfortran.mpi.exe
INPUTS=../inputs_paper_equilibrium_1d_hist_SODEs
RUNNAME=Nav10_RD_Avg1
NRUN=16

###########
# OPTIONS #
###########

# cross section
OPT1="--prob_hi_z 10."
#OPT1="--prob_hi_z 5."

# Schlogl (reaction-diffusion) or diffusion-only
OPT2="--nreactions 4"
#OPT2="--nreactions 0"

# averaging function "tilde n"(n1,n2)
# compare avg_type=1,2,3 (which mean)
# or avg_type=10,1,11,12 (which Heaviside)
OPT3="--avg_type 1"

# which random number generator 
OPTMISC1="--use_bl_rng F --seed 0"
#OPTMISC1="--use_bl_rng T --seed_diffusion 0 --seed_reaction 0 --seed_init 0"

# for testing with large dt
#OPTMISC2="--fixed_dt 0.1  --max_step 110000  --print_int 10000  --hydro_grid_int 10  --n_steps_skip 10000"
#OPTMISC2="--fixed_dt 0.01 --max_step 1100000 --print_int 100000 --hydro_grid_int 100 --n_steps_skip 100000"

OPTS="$OPT1 $OPT2 $OPT3 $OPT4 $OPTMISC1 $OPTMISC2"

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

