#!/bin/bash

RUNNAME=TEST
EXEC=../../main.Linux.gfortran.mpi.exe
NPROC=4

INPUTS=inputs_ExMid
#INPUTS=inputs_RDME
#OPTS="--temporal_integrator 2 --reaction_type 1 --diffusion_type 1 --use_Poisson_rng 0"

OPTS_RESTART="--seed_diffusion -1 --seed_reaction -1"
#OPTS_RESTART="--seed_diffusion 0 --seed_reaction 0"

MAXSTEP=2000
CHKPT=1000

if [ -d $RUNNAME ]
then
  rm -rf $RUNNAME
fi

mkdir $RUNNAME
cd $RUNNAME

echo "##### 1st RUN #####"
mpiexec -n $NPROC ../$EXEC ../$INPUTS --max_step $MAXSTEP --chk_int $CHKPT --plot_int $MAXSTEP $OPTS

PLOTFILE=`printf "plt%08d" $MAXSTEP`
mv $PLOTFILE RUN1

echo "##### 2nd RUN #####"
mpiexec -n $NPROC ../$EXEC ../$INPUTS --max_step $MAXSTEP --chk_int 0 --restart $CHKPT --plot_int $MAXSTEP $OPTS $OPTS_RESTART

PLOTFILE=`printf "plt%08d" $MAXSTEP`
mv $PLOTFILE RUN2

DiffSameGrid22d infile1=RUN1 infile2=RUN2 
