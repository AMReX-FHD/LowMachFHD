#!/bin/bash 

####################
# slurm parameters #
####################

PARTITION=regular
NNODE=3
NPROC=64
TIMELIMIT=3:30:00
SLURMSCR=../slurm.sh

############
# run info #
############

EXEC=../main.Linux.gfortran.mpi.exe
INPUTSFILE=inputs_asymm_convect_fingers

OPTS=

#########
# check #
#########

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

sbatch -p $PARTITION -N $NNODE -t $TIMELIMIT -J fingerDD $SLURMSCR $NPROC $EXEC $INPUTSFILE $OPTS
