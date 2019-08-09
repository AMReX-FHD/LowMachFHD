#!/bin/bash 

####################
# slurm parameters #
####################

PARTITION=debug
NNODE=11
NPROC=256
TIMELIMIT=0:30:00
SLURMSCR=../slurm.sh

############
# run info #
############

EXEC=../main.Linux.gfortran.mpi.exe
INPUTSFILE=inputs_abfinger_ions_3d

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

sbatch -p $PARTITION -N $NNODE -t $TIMELIMIT -J abfinger $SLURMSCR $NPROC $EXEC $INPUTSFILE $OPTS
