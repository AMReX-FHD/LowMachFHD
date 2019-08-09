#!/bin/bash 

####################
# slurm parameters #
####################

PARTITION=regular
NNODE=1
NPROC=16
TIMELIMIT=05:30:00
SLURMSCR=../slurm.sh

############
# run info #
############

EXEC=../main.Linux.gfortran.mpi.exe
INPUTSFILE=inputs_dimerization_giant_fluct

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

sbatch -p $PARTITION -N $NNODE -t $TIMELIMIT -J GF $SLURMSCR $NPROC $EXEC $INPUTSFILE $OPTS
