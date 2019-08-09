#!/bin/bash 

####################
# slurm parameters #
####################

PARTITION=regular
NNODE=1
NPROC=24
TIMELIMIT=1:00:00

###################
# norm parameters #
###################

SLURMSCR=./norm.sh
RUNNAME=full

sbatch -p $PARTITION -N $NNODE -t $TIMELIMIT -J norm $SLURMSCR $NPROC $RUNNAME
