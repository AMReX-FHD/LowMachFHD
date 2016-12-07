#!/bin/bash 

####################
# slurm parameters #
####################

PARTITION=debug
#PARTITION=regular
NNODE=3
NPROC=64
TIMELIMIT=00:30:00
#TIMELIMIT=01:00:00
SLURMSCR=slurm.sh

############
# run info #
############

RUNNAME=TEST
EXEC=../../main.Linux.Intel.mpi.exe
INPUTSFILE=inputs_Lemarchand_AB_bubble_3d

OPT1="--cross_section 1000."

OPT2="--initial_variance 1. --integer_populations T"       # noisy init cond (depending on cs)
#OPT2="--initial_variance 0. --integer_populations F"       # smooth init cond (no fluctuations)

OPT3="--variance_coef_mass 1. --use_Poisson_rng 1"         # diff fluct on, react fluct on
#OPT3="--variance_coef_mass 0. --use_Poisson_rng -1"        # diff fluct off, react fluct off

OPT4="--fixed_dt 0.25 --max_step 1000 --plot_int 100 --print_int 4"

OPTS="$OPT1 $OPT2 $OPT3 $OPT4"

if [ -d $RUNNAME ]
then
  echo "directory $RUNNAME already exists..."
  sleep 3s
else
  mkdir $RUNNAME
fi

cp $INPUTSFILE $RUNNAME
cp $0 $RUNNAME

echo "cd $RUNNAME"
cd $RUNNAME

sbatch -p $PARTITION -N $NNODE -t $TIMELIMIT -J $RUNNAME ../$SLURMSCR $NPROC ../$EXEC $INPUTSFILE $OPTS
