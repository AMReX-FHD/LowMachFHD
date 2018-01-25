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

# initial conditions
OPT1="--initial_variance_mass 1. --integer_populations T"   # noisy init cond
#OPT1="--initial_variance_mass 0. --integer_populations F"  # smooth init cond (no fluctuations)

# stochastic diffusion and reaction
OPT2="--variance_coef_mass 1. --use_Poisson_rng 1"          # diff fluct on, react fluct on
#OPT2="--variance_coef_mass 0. --use_Poisson_rng -1"        # diff fluct off, react fluct off

# timestep
OPT3="--fixed_dt 0.25 --max_step 1000 --plot_int 20 --print_int 4"

# which random number generator 
OPTMISC1="--use_bl_rng F --seed 0"
#OPTMISC1="--use_bl_rng T --seed_diffusion 0 --seed_reaction 0 --seed_init_mass 0"

# max_grid_size
OPTMISC2="--max_grid_size_x 64 --max_grid_size_y 64 --max_grid_size_z 64"
#OPTMISC2="--max_grid_size_x 32 --max_grid_size_y 32 --max_grid_size_z 32"

OPTS="$OPT1 $OPT2 $OPT3 $OPTMISC1 $OPTMISC2"

if [ -d $RUNNAME ]
then
  echo "directory $RUNNAME already exists..."
  exit
fi

mkdir $RUNNAME

cp $INPUTSFILE $RUNNAME
cp $0 $RUNNAME

echo "cd $RUNNAME"
cd $RUNNAME

sbatch -p $PARTITION -N $NNODE -t $TIMELIMIT -J $RUNNAME ../$SLURMSCR $NPROC ../$EXEC $INPUTSFILE $OPTS
