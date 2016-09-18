#!/bin/bash

RUNNAME=TEST

EXEC=../../main.Linux.gfortran.mpi.exe
INPUTSFILE=inputs_Schlogl_monomodal_1d
PYSCR=Sk_xi.py

# OPT1: timestep and sampling 
MAXSTEP=60000
#MAXSTEP=510000
#MAXSTEP=1010000
NSTEPSSKIP=10000
HYDROGRIDINT=1
OPT1="--max_step $MAXSTEP --n_steps_skip $NSTEPSSKIP --hydro_grid_int $HYDROGRIDINT"

# OPT2: which RNG  
OPT2="--use_bl_rng F --seed 0"
#OPT2="--use_bl_rng T --seed_diffusion 0 --seed_reaction 0 --seed_init 0"

# OPT3: which numerical scheme
OPT3="--temporal_integrator -2 --use_Poisson_rng 1"	# ExMidTau
#OPT3="--temporal_integrator -2 --use_Poisson_rng 2"	# ExMidSSA
#OPT3="--temporal_integrator 2 --diffusion_type 3 --reaction_type 1 --use_Poisson_rng 1"	# MN/2+Tau2nd+MN/2
#OPT3="--temporal_integrator 2 --diffusion_type 3 --reaction_type 0 --use_Poisson_rng 2"	# MN/2+SSA+MN/2
#OPT3="--temporal_integrator -4 --use_Poisson_rng 1"	# ImMidTau
#OPT3="--temporal_integrator -4 --use_Poisson_rng 2"	# ImMidSSA
#OPT3="--temporal_integrator 2 --diffusion_type 1 --reaction_type 1 --use_Poisson_rng 1"	# CN/2+Tau2nd+CN/2
#OPT3="--temporal_integrator 2 --diffusion_type 1 --reaction_type 0 --use_Poisson_rng 2"	# CN/2+SSA+CN/2

# OPT4: alpha and beta
RM=0.0285714	# alpha=0.05
#RM=0.0571429	# alpha=0.1
#RM=0.285714	# alpha=0.5 
#DFICK=0.1	# beta=0.1
#DFICK=0.2	# beta=0.2
DFICK=0.25	# beta=0.25
#DFICK=2.5	# beta=2.5
#DFICK=5.	# beta=5
#DFICK=10.	# beta=10
OPT4="--rate_multiplier $RM --D_Fick_1 $DFICK"

OPTS="$OPT1 $OPT2 $OPT3 $OPT4"

if [ ! -f $EXEC ]
then
  echo "executable $EXEC does not exist..."
  exit
fi

if [ -d $RUNNAME ]
then
  echo "directory $RUNNAME already exists..."
  sleep 3s
else
  echo "directory $RUNNAME has been created."
  mkdir $RUNNAME 
fi

cd $RUNNAME

cp ../$INPUTSFILE .
cp ../run.sh .

echo "mpiexec -n 1 ../$EXEC ../$INPUTSFILE $OPTS" > scr_out
echo "mpiexec -n 1 ../$EXEC ../$INPUTSFILE $OPTS | tee -a scr_out"
mpiexec -n 1 ../$EXEC ../$INPUTSFILE $OPTS | tee -a scr_out

python ../$PYSCR $RM $DFICK
