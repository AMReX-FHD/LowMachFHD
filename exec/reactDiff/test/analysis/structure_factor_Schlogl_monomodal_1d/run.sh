#!/bin/bash

RUNNAME=ExMid_CLE
EXEC=../../main.Linux.gfortran.mpi.exe
INPUTSFILE=inputs_Schlogl_monomodal_1d
PYSCR=Sk_xi.py

# timestep and sampling 
MAXSTEP=60000
NSTEPSSKIP=10000
HYDROGRIDINT=1
OPT1="--max_step $MAXSTEP --n_steps_skip $NSTEPSSKIP --hydro_grid_int $HYDROGRIDINT"

# which RNG  
OPT2="--use_bl_rng F --seed 0"
#OPT2="--use_bl_rng T --seed_diffusion 0 --seed_reaction 0 --seed_init 0"

# which numerical scheme
OPT3="--temporal_integrator -2 --use_Poisson_rng 0"	# unsplit explicit midpoint with CLE
#OPT3="--temporal_integrator -2 --use_Poisson_rng 2"	# unsplit explicit midpoint with SSA

# alpha and beta
RM=0.0285714	# alpha=0.05
#RM=0.0571429	# alpha=0.1
#RM=0.285714	# alpha=0.5 
#DFICK=0.1	# beta=0.1
DFICK=0.2	# beta=0.2
#DFICK=2.5	# beta=2.5
#DFICK=5.	# beta=5
#DFICK=10.	# beta=10
OPT4="--rate_multiplier $RM --D_Fick_1 $DFICK"

OPTS="$OPT1 $OPT2 $OPT3 $OPT4"

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

echo "../$EXEC ../$INPUTSFILE $OPTS | tee scr_out"
../$EXEC ../$INPUTSFILE $OPTS | tee scr_out

python ../$PYSCR $RM $DFICK
