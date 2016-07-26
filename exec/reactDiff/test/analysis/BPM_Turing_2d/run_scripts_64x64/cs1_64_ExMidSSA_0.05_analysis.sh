#!/bin/bash

RUNNAME=cs1_64_ExMidSSA_0.05
RUNNO=$1
EXEC=../../main.Linux.gfortran.mpi.exe
INPUTSFILE=inputs_BPM_Turing_2d
NPROC=1

# OPT0: cross_section
OPT0="--cross_section 1."

# OPT1: resolution
OPT1="--n_cells_x 64  --n_cells_y 64  --max_grid_size_x 32 --max_grid_size_y 32"

# OPT2: which numerical scheme
OPT2="--temporal_integrator -2 --use_Poisson_rng 2"    # unsplit explicit midpoint with SSA

# OPT3: time step 
# for quantitative analysis
OPT3="--fixed_dt 0.05  --max_step 150000  --print_int 200  --n_steps_write_avg 20  --plot_int 0"

# OPT4: which RNG  
OPT4="--use_bl_rng F --seed 0"
#OPT4="--use_bl_rng T --seed_diffusion 0 --seed_reaction 0 --seed_init 0"

OPTS="$OPT0 $OPT1 $OPT2 $OPT3 $OPT4"

if [ ! -f $EXEC ]
then
  echo "executable $EXEC does not exist..."
  exit
fi

if [ "$RUNNO" = "plot" ]
then
  RUNNAME=${RUNNAME}_plot
  RUNNO=
elif [ "$RUNNO" != "" ]
then 
  RUNNAME=${RUNNAME}_RUN$RUNNO
fi

if [ -d $RUNNAME ]
then
  echo "directory $RUNNAME already exists..."
  sleep 3s
else
  echo "directory $RUNNAME has been created."
  mkdir $RUNNAME
fi

cp $0 $RUNNAME
cp $INPUTSFILE $RUNNAME

cd $RUNNAME

echo "mpiexec -np $NPROC ../$EXEC ../$INPUTSFILE $OPTS" > scr_out
echo "mpiexec -np $NPROC ../$EXEC ../$INPUTSFILE $OPTS >> scr_out"
mpiexec -np $NPROC ../$EXEC ../$INPUTSFILE $OPTS >> scr_out
