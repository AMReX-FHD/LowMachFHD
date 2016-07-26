#!/bin/bash

RUNNAME=TEST
RUNNO=$1
EXEC=../../main.Linux.gfortran.mpi.exe
INPUTSFILE=inputs_BPM_Turing_2d
NPROC=4

# OPT0: cross_section
OPT0="--cross_section 1."
#OPT0="--cross_section 2."
#OPT0="--cross_section 10."

# OPT1: resolution
OPT1="--n_cells_x 64  --n_cells_y 64  --max_grid_size_x 32 --max_grid_size_y 32"
#OPT1="--n_cells_x 128 --n_cells_y 128 --max_grid_size_x 32 --max_grid_size_y 32"
#OPT1="--n_cells_x 256 --n_cells_y 256 --max_grid_size_x 32 --max_grid_size_y 32"

# OPT2: which numerical scheme
#OPT2="--temporal_integrator 2  --diffusion_type 1 --reaction_type 0 --use_Poisson_rng 2"        # CN/2+SSA+CN/2
#OPT2="--temporal_integrator 2  --diffusion_type 3 --reaction_type 0 --use_Poisson_rng 2"        # MN/2+SSA+MN/2
OPT2="--temporal_integrator -2 --use_Poisson_rng 1"    # unsplit explicit midpoint with tau-leaping 
#OPT2="--temporal_integrator -2 --use_Poisson_rng 2"    # unsplit explicit midpoint with SSA
#OPT2="--temporal_integrator -4 --use_Poisson_rng 1"    # unsplit implicit midpoint with tau-leaping 
#OPT2="--temporal_integrator -4 --use_Poisson_rng 2"    # unsplit implicit midpoint with SSA
#OPT2="--temporal_integrator -4 --use_Poisson_rng -1 --variance_coef_mass 0. --include_discrete_LMA_correction F"    # deterministic 

# OPT3: time step 
# for movie
#OPT3="--fixed_dt 0.025 --max_step 400000  --print_int 400  --n_steps_write_avg 40  --plot_int 4000"
#OPT3="--fixed_dt 0.05  --max_step 200000  --print_int 200  --n_steps_write_avg 20  --plot_int 2000"
#OPT3="--fixed_dt 0.1   --max_step 100000  --print_int 100  --n_steps_write_avg 10  --plot_int 1000"
#OPT3="--fixed_dt 0.2   --max_step  50000  --print_int  50  --n_steps_write_avg  5  --plot_int  500"
#OPT3="--fixed_dt 0.25  --max_step  40000  --print_int  40  --n_steps_write_avg  4  --plot_int  400"
OPT3="--fixed_dt 0.5   --max_step  20000  --print_int  20  --n_steps_write_avg  2  --plot_int  200"
#OPT3="--fixed_dt 1.    --max_step  10000  --print_int  10  --n_steps_write_avg  1  --plot_int  100"
#OPT3="--fixed_dt 2.    --max_step   5000  --print_int   5  --n_steps_write_avg  1  --plot_int   50"
# for quantitative analysis
#OPT3="--fixed_dt 0.025 --max_step 300000  --print_int 400  --n_steps_write_avg 40  --plot_int 0"
#OPT3="--fixed_dt 0.05  --max_step 150000  --print_int 200  --n_steps_write_avg 20  --plot_int 0"
#OPT3="--fixed_dt 0.1   --max_step  75000  --print_int 100  --n_steps_write_avg 10  --plot_int 0"
#OPT3="--fixed_dt 0.2   --max_step  37500  --print_int  50  --n_steps_write_avg  5  --plot_int 0"
#OPT3="--fixed_dt 0.25  --max_step  30000  --print_int  40  --n_steps_write_avg  4  --plot_int 0"
#OPT3="--fixed_dt 0.5   --max_step  15000  --print_int  20  --n_steps_write_avg  2  --plot_int 0"
#OPT3="--fixed_dt 1.    --max_step   7500  --print_int  10  --n_steps_write_avg  1  --plot_int 0"
#OPT3="--fixed_dt 2.    --max_step   3750  --print_int   5  --n_steps_write_avg  1  --plot_int 0"

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

if [ "$RUNNO" = "" ]
then
  echo "mpiexec -np $NPROC ../$EXEC ../$INPUTSFILE $OPTS" > scr_out
  echo "mpiexec -np $NPROC ../$EXEC ../$INPUTSFILE $OPTS | tee -a screen_out" 
  mpiexec -np $NPROC ../$EXEC ../$INPUTSFILE $OPTS | tee -a scr_out
else
  echo "mpiexec -np $NPROC ../$EXEC ../$INPUTSFILE $OPTS" > scr_out
  echo "mpiexec -np $NPROC ../$EXEC ../$INPUTSFILE $OPTS >> scr_out"
  mpiexec -np $NPROC ../$EXEC ../$INPUTSFILE $OPTS >> scr_out
fi
