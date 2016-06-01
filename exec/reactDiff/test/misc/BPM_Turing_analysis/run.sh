#!/bin/bash

RUNNAME=64x64_ExMid_dt0.25
RUNNO=$1

### method ###
#INPUTSFILE=inputs_RDME
INPUTSFILE=inputs_ExMid
#INPUTSFILE=inputs_ImMid

### resolution ###
OPT1="--n_cells_x 64  --n_cells_y 64  --max_grid_size_x 64 --max_grid_size_y 64"
#OPT1="--n_cells_x 128 --n_cells_y 128 --max_grid_size_x 32 --max_grid_size_y 32"
#OPT1="--n_cells_x 256 --n_cells_y 256 --max_grid_size_x 32 --max_grid_size_y 32"

### time step ###
# for movie
#OPT2="--fixed_dt 0.025 --max_step 400000  --print_int 400  --hydro_grid_int 0  --n_steps_write_avg 40  --plot_int 4000"
#OPT2="--fixed_dt 0.05  --max_step 200000  --print_int 200  --hydro_grid_int 0  --n_steps_write_avg 20  --plot_int 2000"
#OPT2="--fixed_dt 0.1   --max_step 100000  --print_int 100  --hydro_grid_int 0  --n_steps_write_avg 10  --plot_int 1000"
#OPT2="--fixed_dt 0.2   --max_step  50000  --print_int  50  --hydro_grid_int 0  --n_steps_write_avg  5  --plot_int  500"
#OPT2="--fixed_dt 0.25  --max_step  40000  --print_int  40  --hydro_grid_int 0  --n_steps_write_avg  4  --plot_int  400"
#OPT2="--fixed_dt 0.5   --max_step  20000  --print_int  20  --hydro_grid_int 0  --n_steps_write_avg  2  --plot_int  200"
#OPT2="--fixed_dt 1.    --max_step  10000  --print_int  10  --hydro_grid_int 0  --n_steps_write_avg  1  --plot_int  100"
# for quantitative analysis
#OPT2="--fixed_dt 0.025 --max_step 300000  --print_int 400  --hydro_grid_int 0  --n_steps_write_avg 40  --plot_int 0"
#OPT2="--fixed_dt 0.05  --max_step 150000  --print_int 200  --hydro_grid_int 0  --n_steps_write_avg 20  --plot_int 0"
#OPT2="--fixed_dt 0.1   --max_step  75000  --print_int 100  --hydro_grid_int 0  --n_steps_write_avg 10  --plot_int 0"
#OPT2="--fixed_dt 0.2   --max_step  37500  --print_int  50  --hydro_grid_int 0  --n_steps_write_avg  5  --plot_int 0"
OPT2="--fixed_dt 0.25  --max_step  30000  --print_int  40  --hydro_grid_int 0  --n_steps_write_avg  4  --plot_int 0"
#OPT2="--fixed_dt 0.5   --max_step  15000  --print_int  20  --hydro_grid_int 0  --n_steps_write_avg  2  --plot_int 0"
#OPT2="--fixed_dt 1.    --max_step   7500  --print_int  10  --hydro_grid_int 0  --n_steps_write_avg  1  --plot_int 0"

### cross_section ###
#OPT3="--cross_section 0.5"
OPT3="--cross_section 1."
#OPT3="--cross_section 2."

OPTS="$OPT1 $OPT2 $OPT3"
EXEC=../../main.Linux.gfortran.mpi.exe
NPROC=1

if [ "$RUNNO" != "" ]
then
  RUNNAME=${RUNNAME}_RUN$RUNNO
fi

if [ -d $RUNNAME ]
then
  echo "directory $RUNNAME already exists..."
  sleep 3s
else
  mkdir $RUNNAME
fi

cp run.sh $RUNNAME
cd $RUNNAME

if [ "$RUNNO" = "" ]
then
  mpiexec -np $NPROC ../$EXEC ../$INPUTSFILE $OPTS | tee screen_out
else
  echo "mpiexec -np $NPROC ../$EXEC ../$INPUTSFILE $OPTS > screen_out &"
  mpiexec -np $NPROC ../$EXEC ../$INPUTSFILE $OPTS > screen_out &
fi
