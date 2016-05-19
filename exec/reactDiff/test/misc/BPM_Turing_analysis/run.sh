#!/bin/bash

RUNNAME=TEST
RUNNO=
INPUTSFILE=inputs-256x256-RDME
#INPUTSFILE=inputs-256x256-EXMID
#INPUTSFILE=inputs-256x256-IMMID

#OPT1="--fixed_dt 0.025 --max_step 400000  --print_int 400  --hydro_grid_int 0  --n_steps_write_avg 40  --plot_int 0"
OPT1="--fixed_dt 0.05  --max_step 200000  --print_int 200  --hydro_grid_int 0  --n_steps_write_avg 20  --plot_int 0"
#OPT1="--fixed_dt 0.1   --max_step 100000  --print_int 100  --hydro_grid_int 0  --n_steps_write_avg 10  --plot_int 0"
#OPT1="--fixed_dt 0.2   --max_step  50000  --print_int  50  --hydro_grid_int 0  --n_steps_write_avg  5  --plot_int 0"
#OPT1="--fixed_dt 0.25  --max_step  40000  --print_int  40  --hydro_grid_int 0  --n_steps_write_avg  4  --plot_int 0"
#OPT1="--fixed_dt 0.5   --max_step  20000  --print_int  20  --hydro_grid_int 0  --n_steps_write_avg  2  --plot_int 0"
#OPT1="--fixed_dt 1.    --max_step  20000  --print_int  10  --hydro_grid_int 0  --n_steps_write_avg  1  --plot_int 0"

EXEC=../../main.Linux.gfortran.debug.mpi.exe
NPROC=4

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

cd $RUNNAME

if [ "$RUNNO" = "" ]
then
  mpiexec -np $NPROC ../$EXEC ../$INPUTSFILE $OPT1 | tee screen_out
else
  echo "mpiexec -np $NPROC ../$EXEC ../$INPUTSFILE $OPT1 > screen_out &"
  mpiexec -np $NPROC ../$EXEC ../$INPUTSFILE $OPT1 > screen_out &
fi
