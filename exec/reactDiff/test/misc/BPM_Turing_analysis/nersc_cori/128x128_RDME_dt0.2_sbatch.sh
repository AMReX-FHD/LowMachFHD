#!/bin/bash -l

#SBATCH -p regular    
#SBATCH -N 1
#SBATCH -t 00:8:00 
#SBATCH -J RDME128

# for nersc/cori (queue: debug/regular)
# sbatch (script)
# squeue -u (id)
# scancel (job id)

RUNNAME=RDME128
RUNNO=$1
INPUTSFILE=inputs-128x128-RDME

#OPT1="--fixed_dt 0.2  --max_step 50000  --print_int 50  --hydro_grid_int 0  --n_steps_write_avg 5  --plot_int 500"
OPT1="--fixed_dt 0.2  --max_step 50000  --print_int 50  --hydro_grid_int 0  --n_steps_write_avg 5  --plot_int 0     --seed $RANDOM"

EXEC=../main.Linux.Intel.mpi.exe

if [ "$RUNNO" != "" ]
then
  RUNNAME=${RUNNAME}_RUN$RUNNO
fi

if [ -d $RUNNAME ]
then
  echo "directory $RUNNAME already exists..."
else
  mkdir $RUNNAME
fi

echo "cd $RUNNAME"
cd $RUNNAME
echo "srun -n 16 -e run_err -o scr_out ../$EXEC ../$INPUTSFILE $OPT1"
srun -n 16 -e run_err -o scr_out ../$EXEC ../$INPUTSFILE $OPT1 
