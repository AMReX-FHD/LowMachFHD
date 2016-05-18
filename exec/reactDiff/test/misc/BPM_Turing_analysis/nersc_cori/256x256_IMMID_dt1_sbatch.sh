#!/bin/bash -l

#SBATCH -p regular    
#SBATCH -N 1
#SBATCH -t 00:15:00 
#SBATCH -J IMMID256   

# for nersc/cori (queue: debug/regular)
# sbatch (script)
# squeue -u (id)
# scancel (job id)

RUNNAME=IMMID256
RUNNO=$1
INPUTSFILE=inputs-256x256-IMMID

#OPT1="--fixed_dt 1.  --max_step 10000  --print_int 10  --hydro_grid_int 0  --n_steps_write_avg 1  --plot_int 100"
OPT1="--fixed_dt 1.  --max_step 10000  --print_int 10  --hydro_grid_int 0  --n_steps_write_avg 1  --plot_int 0     --seed $RANDOM"

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
echo "srun -n 64 -e run_err -o scr_out ../$EXEC ../$INPUTSFILE $OPT1"
srun -n 64 -e run_err -o scr_out ../$EXEC ../$INPUTSFILE $OPT1 
