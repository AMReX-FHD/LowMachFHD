#!/bin/bash -l

#SBATCH -p regular    
#SBATCH -N 1
#SBATCH -t 00:50:00 
#SBATCH -J EXMID256   

# for nersc/cori (queue: debug/regular)
# sbatch (script)
# squeue -u (id)
# scancel (job id)

RUNNAME=EXMID256
RUNNO=$1
INPUTSFILE=inputs-256x256-EXMID

#OPT1="--fixed_dt 0.025  --max_step 400000  --print_int 400  --hydro_grid_int 0  --n_steps_write_avg 40  --plot_int 4000"
OPT1="--fixed_dt 0.025  --max_step 400000  --print_int 400  --hydro_grid_int 0  --n_steps_write_avg 40  --plot_int 0     --seed $RANDOM"

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
