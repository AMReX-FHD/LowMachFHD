#!/bin/bash -l

#SBATCH -p regular    
#SBATCH -N 1
#SBATCH -t 00:45:00 
#SBATCH -J EXMID256   

# for nersc/cori (queue: debug/regular)
# sbatch (script)
# squeue -u (id)
# scancel (job id)

RUNNAME=256x256_ExMid_dt0.025
RUNNO=$1
INPUTSFILE=inputs_ExMid

OPT1="--n_cells_x 256 --n_cells_y 256 --max_grid_size_x 32 --max_grid_size_y 32"
#OPT2="--fixed_dt 0.025 --max_step 400000  --print_int 400  --hydro_grid_int 0  --n_steps_write_avg 40  --plot_int 4000"
OPT2="--fixed_dt 0.025 --max_step 300000  --print_int 400  --hydro_grid_int 0  --n_steps_write_avg 40  --plot_int 0     --seed $RANDOM"
OPT3="--cross_section 2."
OPTS="$OPT1 $OPT2 $OPT3"

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
echo "srun -n 64 -e run_err -o scr_out ../$EXEC ../$INPUTSFILE $OPTS"
srun -n 64 -e run_err -o scr_out ../$EXEC ../$INPUTSFILE $OPTS 
