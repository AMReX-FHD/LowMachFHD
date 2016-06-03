#!/bin/bash -l

#SBATCH -p regular
#SBATCH -N 1
#SBATCH -t 00:10:00 
#SBATCH -J IMMID128  

# for nersc/cori (queue: debug/regular)
# sbatch (script)
# squeue -u (id)
# scancel (job id)

RUNNAME=128x128_ImMid_dt0.5
RUNNO=$1
INPUTSFILE=inputs_ImMid

OPT1="--n_cells_x 128 --n_cells_y 128 --max_grid_size_x 32 --max_grid_size_y 32"
#OPT2="--fixed_dt 0.5 --max_step 20000  --print_int 20  --hydro_grid_int 0  --n_steps_write_avg 2  --plot_int 200"
OPT2="--fixed_dt 0.5 --max_step 15000  --print_int 20  --hydro_grid_int 0  --n_steps_write_avg 2  --plot_int 0     --seed $RANDOM"
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
echo "srun -n 16 -e run_err -o scr_out ../$EXEC ../$INPUTSFILE $OPTS"
srun -n 16 -e run_err -o scr_out ../$EXEC ../$INPUTSFILE $OPTS 
