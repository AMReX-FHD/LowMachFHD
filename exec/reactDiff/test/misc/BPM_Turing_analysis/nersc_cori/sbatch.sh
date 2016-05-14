#!/bin/bash -l

#SBATCH -p debug    
#SBATCH -N 1
#SBATCH -t 00:30:00 
#SBATCH -J test   

# for nersc/cori (queue: debug/regular)
# sbatch (script)
# squeue -u (id)
# scancel (job id)

RUNNAME=TEST
RUNNO=$1
INPUTSFILE=inputs-256x256-RDME
#INPUTSFILE=inputs-256x256-EXMID
#INPUTSFILE=inputs-256x256-IMMID

OPT1="--fixed_dt 0.03   --max_step 3000 --print_int 300 --hydro_grid_int -30 --n_steps_write_avg -30"
#OPT1="--fixed_dt 0.015  --max_step 680000 --print_int 600 --hydro_grid_int -60 --n_steps_write_avg -60"
#OPT1="--fixed_dt 0.03   --max_step 340000 --print_int 300 --hydro_grid_int -30 --n_steps_write_avg -30"
#OPT1="--fixed_dt 0.06   --max_step 170000 --print_int 150 --hydro_grid_int -15 --n_steps_write_avg -15"
#OPT1="--fixed_dt 0.09   --max_step 120000 --print_int 100 --hydro_grid_int -10 --n_steps_write_avg -10"
#OPT1="--fixed_dt 0.15   --max_step  68000 --print_int  60 --hydro_grid_int  -6 --n_steps_write_avg  -6"

EXEC=../../../main.Linux.Intel.mpi.exe

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
