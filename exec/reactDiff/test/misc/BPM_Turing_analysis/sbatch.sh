#!/bin/bash -l

# for nersc/cori (queue: debug/regular)
# sbatch (script)
# squeue -u (id)
# scancel (job id)

#SBATCH -p debug    
#SBATCH -N 2       
#SBATCH -t 00:30:00 
#SBATCH -J test   

RUNNAME=TEST
NRUN=16
INPUTSFILE=inputs-256x256-RDME
#INPUTSFILE=inputs-256x256-EXMID
#INPUTSFILE=inputs-256x256-IMMID
EXEC=../../main.Linux.Intel.mpi.exe

OPT1="--fixed_dt 0.03   --max_step 3000 --print_int 300 --hydro_grid_int -30 --n_steps_write_avg -30"
#OPT1="--fixed_dt 0.015  --max_step 680000 --print_int 600 --hydro_grid_int -60 --n_steps_write_avg -60"
#OPT1="--fixed_dt 0.03   --max_step 340000 --print_int 300 --hydro_grid_int -30 --n_steps_write_avg -30"
#OPT1="--fixed_dt 0.06   --max_step 170000 --print_int 150 --hydro_grid_int -15 --n_steps_write_avg -15"
#OPT1="--fixed_dt 0.09   --max_step 120000 --print_int 100 --hydro_grid_int -10 --n_steps_write_avg -10"
#OPT1="--fixed_dt 0.15   --max_step  68000 --print_int  60 --hydro_grid_int  -6 --n_steps_write_avg  -6"

mkdir $RUNNAME
cd $RUNNAME

for ((i=1;i<=NRUN;i++))
do
  mkdir ${RUNNAME}_RUN$i
  cd ${RUNNAME}_RUN$i
  srun -n 4 ../../$EXEC ../../$INPUTSFILE $OPT1 > screen_out &
  sleep 5s
  cd ..
done

wait
