#!/bin/bash -l

#SBATCH -p debug
#SBATCH -N 1 
#SBATCH -t 00:30:00 
#SBATCH -J AB2D 

# for nersc/cori (queue: debug/regular)
# sbatch (script)
# squeue -u (id)
# scancel (job id)

RUNNAME=RUN_stripe_2d
#RUNNAME=RUN_bubble_2d
RUNNO=$1
INPUTSFILE=inputs_Lemarchand_AB_stripe_2d
#INPUTSFILE=inputs_Lemarchand_AB_bubble_2d

EXEC=../../main.Linux.Intel.mpi.exe

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

cp $INPUTSFILE $RUNNAME

echo "cd $RUNNAME"
cd $RUNNAME
echo "srun -n 64 -e run_err -o scr_out ../$EXEC ../$INPUTSFILE"
srun -n 64 -e run_err -o scr_out ../$EXEC ../$INPUTSFILE  
