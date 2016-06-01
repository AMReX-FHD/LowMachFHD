#!/bin/bash -l

#SBATCH -p debug
#SBATCH -N 8 
#SBATCH -t 00:30:00 
#SBATCH -J AB3D

# for nersc/cori (queue: debug/regular)
# sbatch (script)
# squeue -u (id)
# scancel (job id)

RUNNAME=RUN_stripe_3d
#RUNNAME=RUN_bubble_3d
RUNNO=$1
INPUTSFILE=inputs_Lemarchand_AB_stripe_3d
#INPUTSFILE=inputs_Lemarchand_AB_bubble_3d

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
echo "srun -n 512 -e run_err -o scr_out ../$EXEC ../$INPUTSFILE"
srun -n 512 -e run_err -o scr_out ../$EXEC ../$INPUTSFILE  
