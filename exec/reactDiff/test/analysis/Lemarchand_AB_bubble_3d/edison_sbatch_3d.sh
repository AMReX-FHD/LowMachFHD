#!/bin/bash -l

#SBATCH -p regular 
#SBATCH -N 3 
#SBATCH -t 01:00:00 
#SBATCH -J bubble3D

# for nersc/edison (queue: debug/regular)
# sbatch (this script)
# sqs
# scancel (job id)

RUNNAME=RUN_bubble_3d
INPUTSFILE=inputs_Lemarchand_AB_bubble_3d

EXEC=../../main.Linux.Intel.mpi.exe

if [ -d $RUNNAME ]
then
  echo "directory $RUNNAME already exists..."
  sleep 3s
else
  mkdir $RUNNAME
fi

cp $INPUTSFILE $RUNNAME
cp $0 $RUNNAME

echo "cd $RUNNAME"
cd $RUNNAME
srun -n 72 ../$EXEC ../$INPUTSFILE  
