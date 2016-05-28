#!/bin/bash -l

#SBATCH -p debug
#SBATCH -N 8 
#SBATCH -t 00:30:00 
#SBATCH -J TEST

# for nersc/cori (queue: debug/regular)
# sbatch (script)
# squeue -u (id)
# scancel (job id)

RUNNAME=TEST
RUNNO=$1
INPUTSFILE=inputs_Lemarchand_AB_bubble_3d

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

cp $INPUTSFILE $RUNNAME
cp sbatch.sh $RUNNAME

echo "cd $RUNNAME"
cd $RUNNAME
echo "srun -n 512 -e run_err -o scr_out ../$EXEC ../$INPUTSFILE"
srun -n 512 -e run_err -o scr_out ../$EXEC ../$INPUTSFILE  
