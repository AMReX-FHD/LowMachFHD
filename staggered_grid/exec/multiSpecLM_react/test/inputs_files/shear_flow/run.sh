#/bin/bash

EXEC=../../main.Linux.gfortran.mpi.exe
INPUTS=../inputs_shear_flow
RUNNAME=$1

#######
# RUN #
#######

if [ "$RUNNAME" == "" ]
then
  RUNNAME=TEST
fi

echo $RUNNAME

# check executable and inputs file

if [ ! -f $EXEC ]
then
  echo "ERROR: $EXEC does not exist"
  exit
fi

if [ ! -f $INPUTS ]
then
  echo "ERROR: $INPUTS does not exist"
  exit
fi

if [ -d $RUNNAME ]
then
  echo "ERROR: $RUNNAME already exists"
  exit
fi

mkdir $RUNNAME
cp $INPUTS $RUNNAME
cp $0 $RUNNAME

cd $RUNNAME
mpiexec -n 4 ../$EXEC ../$INPUTS $OPTS | tee scr_out
