#!/bin/bash

# in order to execute this script, you need to compile fcoarsen:
# BoxLib/Tools/Postprocessing/F_Src/fcoarsen.f90

EXEC=fcoarsen
FACTOR=2
RUNNAME=TEST
RUNNAME2=${RUNNAME}_factor$FACTOR

if [ -d $RUNNAME2 ]
then
  echo "$RUNNAME2 already exists..."
  exit
else
  mkdir $RUNNAME2
fi

cd $RUNNAME
PLOTFILES=`ls -d plt*`
cd ..

for pf in $PLOTFILES
do
  echo "$EXEC --coarsen $FACTOR --infile $RUNNAME/$pf --outfile $RUNNAME2/$pf"
  $EXEC --coarsen $FACTOR --infile $RUNNAME/$pf --outfile $RUNNAME2/$pf
done
