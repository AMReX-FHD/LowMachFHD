#!/bin/bash

RUNNAME=TEST
NRUN=16

PYSCR1=cond_exp.py
PYSCR2=cond_exp_stat.py

################################################################################

for ((i=1;i<=NRUN;i++))
do
  dir=$RUNNAME$i
  cd $dir
  echo "***** $dir"
  python ../$PYSCR1 $NCELL $DV $AVGTYPE | tee res.face_avg
  cd ..
  echo 
done

python $PYSCR2 $RUNNAME $NRUN
