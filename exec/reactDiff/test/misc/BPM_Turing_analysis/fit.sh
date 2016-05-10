#!/bin/bash

RUNNAME=TEST_RUN
NRUN=16

for ((i=1;i<=$NRUN;i++))
do
  cd $RUNNAME$i
  python ../fit.py
  cd ..
done
