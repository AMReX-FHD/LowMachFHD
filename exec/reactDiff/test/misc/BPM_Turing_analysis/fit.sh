#!/bin/bash

RUNNAME=TEST_RUN
NRUN=16
OUTPUT=$RUNNAME.fit

for ((i=1;i<=$NRUN;i++))
do
  cd $RUNNAME$i
  python ../fit.py
  tail -1 res.fit >> ../$OUTPUT
  cd ..
done
