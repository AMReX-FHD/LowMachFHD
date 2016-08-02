#!/bin/bash

RUNNAME=$1
NRUN=$2
OUTPUT=$RUNNAME.fit

for ((i=1;i<=$NRUN;i++))
do
  cd ${RUNNAME}_RUN$i
  python ../fit.py
  tail -1 res.fit >> ../$OUTPUT
  cd ..
done
