#!/bin/bash

RUNNAME=cs10_Nc64_ExMidTau_dt0.5
NRUN=4
OUTPUT=$RUNNAME.fit

for ((i=1;i<=$NRUN;i++))
do
  cd ${RUNNAME}_RUN$i
  python ../fit.py
  tail -1 res.fit >> ../$OUTPUT
  cd ..
done
