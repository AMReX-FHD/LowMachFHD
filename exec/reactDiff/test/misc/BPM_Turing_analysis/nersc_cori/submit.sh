#!/bin/bash

NRUN=16

date >> res.submit
for ((i=1;i<=$NRUN;i++))
do
  echo "sbatch sbatch.sh $i"
  sbatch sbatch.sh $i | tee -a res.submit
done 
