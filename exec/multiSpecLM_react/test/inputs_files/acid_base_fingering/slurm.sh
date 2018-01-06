#!/bin/bash -l

cnt=1
for arg in $*
do
  if [ $cnt == 1 ]
  then
    NPROC=$arg
  elif [ $cnt == 2 ]
  then
    EXEC=$arg
  elif [ $cnt  == 3 ]
  then
    INPUTSFILE=$arg
  else
    OPTS="$OPTS $arg"
  fi 
  ((cnt=cnt+1))
done

echo "srun -n $NPROC $EXEC $INPUTSFILE $OPTS"
srun -n $NPROC $EXEC $INPUTSFILE $OPTS
