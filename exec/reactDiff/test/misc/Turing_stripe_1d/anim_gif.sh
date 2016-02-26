#!/bin/bash

EXEC=../../main.Linux.gfortran.debug.mpi.exe
INPUTSFILE=inputs_Turing_stripe_1d

DATAFILEPREFIX=vstat
TIMEMIN=0
TIMEINCR=1000
TIMEMAX=20000

GNUPLOTSCRIPT="vstat.plt"
GNUPLOTOUTPUTEXT="png"
GIFDELAY=50
GIFOUTPUT=${DATAFILEPREFIX}.gif

# if the script is executed by "./anim_gif.sh run", run the reactdiff program
# otherwise, skip this 
if [ $# -eq 1 ]
then
  if [ $1 == run ]
  then
  $EXEC $INPUTSFILE
  fi
fi

# generate figures from datafiles
for i in $(seq -f "%08g" $TIMEMIN $TIMEINCR $TIMEMAX)
do
  datafile=$DATAFILEPREFIX$i
  gnuplot -e "datafile='$datafile'" $GNUPLOTSCRIPT
  echo "$datafile generated"
done

# generate an animated gif
convert -delay $GIFDELAY $DATAFILEPREFIX*.$GNUPLOTOUTPUTEXT $GIFOUTPUT
echo "$GIFOUTPUT generated"
animate $GIFOUTPUT &
