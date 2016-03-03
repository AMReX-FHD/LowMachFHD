#!/bin/bash

RUNNAME=TEST

EXEC=../../main.Linux.gfortran.debug.mpi.exe
INPUTSFILE=inputs_AB
GNUPLOTSCRIPT="AB_vstat.plt"

DATAFILEPREFIX=vstat
TIMEMIN=0
TIMEINCR=1000
TIMEMAX=20000

GNUPLOTOUTPUTEXT="png"
GIFDELAY=50
GIFOUTPUT=AB.gif

# if the script is executed by "./anim_gif.sh run", run the reactdiff program
# otherwise, skip this 
if [ $# -eq 1 ]
then
  if [ $1 == run ]
  then
    if [ -d $RUNNAME ]
    then
      echo "** directory $RUNNAME already present **"
      sleep 3s
    else
      mkdir $RUNNAME
    fi

    cd $RUNNAME
    ../$EXEC ../$INPUTSFILE
    cd ..
  fi
fi

# generate figures from datafiles
if [ -d $RUNNAME ]
then 
  cd  $RUNNAME
else
  echo "directory $RUNNAME does not exist"
  exit
fi

for i in $(seq -f "%08g" $TIMEMIN $TIMEINCR $TIMEMAX)
do
  datafile=$DATAFILEPREFIX$i
  gnuplot -e "datafile='$datafile'" ../$GNUPLOTSCRIPT
  echo "$datafile generated"
done

# generate an animated gif
convert -delay $GIFDELAY $DATAFILEPREFIX*.$GNUPLOTOUTPUTEXT $GIFOUTPUT
echo "$GIFOUTPUT generated"
animate $GIFOUTPUT &
