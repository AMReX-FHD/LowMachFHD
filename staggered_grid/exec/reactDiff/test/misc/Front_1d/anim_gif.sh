#!/bin/bash


EXEC=../../main.Linux.gfortran.debug.mpi.exe

RUNNAME=TEST

#INPUTSFILE=inputs_Schlogl_Dirichlet
#INPUTSFILE=inputs_Schlogl_Dirichlet_det
#INPUTSFILE=inputs_Schlogl_PBC
#INPUTSFILE=inputs_Schlogl_PBC_det

#GNUPLOTSCRIPT=Schlogl.plt
#GIFOUTPUT=Schlogl.gif

INPUTSFILE=inputs_AB_Dirichlet
#INPUTSFILE=inputs_AB_Dirichlet_det
#INPUTSFILE=inputs_AB_PBC
#INPUTSFILE=inputs_AB_PBC_det

GNUPLOTSCRIPT=AB.plt
GIFOUTPUT=AB.gif

DATAFILEPREFIX=vstat

GNUPLOTOUTPUTEXT="png"
GIFDELAY=50

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
  echo "try ./anim_gif.sh run"
  exit
fi

DATAFILES=`ls $DATAFILEPREFIX????????`

for datafile in $DATAFILES 
do
  gnuplot -e "datafile='$datafile'" ../$GNUPLOTSCRIPT
  echo "$datafile generated"
done

# generate an animated gif
convert -delay $GIFDELAY $DATAFILEPREFIX*.$GNUPLOTOUTPUTEXT $GIFOUTPUT
echo "$GIFOUTPUT generated"
animate $GIFOUTPUT &
