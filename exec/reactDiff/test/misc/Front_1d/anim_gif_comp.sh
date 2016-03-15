#!/bin/bash

#RUNNAME1=Schlogl_Dirichlet
#RUNNAME2=Schlogl_Dirichlet_det
#OUTDIR=Schlogl_Dirichlet_comp

RUNNAME1=Schlogl_PBC
RUNNAME2=Schlogl_PBC_det
OUTDIR=Schlogl_PBC_comp

GNUPLOTSCRIPT="Schlogl_comp.plt"

DATAFILEPREFIX=vstat

GNUPLOTOUTPUTEXT="png"
GIFDELAY=50
GIFOUTPUT=Schlogl_comp.gif

# check output directory 
if [ -d $RUNNAME1 ]
then
  cd $RUNNAME1
  DATAFILES=`ls $DATAFILEPREFIX????????`
  cd ..
else
  echo "directory $RUNNAME1 does not exist..."
  exit
fi

if [ ! -d $RUNNAME2 ]
then 
  echo "directory $RUNNAME2 does not exist..."
  exit
fi

if [ -d $OUTDIR ]
then 
  echo "directory $OUTDIR already exists..."
  sleep 3s
else
  mkdir $OUTDIR
fi

cd $OUTDIR

# make plots
for datafile in $DATAFILES 
do
  datafile1=../$RUNNAME1/$datafile
  datafile2=../$RUNNAME2/$datafile
  gnuplot -e "datafile='$datafile'" -e "datafile1='$datafile1'" -e "datafile2='$datafile2'" ../$GNUPLOTSCRIPT
  echo "$datafile.$GNUPLOTOUTPUTEXT generated"
done

# generate an animated gif
convert -delay $GIFDELAY $DATAFILEPREFIX*.$GNUPLOTOUTPUTEXT $GIFOUTPUT
echo "$GIFOUTPUT generated"
animate $GIFOUTPUT &
