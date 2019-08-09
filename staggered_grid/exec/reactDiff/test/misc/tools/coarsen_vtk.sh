#!/bin/bash

# compile coarsen_vtk/coarsen.f90 to obtain an excutable
# (e.g. coarsen_vtk/coarsen_vtk_256_4) 

RUNNAME=TEST256
RUNNAME2=TEST256_factor4

COARSEN_VTK_DIR=`pwd`/coarsen_vtk
COARSEN_VTK_SCR=$COARSEN_VTK_DIR/coarsen_vtk_256_4_single.sh

if [ -d $RUNNAME2 ]
then
  echo "$RUNNAME2 already exists..."
  exit
else
  mkdir $RUNNAME2
fi

cd $RUNNAME
VTKFILES=`ls *.vtk`
cd ..

for vtk in $VTKFILES
do
  echo "$COARSEN_VTK_SCR $RUNNAME/$vtk $RUNNAME2/$vtk $COARSEN_VTK_DIR"
  $COARSEN_VTK_SCR $RUNNAME/$vtk $RUNNAME2/$vtk $COARSEN_VTK_DIR
done
