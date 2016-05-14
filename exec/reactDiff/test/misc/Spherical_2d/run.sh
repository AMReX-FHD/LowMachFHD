#!/bin/bash

RUNNAME=DET

EXEC=../../main.Linux.gfortran.debug.mpi.exe
#INPUTS=inputs_Schlogl_2d
#INPUTS=inputs_Schlogl_2d_det
#INPUTS=inputs_Schlogl_2d_rate2
#INPUTS=inputs_Schlogl_2d_rate2_det
INPUTS=inputs_Schlogl_2d_rate2_det_model_file

mkdir $RUNNAME
cd $RUNNAME

mpiexec -np 4 ../$EXEC ../$INPUTS | tee screen_out
