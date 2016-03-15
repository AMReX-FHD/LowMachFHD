#!/bin/bash

# run reactdiff program with inputs file
EXEC=../../main.Linux.gfortran.debug.mpi.exe
INPUTSFILE=inputs_Schlogl_linearized_eqn_analysis_1d
$EXEC $INPUTSFILE

# generate a figure by using gnuplot and show it
GNUPLOTSCR=Sk.plt
GNUPLOTOUTPUT=Sk.eps
gnuplot $GNUPLOTSCR
gv $GNUPLOTOUTPUT &
