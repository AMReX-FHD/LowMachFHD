#!/bin/bash

# before running this script, you should run reactdiff program with "inputs_Schlogl_1d_analysis" and obtain "Schlogl.S_k.pair=1.Re.dat" in test directory.

DATA=../../Schlogl.S_k.pair=1.Re.dat
NEWNAME=res
cp $DATA $NEWNAME

GNUPLOTSCR=Sk.plt
gnuplot $GNUPLOTSCR
