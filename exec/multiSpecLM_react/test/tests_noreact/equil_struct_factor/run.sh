#!/bin/bash

EXEC=../../main.Linux.gfortran.mpi.exe
INPUTS=inputs_3species_structure_factors_2d

OPTS="--max_step 12000  --print_int 100  --n_steps_skip 2000"
#OPTS="--max_step 120000 --print_int 1000 --n_steps_skip 20000"

mpiexec -n 4 $EXEC $INPUTS $OPTS | tee scr_out
