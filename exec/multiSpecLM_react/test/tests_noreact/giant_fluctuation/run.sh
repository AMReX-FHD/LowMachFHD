#!/bin/bash

EXEC=../../main.Linux.gfortran.mpi.exe
INPUTS=inputs_3species_diffbarrier_2d

OPT1="--max_step 10000  --print_int 10  --n_steps_skip 2500  --stats_int 1000"
#OPT1="--max_step 100000 --print_int 100 --n_steps_skip 25000 --stats_int 10000"

OPT2="--max_grid_size_x 32 --max_grid_size_y 16"

OPTS="$OPT1 $OPT2"

mpiexec -n 4 $EXEC $INPUTS $OPTS | tee scr_out
