#!/bin/bash

###########
# OPTIONS #
###########

RUNNAME=TEST

#dV=5.
dV=10.
#dV=100.
#dV=1000.

#OPT1="--nreactions 0"
OPT1="--nreactions 4 --rate_multiplier 0.1"

#OPT2="--fixed_dt 2.   --max_step 5500   --print_int 500    --hydro_grid_int 1   --n_steps_skip 500"
#OPT2="--fixed_dt 1.   --max_step 6000   --print_int 500    --hydro_grid_int 1   --n_steps_skip 1000"
OPT2="--fixed_dt 0.5  --max_step 12000  --print_int 1000   --hydro_grid_int 2   --n_steps_skip 2000"
#OPT2="--fixed_dt 0.25 --max_step 24000  --print_int 2000   --hydro_grid_int 4   --n_steps_skip 4000"
#OPT2="--fixed_dt 0.1  --max_step 60000  --print_int 5000   --hydro_grid_int 10  --n_steps_skip 10000"

#OPT3="--temporal_integrator -2 --avg_type 1 --midpoint_stoch_flux_type 3"
OPT3="--temporal_integrator -4 --avg_type 1 --midpoint_stoch_flux_type 3"
#OPT3="--temporal_integrator 0 --diffusion_type 3 --reaction_type 1"

OPT4="--use_Poisson_rng 1"

OPT5="--initial_variance 1."
#OPT5="--initial_variance 1.  --integer_populsations T"
#OPT5="--initial_variance 0.  --integer_populations T"
#OPT5="--initial_variance -1. --integer_populations F"

#######
# RUN #
#######

EXEC=../../main.Linux.gfortran.debug.mpi.exe
INPUTS=inputs_Schlogl_hist_2d
OPTS="--cross_section $dV $OPT1 $OPT2 $OPT3 $OPT4 $OPT5"

if [ -d $RUNNAME ]; then echo "Warning: $RUNNAME already present"; else mkdir $RUNNAME; fi

echo "cd $RUNNAME"
cd $RUNNAME

echo "mpiexec - 4 ../$EXEC ../$INPUTS $OPTS" | tee screen_out; sleep 3s
mpiexec -n 4 ../$EXEC ../$INPUTS $OPTS | tee -a screen_out

cd ..

########
# hist #
########

SCRHISTPY=../../hist_n.py

cd $RUNNAME
python ../$SCRHISTPY 1. $dV yes yes
cd ..

########
# S(k) #
########

SCRSKVISIT=visit_drawSk2d.py

cd $RUNNAME
visit -nowin -cli -s ../$SCRSKVISIT
eog Sk.png &
cd ..
