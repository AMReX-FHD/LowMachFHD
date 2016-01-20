#!/bin/bash

###########
# OPTIONS #
###########

RUNNAME=TEST
NRUN=8

dV=5.

OPT0="--nreactions 0"
#OPT0="--nreactions 4 --rate_multiplier 0.1"

OPT1="--fixed_dt 0.1 --max_step 50000 --print_int 10000 --hydro_grid_int 5 --n_steps_skip 10000"

OPT2="--temporal_integrator -2"
#OPT2="--temporal_integrator 0 --diffusion_type 3 --reaction_type 0"

OPT3="--use_Poisson_rng 1 --avg_type 3"
OPT4="--midpoint_stoch_flux_type 2"

#######
# RUN #
#######

EXEC=../../main.Linux.gfortran.debug.mpi.exe
INPUTS=inputs_Schlogl_hist_1d
OPTS="--cross_section $dV $OPT0 $OPT1 $OPT2 $OPT3 $OPT4"

if [ $NRUN == 1 ]
then
  if [ -d $RUNNAME ]; then echo "Warning: $RUNNAME already present"; else mkdir $RUNNAME; fi

  echo "cd $RUNNAME"
  cd $RUNNAME

  echo "../$EXEC ../$INPUTS $OPTS" | tee screen_out; sleep 3s
  ../$EXEC ../$INPUTS $OPTS | tee -a screen_out

  PYSCR_HIST=../../hist_n.py
  python ../$PYSCR_HIST fort.10 res.hist res.hist_cont res.hist_poiss 1. $dV
else
  for ((i=1;i<=$NRUN;i++))
  do 
    if [ -d $RUNNAME$i ]; then echo "Warning: $RUNNAME$i already present"; else mkdir $RUNNAME$i; fi

    echo "cd $RUNNAME$i"
    cd $RUNNAME$i

    echo "../$EXEC ../$INPUTS $OPTS" | tee screen_out; sleep 3s
    ../$EXEC ../$INPUTS $OPTS | tee -a screen_out

    PYSCR_HIST=../../hist_n.py
    python ../$PYSCR_HIST fort.10 res.hist res.hist_cont res.hist_poiss 1. $dV no

    cd -
  done
fi
