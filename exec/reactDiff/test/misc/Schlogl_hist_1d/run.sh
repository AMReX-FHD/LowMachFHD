#!/bin/bash

###########
# OPTIONS #
###########

RUNNAME=TEST
NRUN=16
USEMPIEXEC=yes

#dV=5.
dV=10.
#dV=100.

OPT1="--nreactions 0"
#OPT1="--nreactions 4 --rate_multiplier 0.1"

#OPT2="--fixed_dt 1. --max_step 22000 --print_int 2000 --hydro_grid_int 1 --n_steps_skip 2000"
#OPT2="--fixed_dt 0.25 --max_step 44000 --print_int 4000 --hydro_grid_int 2 --n_steps_skip 4000"
OPT2="--fixed_dt 0.1 --max_step 110000 --print_int 10000 --hydro_grid_int 5 --n_steps_skip 10000"
#OPT2="--fixed_dt 0.01 --max_step 1100000 --print_int 100000 --hydro_grid_int 50 --n_steps_skip 100000"

OPT3="--temporal_integrator -2 --avg_type 1 --midpoint_stoch_flux_type 2"
#OPT3="--temporal_integrator -4 --avg_type 1 --midpoint_stoch_flux_type 2"
#OPT3="--temporal_integrator 0 --diffusion_type 3 --reaction_type 2"

OPT4="--use_Poisson_rng 1"

#OPT5="--initial_variance 1."
#OPT5="--initial_variance 1.  --integer_populsations T"
#OPT5="--initial_variance 0.  --integer_populations T"
OPT5="--initial_variance -1. --integer_populations F"

#######
# RUN #
#######

EXEC=../../main.Linux.gfortran.debug.mpi.exe
INPUTS=inputs_Schlogl_hist_1d
OPTS="--cross_section $dV $OPT1 $OPT2 $OPT3 $OPT4 $OPT5"
PYSCR_HIST=../../hist_n.py
PYSCR_SUMSK=sum_Sk_1d.py

if [ "$NRUN" == 1 ]
then
  if [ -d $RUNNAME ]; then echo "Warning: $RUNNAME already present"; else mkdir $RUNNAME; fi

  echo "cd $RUNNAME"
  cd $RUNNAME

  echo "../$EXEC ../$INPUTS $OPTS" | tee screen_out; sleep 3s
  ../$EXEC ../$INPUTS $OPTS | tee -a screen_out

  python ../$PYSCR_HIST fort.10 res.hist res.hist_cont res.hist_poiss 1. $dV yes
  python ../$PYSCR_SUMSK 
else
  for ((i=1;i<=$NRUN;i++))
  do 
    if [ -d $RUNNAME$i ]; then echo "Warning: $RUNNAME$i already present"; else mkdir $RUNNAME$i; fi

    echo "cd $RUNNAME$i"
    cd $RUNNAME$i

    if [ "$USEMPIEXEC" == yes ]
    then
      echo "mpiexec -n 1 ../$EXEC ../$INPUTS $OPTS >> screen_out &" | tee screen_out
      mpiexec -n 1 ../$EXEC ../$INPUTS $OPTS &>> screen_out &
    else
      echo "../$EXEC ../$INPUTS $OPTS" | tee screen_out
      ../$EXEC ../$INPUTS $OPTS | tee -a screen_out
    fi

    echo "cd -"
    cd -
    echo "" 
    sleep 3s
  done

  if [ "$USEMPIEXEC" == yes ]
  then 
    echo "waiting for all runs to finish..."
    wait
    echo ""
  fi

  for ((i=1;i<=$NRUN;i++))
  do
    echo "cd $RUNNAME$i"
    cd $RUNNAME$i
    python ../$PYSCR_HIST fort.10 res.hist res.hist_cont res.hist_poiss 1. $dV no
    python ../$PYSCR_SUMSK > res.sum_Sk
    echo "res.sum_Sk generated"
    cd - > /dev/null
  done
fi
