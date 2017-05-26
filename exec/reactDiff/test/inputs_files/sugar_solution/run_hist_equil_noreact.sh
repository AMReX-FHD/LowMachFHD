#!/bin/bash

EXEC=../../main.Linux.gfortran.mpi.exe
INPUTS=../inputs_sugar_solution
RUNNAME=TEST
NRUN=4

###########
# OPTIONS #
###########

# assume no reaction
OPT0="--nreactions 0"

# domain setting
OPT1="--prob_hi_x 3.2e-4 --prob_hi_y 3.2e-4 --prob_hi_z 1.e-5"
OPT1="$OPT1 --n_cells_x 32 --n_cells_y 32 --n_cells_z 1"
OPT1="$OPT1 --max_grid_size_x 16 --max_grid_size_y 16"

# equilibrium composition
OPT2="--prob_type 0"
OPT2="$OPT2 --n_init_in1_1 1.e16 --n_init_in1_2 5.e15 --n_init_in1_3 0."

# timestep
OPT3="--fixed_dt 5.e-7"

# task (plot files | structure factor & histogram)
#OPT4="--max_step 1000    --plot_int 100 --print_int 100 --hydro_grid_int 0 --n_steps_skip 0"
OPT4="--max_step 10000   --plot_int 0   --print_int 100 --hydro_grid_int 1 --n_steps_skip 1000 --histogram_unit 10"

# which algorithm (-2=unsplit explicit midpoint)
OPT6="--temporal_integrator -2"

# relevant options
OPT7="--avg_type 1"
OPT7="$OPT7 --midpoint_stoch_flux_type 2"

# misc options
OPTMISC1=
OPTMISC2=

OPTS="$OPT0 $OPT1 $OPT2 $OPT3 $OPT4 $OPT5 $OPT6 $OPT7 $OPTMISC1 $OPTMISC2"

#######
# RUN #
#######

# check executable and inputs file

if [ ! -f $EXEC ]
then
  echo "ERROR: $EXEC does not exist"
  exit
fi

if [ ! -f $INPUTS ]
then
  echo "ERROR: $INPUTS does not exist"
  exit
fi

# single run

if [ $NRUN == 1 ]
then
  if [ -d $RUNNAME ]
  then
    echo "ERROR: $RUNNAME already exists"
    exit
  fi

  mkdir $RUNNAME
  cp $INPUTS $RUNNAME
  cp $0 $RUNNAME

  cd $RUNNAME
  mpiexec -n 4 ../$EXEC ../$INPUTS $OPTS | tee scr_out

  grep n_avg scr_out > res.n_avg

  if [ -f sucrose.S_k.pair=1.Re.dat ]
  then
    ../../../../../multiSpecLM_react/test/analysis/Sk/Sk_2d.sh sucrose
  fi

  if [ -f fort.10 ]
  then
    echo "try something like (on $RUNNAME directory):"
    echo " python ../../../../../multiSpecLM_react/test/analysis/hist/hist_N.py 1 1.e15 10 -2 30 yes yes"
    echo " python ../../../../../multiSpecLM_react/test/analysis/hist/hist_N.py 2 1.e15  5 -2 20 yes yes"
    echo " python ../../../../../multiSpecLM_react/test/analysis/hist/hist_N_near_zero.py 1 1.e15 101 -1 2 yes yes"
    echo " python ../../../../../multiSpecLM_react/test/analysis/hist/hist_N_near_zero.py 2 1.e15 101 -1 2 yes yes"
  fi 

  exit
fi

# multiple run

for ((i=1;i<=$NRUN;i++))
do
  RUNDIR=${RUNNAME}_RUN$i

  if [ -d $RUNDIR ]
  then
    echo "ERROR: $RUNDIR already exists"
    exit
  fi

  echo $RUNDIR 
  mkdir $RUNDIR
  cp $INPUTS $RUNDIR
  cp $0 $RUNDIR

  cd $RUNDIR
  echo "mpiexec -n 1 ../$EXEC ../$INPUTS $OPTS > scr_out &"
  mpiexec -n 1 ../$EXEC ../$INPUTS $OPTS > scr_out &

  sleep 2s
  cd ..
done

echo
echo "try: tail -f $RUNDIR/scr_out"
