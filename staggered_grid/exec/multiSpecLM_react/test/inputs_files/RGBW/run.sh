#!/bin/bash

EXEC=../../main.Linux.gfortran.mpi.exe
INPUTS=../inputs_RGBW
RUNNAME=TEST
NRUN=4

###########
# OPTIONS #
###########

# which algorithm (5=inertial midpoint, 2=overdamped, 0=inertial trapezoidal)
OPT1="--algorithm_type 5"
#OPT1="--algorithm_type 2"
#OPT1="--algorithm_type 0"

# reaction on/off
OPT2="--nreactions 2"
#OPT2="--nreactions 0"

# task (check gmres solver | plot files | structure factor & histogram)
#OPT3="--max_step 10      --plot int 0   --print_int 10  --hydro_grid_int 0 --n_steps_skip 0   --gmres_verbose 1"
#OPT3="--max_step 1000    --plot_int 100 --print_int 100 --hydro_grid_int 0 --n_steps_skip 0"
OPT3="--max_step 1000    --plot_int 0   --print_int 100 --hydro_grid_int 1 --n_steps_skip 100 --histogram_unit 10"

# dt, dz, rate_multiplier, avg_type
OPT4="--fixed_dt 1.e-8 --prob_hi_z 0.01 --rate_multiplier 1. --avg_type 1"

# turn off gmres solver (velocity becomes exactly zero)
#OPT5="--variance_coef_mom 0. --initial_variance 0. --gmres_abs_tol 1.e-10 --mg_abs_tol 1.e-10"

# misc options
OPTMISC=""

OPTS="$OPT1 $OPT2 $OPT3 $OPT4 $OPT5 $OPTMISC"

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

  grep sum scr_out | grep 'i=           1' | awk '{print $8}' > rho1
  grep sum scr_out | grep 'i=           2' | awk '{print $8}' > rho2
  grep sum scr_out | grep 'i=           3' | awk '{print $8}' > rho3
  grep sum scr_out | grep 'i=           4' | awk '{print $8}' > rho4
  grep sum scr_out | grep 'rho=' | awk '{print $5}' > rhotot
  paste rho1 rho2 rho3 rho4 rhotot > res.rho
  rm rho1 rho2 rho3 rho4 rhotot

  if [ -f RGBW.S_k.pair=1.Re.dat ]
  then
    ../../../analysis/Sk/Sk_2d.sh RGBW 
  fi

  if [ -f fort.10 ]
  then
    python ../../../analysis/hist/hist_N.py 1 3.e-9 10. -2 30 yes yes
    python ../../../analysis/hist/hist_N_near_zero.py 1 3.e-9 101 -1 2 yes yes
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
