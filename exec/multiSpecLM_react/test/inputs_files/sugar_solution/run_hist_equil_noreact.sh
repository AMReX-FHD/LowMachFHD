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
OPT2="--prob_type 2"
OPT2="$OPT2 --c_init1_1 2.991559e-06 --c_init1_2 2.841983e-06 --c_init1_3 0. --c_init1_4 9.9999416645757710e-01"
OPT2="$OPT2 --c_init2_1 2.991559e-06 --c_init2_2 2.841983e-06 --c_init2_3 0. --c_init2_4 9.9999416645757710e-01"

# timestep
OPT3="--fixed_dt 5.e-7"

# task (check gmres solver | plot files | structure factor & histogram)
#OPT4="--max_step 10      --plot int 0   --print_int 10  --hydro_grid_int 0 --n_steps_skip 0    --gmres_verbose 1"
#OPT4="--max_step 1000    --plot_int 100 --print_int 100 --hydro_grid_int 0 --n_steps_skip 0"
OPT4="--max_step 10000   --plot_int 0   --print_int 100 --hydro_grid_int 1 --n_steps_skip 1000 --histogram_unit 10"

# turn off gmres solver (velocity becomes exactly zero)
OPT5="--variance_coef_mom 0. --initial_variance 0. --gmres_abs_tol 1.e-10 --mg_abs_tol 1.e-10"

# which algorithm (5=inertial midpoint, 2=overdamped, 0=inertial trapezoidal)
OPT6="--algorithm_type 5"
#OPT6="--algorithm_type 2"
#OPT6="--algorithm_type 0"

# relevant options
OPT7="--avg_type 1"
OPT7="$OPT7 --midpoint_stoch_mass_flux_type 2"

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

  grep sum scr_out | grep 'i=           1' | awk '{print $8}' > rho1
  grep sum scr_out | grep 'i=           2' | awk '{print $8}' > rho2
  grep sum scr_out | grep 'i=           3' | awk '{print $8}' > rho3
  grep sum scr_out | grep 'i=           4' | awk '{print $8}' > rho4
  grep sum scr_out | grep 'rho=' | awk '{print $5}' > rhotot
  paste rho1 rho2 rho3 rho4 rhotot > res.rho
  rm rho1 rho2 rho3 rho4 rhotot

  if [ -f sucrose.S_k.pair=1.Re.dat ]
  then
    ../../../analysis/Sk/Sk_2d.sh sucrose
  fi

  if [ -f fort.10 ]
  then
    echo "try something like (on $RUNNAME directory):"
    echo " python ../../../analysis/hist/hist_N.py 1 2.991559e-07 10 -2 30 yes yes"
    echo " python ../../../analysis/hist/hist_N.py 2 5.683967e-07  5 -2 20 yes yes"
    echo " python ../../../analysis/hist/hist_N_near_zero.py 1 2.991559e-07 101 -1 2 yes yes"
    echo " python ../../../analysis/hist/hist_N_near_zero.py 2 5.683967e-07 101 -1 2 yes yes"
    echo " python ../../../analysis/hist/hist_N_near_zero.py 3 2.991559e-07 101 -1e-11 1e-11 yes yes"
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
