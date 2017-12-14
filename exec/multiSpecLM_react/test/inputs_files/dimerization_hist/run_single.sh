#!/bin/bash

EXEC=../../main.Linux.gfortran.mpi.exe
INPUTS=../inputs_dimerization_hist
RUNNAME=TEST
NPROC=4

###########
# OPTIONS #
###########

# number of atoms per cell
#OPT1="--rho0 20. --rate_const_1 1. --rate_const_2 1.290564"    # N0=20
OPT1="--rho0 40. --rate_const_1 1. --rate_const_2 1.311528"     # N0=40
#OPT1="--rho0 60. --rate_const_1 1. --rate_const_2 1.318704"    # N0=60

# reaction on/off
#OPT2="--rate_multiplier 0."            # reaction off
#OPT2="--rate_multiplier 0.421875"      # N0=20, d = sqrt(10)*dx, k1 = 27*N0/1280
OPT2="--rate_multiplier 0.84375"        # N0=40
#OPT2="--rate_multiplier 1.265625"      # N0=60

# advection type
OPT3="--advection_type 0"

# temperature
OPT4="--T_init_1 1.e4 --T_init_2 1.e4"

# time step size and task (check gmres solver | plot files | histogram)
#OPT5="--fixed_dt 1.e-2 --max_step 10     --plot int 0   --print_int 10   --hydro_grid_int 0   --gmres_verbose 1"
#OPT5="--fixed_dt 1.e-2 --max_step 10000  --plot_int 100 --print_int 100  --hydro_grid_int 0"
OPT5="--fixed_dt 1.e-2 --max_step 100000  --plot_int 0   --print_int 1000 --hydro_grid_int 100 --n_steps_skip 20000 --histogram_unit 10"

# Stratonovich or Ito
OPT6="--midpoint_stoch_mass_flux_type 1"
#OPT6="--midpoint_stoch_mass_flux_type 2"

OPTMISC1="--max_grid_size_x 16 --max_grid_size_y 16"
#OPTMISC2="--use_Poisson_rng 0"

OPTS="$OPT1 $OPT2 $OPT3 $OPT4 $OPT5 $OPT6 $OPTMISC1 $OPTMISC2"

#######
# RUN #
#######

# check executable / inputs file / run directory

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

if [ -d $RUNNAME ]
then
  echo "ERROR: $RUNNAME already exists"
  exit
fi

mkdir $RUNNAME
cp $INPUTS $RUNNAME
cp $0 $RUNNAME

cd $RUNNAME
mpiexec -n $NPROC ../$EXEC ../$INPUTS $OPTS | tee scr_out

grep sum scr_out | grep 'i=           1' | awk 'NR>2{print $8}' > rho1
grep sum scr_out | grep 'i=           2' | awk 'NR>2{print $8}' > rho2
grep sum scr_out | grep 'rho=' | awk 'NR>2{print $5}' > rhotot
paste rho1 rho2 rhotot > res.rho
rm rho1 rho2 rhotot

if [ -f fort.10 ]
then
  echo "try: (N0=20) cd $RUNNAME; python ../../../analysis/hist/hist_N.py 2 2. 5 0 10 yes yes"
  echo "     gnuplot> set log y; plot \"../res.hist_CME_20\" u ((20-\$1)/2):2 w l,\"res.hist\" w l"
  echo "     (N0=40) cd $RUNNAME; python ../../../analysis/hist/hist_N.py 2 2. 10 0 20 yes yes"
  echo "     gnuplot> set log y; plot \"../res.hist_CME_40\" u ((40-\$1)/2):2 w l,\"res.hist\" w l"
  echo "     (N0=60) cd $RUNNAME; python ../../../analysis/hist/hist_N.py 2 2. 15 0 30 yes yes"
  echo "     gnuplot> set log y; plot \"../res.hist_CME_60\" u ((60-\$1)/2):2 w l,\"res.hist\" w l"
fi
