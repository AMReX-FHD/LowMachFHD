#!/bin/bash

EXEC=../../main.Linux.gfortran.mpi.exe
INPUTS=../inputs_SU
RUNNAME=TEST

###########
# OPTIONS #
###########

# cross section
OPT1="--prob_hi_z 1e-5"

# timestep size
OPT2="--fixed_dt 1.e-8"

# task (check gmres solver | plot files | structure factor )
#OPT3="--max_step 10      --plot int 0   --print_int 10  --hydro_grid_int 0 --n_steps_skip 0    --gmres_verbose 1"
OPT3="--max_step 1000    --plot_int 100 --print_int 100 --hydro_grid_int 0 --n_steps_skip 0"
#OPT3="--max_step 10000   --plot_int 0   --print_int 100 --hydro_grid_int 1 --n_steps_skip 1000"

# turn off gmres solver (velocity becomes exactly zero)
#OPT4="--variance_coef_mom 0. --initial_variance 0. --gmres_abs_tol 1.e-10 --mg_abs_tol 1.e-10"

OPTMISC1=
OPTMISC2=

OPTS="$OPT1 $OPT2 $OPT3 $OPT4 $OPTMISC1 $OPTMISC2"

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
mpiexec -n 4 ../$EXEC ../$INPUTS $OPTS | tee scr_out

grep sum scr_out | grep 'i=           1' | awk 'NR>2{print $8}' > rho1
grep sum scr_out | grep 'i=           2' | awk 'NR>2{print $8}' > rho2
grep sum scr_out | grep 'rho=' | awk 'NR>2{print $5}' > rhotot
paste rho1 rho2 rhotot > res.rho
rm rho1 rho2 rhotot

if [ -f SU.S_k.pair=1.Re.dat ]
then
  ../../../analysis/Sk/Sk_2d.sh SU 
fi
