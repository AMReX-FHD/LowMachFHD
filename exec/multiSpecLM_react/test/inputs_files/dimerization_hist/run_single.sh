#!/bin/bash

EXEC=../../main.Linux.gfortran.mpi.exe
INPUTS=../inputs_dimerization_hist
RUNNAME=TEST
NPROC=4

###########
# OPTIONS #
###########

# cross section
OPT1="--prob_hi_z 40."

# reaction
#OPT2="--rate_multiplier 0."
OPT2="--rate_multiplier 0.1"

# timestep size
OPT3="--fixed_dt 1.e-2"

# task (check gmres solver | plot files | structure factor )
#OPT4="--max_step 10      --plot int 0   --print_int 10   --hydro_grid_int 0   --gmres_verbose 1"
#OPT4="--max_step 1000    --plot_int 100 --print_int 10   --hydro_grid_int 0"
OPT4="--max_step 100000  --plot_int 0   --print_int 1000 --hydro_grid_int 100 --n_steps_skip 20000 --histogram_unit 10"

OPTMISC1="--max_grid_size_x 16 --max_grid_size_y 16"
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
mpiexec -n $NPROC ../$EXEC ../$INPUTS $OPTS | tee scr_out

grep sum scr_out | grep 'i=           1' | awk 'NR>2{print $8}' > rho1
grep sum scr_out | grep 'i=           2' | awk 'NR>2{print $8}' > rho2
grep sum scr_out | grep 'rho=' | awk 'NR>2{print $5}' > rhotot
paste rho1 rho2 rhotot > res.rho
rm rho1 rho2 rhotot

if [ -f fort.10 ]
then
  echo "try: cd $RUNNAME; python ../../../analysis/hist/hist_N.py 2 0.05 10 0 20 yes yes"
  echo "     gnuplot> plot \"../res.hist_CME_40\" u ((40-\$1)/2):2 w l,\"res.hist\" w l"
fi
