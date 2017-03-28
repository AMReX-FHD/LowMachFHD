#!/bin/bash

EXEC=../../main.Linux.gfortran.mpi.exe
INPUTS=../../inputs_SU

RUNNAME=TEST

# Dbar: should be changed in the inputs file

# rate_multiplier
OPT1="--rate_multiplier 100."

# grid 
OPT2="--prob_hi_x 6.4e-5 --prob_hi_y 6.4e-5 --n_cells_x 64 --n_cells_y 64 --max_grid_size_x 32 --max_grid_size_y 32"
#OPT2="--prob_hi_x 3.2e-5 --prob_hi_y 3.2e-5 --n_cells_x 64 --n_cells_y 64 --max_grid_size_x 32 --max_grid_size_y 32"
#OPT2="--prob_hi_x 6.4e-5 --prob_hi_y 6.4e-5 --n_cells_x 128 --n_cells_y 128 --max_grid_size_x 64 --max_grid_size_y 64"

# dz and stoch_mom_flux
DZ=1.e-5
VC=`python -c "print 1./$DZ"`
OPT3="--variance_coef_mom $VC --variance_coef_mass $VC --cross_section $DZ"
#OPT3="--variance_coef_mom 0. --variance_coef_mass $VC --cross_section $DZ"

# timestep 
OPT4="--max_step 10000 --plot_int 1000 --print_int 100 --hydro_grid_int 1 --n_steps_skip 1000 --fixed_dt 1.e-9"
#OPT4="--max_step 40000 --plot_int 4000 --print_int 400 --hydro_grid_int 4 --n_steps_skip 4000 --fixed_dt 0.25e-9"

# scheme 
OPT5="--algorithm_type 5 --include_reactions T"
#OPT5="--algorithm_type 5 --include_reactions F"
#OPT5="--algorithm_type 0 --include_reactions F"

#####

if [ ! -f $EXEC ]
then
  echo "ERROR: $EXEC does not exist"
  exit
fi

RUNNAME=RUN_$RUNNAME
if [ -d $RUNNAME ]
then
  echo "ERROR: $RUNNAME exists"
  exit
else
  mkdir $RUNNAME
  cp $0 $RUNNAME
  cp $INPUTS $RUNNAME
  cd $RUNNAME
fi

OPTS="$OPT1 $OPT2 $OPT3 $OPT4 $OPT5"

echo "mpiexec -n 4 ../$EXEC ../$INPUTS $OPTS | tee scr_out"
mpiexec -n 4 ../$EXEC ../$INPUTS $OPTS | tee scr_out

#####

grep sum scr_out | grep 'i=           1' | awk '{print $8}' > rho1
grep sum scr_out | grep 'i=           2' | awk '{print $8}' > rho2
grep sum scr_out | grep 'rho=' | awk '{print $5}' > rhotot
paste rho1 rho2 rhotot > res.rho
rm rho1 rho2 rhotot

if [ -f SU.S_k.pair=1.Re.dat ]
then
  sed -i 's/0.00000000/#0.00000000/g' SU.S_k.pair=*.Re.dat
fi
