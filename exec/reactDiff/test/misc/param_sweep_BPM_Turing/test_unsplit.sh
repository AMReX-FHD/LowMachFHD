#!/bin/bash

###############################################################################

PROG=../../main.Linux.gfortran.debug.mpi.exe
INPUTS=../../inputs_BPM_Turing

# choose "short_run" for short runs and "run" for actual runs 
JOB=short_run
#JOB=run

# prefix for test run directories
DIR0="RUN"

###############################################################################

# command to be executed on ${DIR0}_*** 
command="mpiexec -n 4 ../$PROG ../$INPUTS"

# basic options to be applied to all runs
# you may add more options (e.g., variance_coef_mass = 0.d0)
case "$JOB" in 
  "short_run" )
    option0="--seed 1 --fixed_dt 0.5d0 --plot_int 200 --print_int 20 --max_step 200" ;;
  "run" )
    option0="--seed 1 --fixed_dt 0.5d0 --plot_int 200 --print_int 20 --max_step 10000" ;;
esac

###############################################################################

# checking files
if [ ! -f "$PROG" ];   then echo "ERROR: reactdiff program not found"; exit; fi
if [ ! -f "$INPUTS" ]; then echo "ERROR: inputs file not found";       exit; fi

# main loops for paramter sweeping

count=0

for tempintgmid in unsplit1 unsplit2_mid1 unsplit2_mid2 unsplit2_mid3
do
  case "$tempintgmid" in
    "unsplit1" )
      option1="--temporal_integrator -1" ;;  # unsplitting forward Euler
    "unsplit2_mid1" )
      option1="--temporal_integrator -2 --midpoint_stoch_flux_type 1" ;;  # unsplitting explicit midpoint | nold
    "unsplit2_mid2" )
      option1="--temporal_integrator -2 --midpoint_stoch_flux_type 2" ;;  # unsplitting explicit midpoint | npred
    "unsplit2_mid3" )
      option1="--temporal_integrator -2 --midpoint_stoch_flux_type 3" ;;  # unsplitting explicit midpoint | 2*npred-nold
  esac

  for rng in rngP rngG rngD 
  do
    case "$rng" in
      "rngP" )
        option2="--use_Poisson_rng 1" ;;   # first-order tau-leaping
      "rngG" )
        option2="--use_Poisson_rng 0" ;;   # first-order CLE
      "rngD" )
        option2="--use_Poisson_rng -1" ;;  # first-order deterministic 
    esac 

    ### RUN ###

    dir=${DIR0}_${tempintgmid}_${rng}
    options="$option0 $option1 $option2" 

    count=$((count+1))
    echo " "
    echo "case $count: $dir"
    mkdir $dir
    cd $dir
    echo $command $options | tee res.scr
    time $command $options >> res.scr
    cd ..

  done  

done 
