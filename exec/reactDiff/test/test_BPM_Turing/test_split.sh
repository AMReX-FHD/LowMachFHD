#!/bin/bash

###############################################################################

PROG=../main.Linux.gfortran.debug.mpi.exe
INPUTS=../inputs_BPM_Turing

# choose "short_run" for short runs and "run" for actual runs 
JOB=short_run
#JOB=run

# prefix for test run directories
DIR0="RUN"

###############################################################################

# command to be executed test/test_BPM_Turing/RUN_***
command="mpiexec -n 4 ../$PROG ../$INPUTS"

# default options
case "$JOB" in 
  "short_run" )
    option0="--seed 1 --fixed_dt 0.5 --plot_int 200 --print_int 20 --max_step 200" ;;
  "run" )
    option0="--seed 1 --fixed_dt 0.5 --plot_int 200 --print_int 20 --max_step 10000" ;;
esac

###############################################################################

count=0

for tempintg in split1 split2 split3
do
  case "$tempintg" in
    "split1" )
      option1="--temporal_integrator 0" ;;  # first-order splitting
    "split2" )
      option1="--temporal_integrator 1" ;;  # Strang splitting 1
    "split3" )
      option1="--temporal_integrator 2" ;;  # Strang splitting 2
  esac

  for diffmid in diff0 diff1 diff2_mid1 diff2_mid2 diff2_mid3
  do
    case "$diffmid" in
      "diff0" )
        option2="--diffusion_type 0" ;;  # explicit trapezoidal predictor/corrector
      "diff1" )
        option2="--diffusion_type 1" ;;  # Crank-Nicolson semi-implicit
      "diff2_mid1" )
        option2="--diffusion_type 2 --midpoint_stoch_flux_type 1" ;;  # explicit midpoint | nold
      "diff2_mid2" )
        option2="--diffusion_type 2 --midpoint_stoch_flux_type 2" ;;  # explicit midpoint | npred
      "diff2_mid3" )
        option2="--diffusion_type 2 --midpoint_stoch_flux_type 3" ;;  # explicit midpoint | 2*npred-nold
    esac

    for reactrng in react0_rngP react0_rngG react0_rngD react1_rngP react1_rngG react1_rngD react2
    do
      case "$reactrng" in
        "react0_rngP" )
          option3="--reaction_type 0 --use_Poisson_rng 1" ;;   # first-order tau-leaping
        "react0_rngG" )
          option3="--reaction_type 0 --use_Poisson_rng 0" ;;   # first-order CLE
        "react0_rngD" )
          option3="--reaction_type 0 --use_Poisson_rng -1" ;;  # first-order deterministic 
        "react1_rngP" )
          option3="--reaction_type 1 --use_Poisson_rng 1" ;;   # second-order tau-leaping
        "react1_rngG" )
          option3="--reaction_type 1 --use_Poisson_rng 0" ;;   # second-order CLE 
        "react1_rngD" )
          option3="--reaction_type 1 --use_Poisson_rng -1" ;;  # second-order deterministic 
        "react2") 
          option3="--reaction_type 2" ;;  # SSA 
      esac 

      ### RUN ###

      dir=${DIR0}_${tempintg}_${diffmid}_${reactrng}
      options="$option0 $option1 $option2 $option3" 

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

done 
