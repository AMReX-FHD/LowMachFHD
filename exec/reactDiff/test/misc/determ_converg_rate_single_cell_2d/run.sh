#!/bin/bash

###############################################################################

PROG=../../main.Linux.gfortran.debug.mpi.exe    # reactdiff program
INPUTS=../../inputs_single_cell_2d              # inputs file for single cell case

RES=res                                   # output file which will contain number densities 
                                          # at the file name T=DT0*MAXSTEP0 obtained from 
                                          # each of NTEST runs with decreasing dt

OPTIONS0="--use_Poisson_rng -1 --print_int 1"
                                          # OPTIONS0 are needed for deterministic convergence rate test 
                                          # so don't comment this out and don't change this!

#OPTIONS="--reaction_type 0"              
OPTIONS="--reaction_type 1"  
                                          # OPTIONS: put any options you like to change from the input file

NTEST=4             # total number of runs 
DT0=0.01            # largest timestep size
                    # dt will be decreased: DT0, DT0/2, ... , DT0/2^(NTEST-1)
MAXSTEP0=100        # number of timesteps for the largest timestep size case
                    # max_step will be increased: MAXSTEP0, 2*MAXSTEP0, ... , MAXSTEP0*2^(NTEST-1)

###############################################################################

# checking files
if [ ! -f "$PROG" ];   then echo "ERROR: reactdiff program not found"; exit; fi
if [ ! -f "$INPUTS" ]; then echo "ERROR: inputs file not found";       exit; fi 
if [ -f "$RES" ];      then rm $RES; fi

# assigning initial values for DT and MAXSTEP
DT=$(bc -l <<< "$DT0")
MAXSTEP=$MAXSTEP0

# execute the program several times with changing DT and MAXSTEP
for ((i=1;i<=$NTEST;i++))
do
  MAXSTEP1=$((MAXSTEP+1))

  echo "##### dt=$DT max_step=$MAXSTEP" >> $RES

  echo "$PROG $INPUTS $OPTIONS0 $OPTIONS --fixed_dt $DT --max_step $MAXSTEP1 | grep n_avg | tail -1 >> $RES"
  $PROG $INPUTS $OPTIONS0 $OPTIONS --fixed_dt $DT --max_step $MAXSTEP1 | grep n_avg | tail -1 >> $RES

  DT=$(bc -l <<< "$DT/2")
  MAXSTEP=$((MAXSTEP*2))
done

# calculate the differences from the previous runs
echo " "
echo "** differences decrease as follows (1..nspecies):"
grep n_avg $RES \
  | awk '{print substr($0,index($0,$3))}' \
  | awk 'NR==1{for(i=1;i<=NF;i++)prev[i]=$i}NR>1{for(i=1;i<=NF;i++){printf("%g\t",$i-prev[i]);prev[i]=$i}printf("\n")}'

# cacluate the ratios for the differences
echo " "
echo "** ratios are as follows (1..nspecies):"
grep n_avg $RES \
  | awk '{print substr($0,index($0,$3))}' \
  | awk 'NR==1{for(i=1;i<=NF;i++)prev[i]=$i}NR>1{for(i=1;i<=NF;i++){printf("%g\t",$i-prev[i]);prev[i]=$i}printf("\n")}' \
  | awk 'NR==1{for(i=1;i<=NF;i++)prev[i]=$i}NR>1{for(i=1;i<=NF;i++){printf("%g\t",prev[i]/$i);prev[i]=$i}printf("\n")}'
