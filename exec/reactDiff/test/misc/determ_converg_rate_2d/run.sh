#!/bin/bash

###############################################################################

PROG=../../main.Linux.gfortran.debug.mpi.exe    
                                       # reactdiff program
INPUTS=../../inputs_prob_type_Gaussian_2d       
                                       # inputs file for Gaussian prob_type 
EXEC_DIFF=DiffSameGridRefined2d
                                       # executable to calculate the difference between two plot files
                                       # (BoxLib/Tools/C_util/Convergence$/DiffSameGridRefined.cpp)

OPTIONS_DCRT="--initial_variance 0.d0 --variance_coef_mass 0.d0 --use_Poisson_rng -1 --include_discrete_LMA_correction F"
                                       # OPTIONS_DCRT: options needed for deterministic convergence rate test 
                                       # so don't comment this out and don't change this!
                                       # initial_variance=0: no fluctuation for n_init
                                       # variance_coef_mass=0: no stochastic flux
                                       # use_Poisson_rng=-1: deterministic chemistry
                                       # include_discrete_LMA_correction=F: no dependence of rate expressions on dv

                                       # OPTIONS: put some other options you like to change from the input file
                                       # 2nd-order reaction-only case
#MY_OPTIONS="--no_diffusion --temporal_integrator 0 --reaction_type 1"
                                       # 2nd-order diffusion-only case              
#MY_OPTIONS="--nreactions 0 --temporal_integrator 0 --diffusion_type 1"              
                                       # 2nd-order reaction-diffusion case (splitting)
MY_OPTIONS="--temporal_integrator 1 --reaction_type 1 --diffusion_type 1"
                                       # 2nd-order reaction-diffusion case (unsplitting)              
#MY_OPTIONS="--temporal_integrator -2"  # since this is an explicit scheme, decrease dt_fixed. see below. 

NTEST=4                                # total number of runs 

###############################################################################

### checking files
if [ ! -f "$PROG" ];   then echo "ERROR: reactdiff program not found"; exit; fi
if [ ! -f "$INPUTS" ]; then echo "ERROR: inputs file not found";       exit; fi
tmp=`ls -d plt0* 2>/dev/null |wc -l`
if [ "$tmp" != 0 ]; then echo "ERROR: remove plot files before running this script"; exit; fi

### read n_cells from a line indicated by N_CELLS 
# count lines containing N_CELLS 
tmp=`grep N_CELLS $INPUTS | wc | awk '{print $1}'`
# if there are more than one line, abort
if [ $tmp != "1" ]; then echo "ERROR: N_CELLS should appear only in a single line in the script."; exit; fi
# get n_cells (get line|omit up to =|get the first two values)
tmp=`grep N_CELLS $INPUTS | awk '{$0=substr($0,index($0,"=")+1);print $0}' | awk '{print $1,$2}'`
NX=`echo $tmp | awk '{print $1}'`  
NY=`echo $tmp | awk '{print $2}'`

### read max_grid_size from a line indicated by MAX_GRID_SIZE 
# count lines containing MAX_GRID_SIZE 
tmp=`grep MAX_GRID_SIZE $INPUTS | wc | awk '{print $1}'`
# if there are more than one line, abort
if [ $tmp != "1" ]; then echo "ERROR: MAX_GRID_SIZE should appear only in a single line in the script."; exit; fi
# get max_grid_size (get line|omit up to =|get the first two values)
tmp=`grep MAX_GRID_SIZE $INPUTS | awk '{$0=substr($0,index($0,"=")+1);print $0}' | awk '{print $1,$2}'`
MGSX=`echo $tmp | awk '{print $1}'`  
MGSY=`echo $tmp | awk '{print $2}'`

### read dt_fixed from a line indicated by FIXED_DT
# count lines containing FIXED_DT
tmp=`grep FIXED_DT $INPUTS | wc | awk '{print $1}'`
# if there are more than one line, abort
if [ $tmp != "1" ]; then echo "ERROR: FIXED_DT should appear only in a single line in the script."; exit; fi
# get fixed_dt (get line|omit up to =|get 1st value|change fortran format (1.d0->1.e0)|(1.D0->1.e0)|decimal)
DT=`grep FIXED_DT $INPUTS | awk '{$0=substr($0,index($0,"=")+1);print $0}' | awk '{print $1}' | sed s/d/e/g | sed s/D/e/g | awk '{printf("%f",$1)}'`

### read max_step from a line indicated by MAX_STEP 
# count lines containing MAX_STEP 
tmp=`grep MAX_STEP $INPUTS | wc | awk '{print $1}'`
# if there are more than one line, abort
if [ $tmp != "1" ]; then echo "ERROR: MAX_STEP should appear only in a single line in the script."; exit; fi
# get max_step (get line|omit up to =|get 1st value)
MAXSTEP=`grep MAX_STEP $INPUTS | awk '{$0=substr($0,index($0,"=")+1);print $0}' | awk '{print $1}'`

### if you want to assign new values, put them here
# NX=16
# NY=16
# MGSX=8
# MGSY=8
# DT=0.0001      # for explicit schemes, use smaller fixed_dt and increase max_step accordingly.
# MAXSTEP=5000

echo "## n_cells=($NX,$NY), max_grid_size=($MGSX,$MGSY), dt_fixed=$DT, max_step=$MAXSTEP"
echo "## OPTIONS_DCRT=$OPTIONS_DCRT"
echo "## MY_OPTIONS=$MY_OPTIONS"
###############################################################################

# a series of runs
# screen output for each run is saved in res_scr.run$i

INPUTS_TMP=inputs_tmp

for ((i=1;i<=$NTEST;i++))
do
  echo "** run$i: n_cells=($NX,$NY), max_grid_size=($MGSX,$MGSY), dt_fixed=$DT, max_step=$MAXSTEP" 

  # replace values for c_cells, max_grid_size
  awk -v NX=$NX -v NY=$NY '{if(match($0,"N_CELLS")){print "n_cells(1:2) = " NX " " NY}else{print $0}}' $INPUTS | \
  awk -v MGSX=$MGSX -v MGSY=$MGSY '{if(match($0,"MAX_GRID_SIZE")){print "max_grid_size(1:2) = " MGSX " " MGSY}else{print $0}}' > $INPUTS_TMP

  # options for new values of fixed_dt, max_step, plot_int
  OPTIONS_ADD="--fixed_dt $DT --max_step $MAXSTEP --plot_int $MAXSTEP"

  mpiexec -n 4 $PROG $INPUTS_TMP $OPTIONS_DCRT $OPTIONS_ADD $MY_OPTIONS > res_scr.run$i

  # update variable values 
  NX=$((NX*2)); NY=$((NY*2))
  MGSX=$((MGSX*2)); MGSY=$((MGSY*2))
  DT=$(bc -l <<< "$DT/2")
  MAXSTEP=$((MAXSTEP*2))
done

rm $INPUTS_TMP

###############################################################################

# compute the difference between final results of successive runs

RES_DIFF=res_diff
echo "** difference **" | tee $RES_DIFF

list=`ls -d plt0*`          # list contains (NTEST+1) plot files.
                            # the first item is plt00000000 and will be skipped. 
for ((i=1;i<$NTEST;i++))
do
  i1=$((i+1))
  i2=$((i+2))
  file1=`echo $list | awk -v i1=$i1 '{print $i1}'`  # (i+1)th plot file in list
  file2=`echo $list | awk -v i2=$i2 '{print $i2}'`  # (i+2)th plot file in list

  # calculate the difference between two plot files 
  res=`$EXEC_DIFF infile1= $file1 reffile= $file2 | tail -2 | head -1 | awk '{for(j=2;j<=NF;j++){printf("%.15f\t",$j)}}'`
  echo $res | tee -a $RES_DIFF
done

# compute the ratio of differences between successive runs

echo "** ratio **" | tee -a $RES_DIFF
awk -v NTEST=$NTEST 'NR==2{for(j=1;j<=NF;j++)prev[j]=$j} NR>2 && NR <= NTEST {for(j=1;j<=NF;j++){printf("%g\t",prev[j]/$j);prev[j]=$j}printf("\n")}' $RES_DIFF | tee tmp
cat tmp >> $RES_DIFF
rm tmp
