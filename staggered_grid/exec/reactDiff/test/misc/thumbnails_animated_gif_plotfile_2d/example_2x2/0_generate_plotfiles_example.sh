###########################################
# obtain plotfiles from a total of 4 runs #
###########################################

### input part ###

PROG=../../../main.Linux.gfortran.debug.mpi.exe # reactdiff program
INPUTS=../../../inputs_BPM_Turing               # inputs file

declare -a RUN_DIR  # directories for plotfiles #  +----+----+
RUN_DIR[1]=run1                                 #  |run3|run4|
RUN_DIR[2]=run2                                 #  +----+----+
RUN_DIR[3]=run3                                 #  |run1|run2|
RUN_DIR[4]=run4                                 #  +----+----+

declare -a RUN_OPT  # options for each run 
RUN_OPT[1]="--variance_coef_mass 1.d0 --use_Poisson_rng -1"
RUN_OPT[2]="--variance_coef_mass 0.d0 --use_Poisson_rng -1"
RUN_OPT[3]="--variance_coef_mass 1.d0 --use_Poisson_rng 1"
RUN_OPT[4]="--variance_coef_mass 0.d0 --use_Poisson_rng 1"

                    # common options for all runs
#COMMON_OPT="--max_step 100 --plot_int 10  --print_int 10"


### exec part ###

for i in {1..4}
do
  # go into run directory
  WORKDIR=`pwd`/${RUN_DIR[$i]}
  if [ -d $WORKDIR ]; then rm -rf $WORKDIR; fi
  mkdir $WORKDIR 
  cd $WORKDIR

  # execute inputs file and obtain plotfiles
  echo "** RUN$i: mpiexec -n 4 ../$PROG ../$INPUTS $COMMON_OPT ${RUN_OPT[$i]}" 
  mpiexec -n 4 ../$PROG ../$INPUTS $COMMON_OPT ${RUN_OPT[$i]} > res_scr

  # get back to main directory and create visit file
  cd ..
  ls -1 $WORKDIR/plt*/Header > run${i}.visit 
done
