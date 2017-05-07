#!/bin/bash

RUNNAME=TEST
NRUN=4

if [ $NRUN == 1 ]
then
  echo "ERROR: NRUN > 1 expected"
  exit
fi

if [ -d $RUNNAME ]
then
  echo "ERROR: $RUNNAME already exists"
  exit
else
  mkdir $RUNNAME
fi

########################################
# run analysis scripts for each RUNDIR #
# copy results to RUNNAME              #
########################################

PYSCR1="../../analysis/hist/hist_N.py"
PYSCR1OPTS="1 3.e-9 10. -2 30 no yes"
PYSCR2="../../analysis/hist/hist_N_near_zero.py"
PYSCR2OPTS="1 3.e-9 51 -1 2 no yes"
PYSCR3="../../analysis/hist/stat_plot.py"

for ((i=1;i<=$NRUN;i++))
do
  RUNDIR=${RUNNAME}_RUN$i

  if [ ! -d $RUNDIR ]
  then
    echo "ERROR: $RUNDIR does not exist"
    exit
  fi

  echo $RUNDIR
  cd $RUNDIR

  # copy scr_out
  cp scr_out ../$RUNNAME/scr_out$i

  # run script hist_N.py and copy results
  echo "=== running $PYSCR1 ==="
  echo "!! check option values: $PYSCR1OPTS"
  python ../$PYSCR1 $PYSCR1OPTS 
  cp res.hist ../$RUNNAME/res.hist$i
  cp res.N_stat ../$RUNNAME/res.N_stat$i
  if [ $i -eq 1 ]
  then
    cp res.hist_cont ../$RUNNAME/res.hist_cont
  fi

  # run script hist_N_near_zero.py and copy results
  echo "=== running $PYSCR2 ==="
  echo "!! check option values: $PYSCR2OPTS"
  python ../$PYSCR2 $PYSCR2OPTS
  cp res.hist_near_zero ../$RUNNAME/res.hist_near_zero$i

  echo 
  cd ..
done

#########################################################
# statistics (sample mean and variance)                 #
# 1st col = quantity_name, k for S(k), or n for hist    #
# 2nd col = sample mean                                 #
# 3rd col = std error (normalized by NRUN               #
# additional NRUN columns from each sample run          #
#########################################################

cd $RUNNAME

## hist

cmd_awk=""
for ((i=1;i<=$NRUN;i++))
do
  cmd_tmp="<(awk 'NR>=2{print \$2}' res.hist$i)"
  cmd_awk="$cmd_awk $cmd_tmp"
done

TMP="tmp.hist_stat"
eval paste $cmd_awk | awk '{sum1=0;sum2=0;for(i=1;i<=NF;i++){sum1+=$i;sum2+=$i*$i}sum1/=NF;sum2/=NF;printf "%e\t%e\t",sum1,sqrt((sum2-sum1*sum1)/NF);print $0}' > $TMP
paste <(awk 'NR>=2{print $1}' res.hist1) <(awk '{print $0}' $TMP) > res.hist_stat
rm $TMP

## N_stat

cmd_awk=""
for ((i=1;i<=$NRUN;i++))
do
  cmd_tmp="<(awk '{print \$2}' res.N_stat$i)"
  cmd_awk="$cmd_awk $cmd_tmp"
done

TMP="tmp.n_stat_stat"
eval paste $cmd_awk | awk '{sum1=0;sum2=0;for(i=1;i<=NF;i++){sum1+=$i;sum2+=$i*$i}sum1/=NF;sum2/=NF;printf "%e\t%e\t",sum1,sqrt((sum2-sum1*sum1)/NF);print $0}' > $TMP
paste <(awk '{print $1}' res.N_stat1) <(awk '{print $0}' $TMP) > res.N_stat_stat
rm $TMP

## hist_near_zero

cmd_awk=""
for ((i=1;i<=$NRUN;i++))
do
  cmd_tmp="<(awk '{print \$2}' res.hist_near_zero$i)"
  cmd_awk="$cmd_awk $cmd_tmp"
done

TMP="tmp.hist_near_zero_stat"
eval paste $cmd_awk | awk '{sum1=0;sum2=0;for(i=1;i<=NF;i++){sum1+=$i;sum2+=$i*$i}sum1/=NF;sum2/=NF;printf "%e\t%e\t",sum1,sqrt((sum2-sum1*sum1)/NF);print $0}' > $TMP
paste <(awk '{print $1}' res.hist_near_zero1) <(awk '{print $0}' $TMP) > res.hist_near_zero_stat
rm $TMP

########
# plot #
########

echo "=== running $PYSCR3 ==="
python ../$PYSCR3
