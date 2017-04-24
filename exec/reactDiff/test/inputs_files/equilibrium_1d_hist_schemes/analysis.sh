#!/bin/bash

RUNNAME=ExMidTau_dt0.1
NRUN=16

NAV=1.
DV=10.
NCELL=64

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

PYSCR1="../../analysis/hist_1d/hist_n.py"
PYSCR2="../../analysis/hist_1d/hist_n_near_zero.py"
PYSCR3="../../analysis/hist_1d/corr.py"
#PYSCR4="../../analysis/hist_1d/face_avg.py"
PYSCR5="../../analysis/hist_1d/stat_plot.py"

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

  # copy S(k) 
  cp Hist1D.S_k.pair\=1.Re.dat ../$RUNNAME/res.Sk$i

  # run script hist_n.py and copy results
  echo "=== running $PYSCR1 ==="
  python ../$PYSCR1 $NAV $DV no yes
  cp res.hist ../$RUNNAME/res.hist$i
  cp res.n_stat ../$RUNNAME/res.n_stat$i
  if [ $i -eq 1 ]
  then
    cp res.hist_poiss ../$RUNNAME/res.hist_poiss
    cp res.hist_cont ../$RUNNAME/res.hist_cont
  fi

  # run script hist_n_near_zero.py and copy results
  echo "=== running $PYSCR2 ==="
  python ../$PYSCR2
  cp res.hist_near_zero ../$RUNNAME/res.hist_near_zero$i

  # run script corr.py and copy results
  echo "=== running $PYSCR3 ==="
  python ../$PYSCR3 $NCELL > res.corr
  cp res.corr ../$RUNNAME/res.corr$i

  # run script face_avg.py and copy results
  #echo "=== running $PYSCR4 ==="
  #python ../$PYSCR4 $NCELL $DV $AVG > res.face_avg
  #cp res.face_avg ../$RUNNAME/res.face_avg$i

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

## S(k)

cmd_awk=""
for ((i=1;i<=$NRUN;i++))
do
  cmd_tmp="<(awk 'NR>=3{print \$2}' res.Sk$i)"
  cmd_awk="$cmd_awk $cmd_tmp"
done

TMP="tmp.Sk_stat"
eval paste $cmd_awk | awk '{sum1=0;sum2=0;for(i=1;i<=NF;i++){sum1+=$i;sum2+=$i*$i}sum1/=NF;sum2/=NF;printf "%e\t%e\t",sum1,sqrt((sum2-sum1*sum1)/NF);print $0}' > $TMP
paste <(awk 'NR>=3{print $1}' res.Sk1) <(awk '{print $0}' $TMP) > res.Sk_stat
rm $TMP

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

## n_stat

cmd_awk=""
for ((i=1;i<=$NRUN;i++))
do
  cmd_tmp="<(awk '{print \$2}' res.n_stat$i)"
  cmd_awk="$cmd_awk $cmd_tmp"
done

TMP="tmp.n_stat_stat"
eval paste $cmd_awk | awk '{sum1=0;sum2=0;for(i=1;i<=NF;i++){sum1+=$i;sum2+=$i*$i}sum1/=NF;sum2/=NF;printf "%e\t%e\t",sum1,sqrt((sum2-sum1*sum1)/NF);print $0}' > $TMP
paste <(awk '{print $1}' res.n_stat1) <(awk '{print $0}' $TMP) > res.n_stat_stat
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

## corr

cmd_awk=""
for ((i=1;i<=$NRUN;i++))
do
  cmd_tmp="<(awk 'NR>=3{print \$2}' res.corr$i)"
  cmd_awk="$cmd_awk $cmd_tmp"
done

TMP="tmp.corr_stat"
eval paste $cmd_awk | awk '{sum1=0;sum2=0;for(i=1;i<=NF;i++){sum1+=$i;sum2+=$i*$i}sum1/=NF;sum2/=NF;printf "%e\t%e\t",sum1,sqrt((sum2-sum1*sum1)/NF);print $0}' > $TMP
paste <(awk 'NR>=3{print $1}' res.corr1) <(awk '{print $0}' $TMP) > res.corr_stat
rm $TMP

## face_avg 

#cmd_awk=""
#for ((i=1;i<=$NRUN;i++))
#do
#  cmd_tmp="<(awk 'NR>=3{print \$2}' res.face_avg$i)"
#  cmd_awk="$cmd_awk $cmd_tmp"
#done

#TMP="tmp.face_avg_stat"
#eval paste $cmd_awk | awk '{sum1=0;sum2=0;for(i=1;i<=NF;i++){sum1+=$i;sum2+=$i*$i}sum1/=NF;sum2/=NF;printf "%e\t%e\t",sum1,sqrt((sum2-sum1*sum1)/NF);print $0}' > $TMP
#paste <(awk 'NR>=3{print $1}' res.face_avg1) <(awk '{print $0}' $TMP) > res.face_avg_stat
#rm $TMP

########
# plot #
########

echo "=== running $PYSCR5 ==="
python ../$PYSCR5

