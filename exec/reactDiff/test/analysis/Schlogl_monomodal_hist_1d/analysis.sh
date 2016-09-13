#!/bin/bash

RUNNAME=TEST
NRUN=4

NAV=1.
DV=5.
NCELL=64
AVG=1

if [ -d $RUNNAME ]
then
  echo "$RUNNAME already exists..."
  sleep 3s
else
  mkdir $RUNNAME
fi

########################################
# run analysis scripts for each RUNDIR #
# copy results to RUNNAME              #
########################################

for ((i=1;i<=$NRUN;i++))
do
  RUNDIR=${RUNNAME}_RUN$i
  echo $RUNDIR
  cd $RUNDIR

  # copy scr_out
  cp scr_out ../$RUNNAME/scr_out$i

  # copy S(k) 
  cp Hist1D.S_k.pair\=1.Re.dat ../$RUNNAME/res.Sk$i

  # run script hist_n.py and copy results
  python ../../../hist_n.py $NAV $DV no yes
  cp res.hist ../$RUNNAME/res.hist$i
  cp res.n_stat ../$RUNNAME/res.n_stat$i
  if [ $i -eq 1 ]
  then
    cp res.hist_poiss ../$RUNNAME/res.hist_poiss
    cp res.hist_cont ../$RUNNAME/res.hist_cont
  fi

  # run script hist_n_near_zero.py and copy results
  python ../hist_n_near_zero.py
  cp res.hist_near_zero ../$RUNNAME/res.hist_near_zero$i

#  # run script face_avg.py and copy results
#  python ../face_avg.py $NCELL $DV $AVG > res.face_avg
#  cp res.face_avg ../$RUNNAME/res.face_avg$i

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

python ../stat_plot.py

