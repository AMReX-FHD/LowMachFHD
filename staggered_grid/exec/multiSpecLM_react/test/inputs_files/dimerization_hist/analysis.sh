#!/bin/bash

RUNNAME=TEST
NRUN=16
#PYARGS="../../../analysis/hist/hist_N.py 2 2. 5 0 10 no yes"
PYARGS="../../../analysis/hist/hist_N.py 2 2. 10 0 20 no yes"
#PYARGS="../../../analysis/hist/hist_N.py 2 2. 15 0 30 no yes"

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

  # run script hist_n.py and copy results
  python $PYARGS  
  cp res.hist ../$RUNNAME/res.hist$i
  cp res.N_stat ../$RUNNAME/res.N_stat$i

  # extract a vtk file of S(k) for sample average
  python ../extract_Sk_vtk.py
  cp Sk_vtk_out0 ../$RUNNAME/Sk_vtk_out0
  cp Sk_vtk_out1 ../$RUNNAME/Sk_vtk_out1.$i
  cp Sk_vtk_out2 ../$RUNNAME/Sk_vtk_out2

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

TMP="tmp.N_stat_stat"
eval paste $cmd_awk | awk '{sum1=0;sum2=0;for(i=1;i<=NF;i++){sum1+=$i;sum2+=$i*$i}sum1/=NF;sum2/=NF;printf "%e\t%e\t",sum1,sqrt((sum2-sum1*sum1)/NF);print $0}' > $TMP
paste <(awk '{print $1}' res.N_stat1) <(awk '{print $0}' $TMP) > res.N_stat_stat
rm $TMP

## Sk vtk

cmd_awk=""
for ((i=1;i<=$NRUN;i++))
do
  cmd_tmp="Sk_vtk_out1.$i"
  cmd_awk="$cmd_awk $cmd_tmp"
done

TMP="Sk_vtk_out1_stat"
FINAL="Sk_stat.vtk"
eval paste $cmd_awk | awk '{sum=0;for(i=1;i<=NF;i++){sum+=$i}sum/=NF;printf "%.12e\n",sum}' > $TMP
cat Sk_vtk_out0 > $FINAL
cat $TMP >> $FINAL
cat Sk_vtk_out2 >> $FINAL
