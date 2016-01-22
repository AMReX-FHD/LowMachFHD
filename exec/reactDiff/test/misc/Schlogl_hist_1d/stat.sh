#!/bin/bash

RUNNAME=TEST
NRUN=4
OUTPUT1=res.hist_stat
OUTPUT2=res.Sk_stat
OUTPUT3=res.n_stat

############################
# copy files from each run #
############################

if [ -d $RUNNAME ]
then
  echo "** Error: directory $RUNNAME already present"
  exit
else
  mkdir $RUNNAME
  cp ${RUNNAME}1/res.hist_poiss $RUNNAME
  cp ${RUNNAME}1/res.hist_cont $RUNNAME
fi

for ((i=1;i<=$NRUN;i++))
do
  cp $RUNNAME$i/screen_out $RUNNAME/screen_out$i
  cp $RUNNAME$i/res.hist $RUNNAME/res.hist$i
  cp $RUNNAME$i/Hist1D.S_k.pair\=1.Re.dat $RUNNAME/res.Sk$i
  cp $RUNNAME$i/res.sum_Sk $RUNNAME/res.sum_Sk$i
done

##############
# statistics #
##############

cd $RUNNAME

# OUTPUT1 (hist)

cmd_awk=""
for ((i=1;i<=$NRUN;i++))
do
  cmd_tmp="<(awk 'NR>1{print \$2}' res.hist$i)"
  cmd_awk="$cmd_awk $cmd_tmp"
done

TMP_OUTPUT1=${OUTPUT1}_tmp
eval paste $cmd_awk | awk '{sum1=0;sum2=0;for(i=1;i<=NF;i++){sum1+=$i;sum2+=$i*$i}sum1/=NF;sum2/=NF;print sum1,(sum2-sum1*sum1)*NF/(NF-1),$0}' > $TMP_OUTPUT1
paste <(awk 'NR>1{print $1}' res.hist1) <(awk '{print $0}' $TMP_OUTPUT1) > $OUTPUT1
rm $TMP_OUTPUT1

# OUTPUT2 (Sk)

cmd_awk=""
for ((i=1;i<=$NRUN;i++))
do
  cmd_tmp="<(awk 'NR>2{print \$2}' res.Sk$i)"
  cmd_awk="$cmd_awk $cmd_tmp"
done

TMP_OUTPUT2=${OUTPUT2}_tmp
eval paste $cmd_awk | awk '{sum1=0;sum2=0;for(i=1;i<=NF;i++){sum1+=$i;sum2+=$i*$i}sum1/=NF;sum2/=NF;print sum1,(sum2-sum1*sum1)*NF/(NF-1),$0}' > $TMP_OUTPUT2
paste <(awk 'NR>2{print $1}' res.Sk1) <(awk '{print $0}' $TMP_OUTPUT2) > $OUTPUT2
rm $TMP_OUTPUT2

# OUTPUT3 (mean and variance of n)

cmd_awk1=""
cmd_awk2=""
for ((i=1;i<=$NRUN;i++))
do
  cmd_tmp1="<(awk 'NR==3{print \$2}' res.sum_Sk$i)"  # mean of n
  cmd_tmp2="<(awk 'NR==5{print \$2}' res.sum_Sk$i)"  # variance of n
  cmd_awk1="$cmd_awk1 $cmd_tmp1"
  cmd_awk2="$cmd_awk2 $cmd_tmp2"
done

eval paste $cmd_awk1 | awk '{sum1=0;sum2=0;for(i=1;i<=NF;i++){sum1+=$i;sum2+=$i*$i}sum1/=NF;sum2/=NF;printf("n_mean= %g (std= %g)\n",sum1,sqrt((sum2-sum1*sum1)*NF/(NF-1)/NF))}' | tee $OUTPUT3
eval paste $cmd_awk2 | awk '{sum1=0;sum2=0;for(i=1;i<=NF;i++){sum1+=$i;sum2+=$i*$i}sum1/=NF;sum2/=NF;printf("n_var=  %g (std= %g)\n",sum1,sqrt((sum2-sum1*sum1)*NF/(NF-1)/NF))}' | tee -a $OUTPUT3

########
# plot #
########

PYSCR_PLOT=stat_plot.py
python ../$PYSCR_PLOT $NRUN
