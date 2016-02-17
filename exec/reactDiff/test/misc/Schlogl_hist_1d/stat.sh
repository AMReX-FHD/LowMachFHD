#!/bin/bash

RUNNAME=TEST
NRUN=16
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
  cp $RUNNAME$i/res.n_stat $RUNNAME/res.n_stat$i
done

##############
# statistics #
##############

# 1st column=quantity
# 2nd=mean
# 3rd=var (not normalized by NRUN)
# 4th=standard error (normalized by NRUN)

cd $RUNNAME

# OUTPUT1 (hist)

cmd_awk=""
for ((i=1;i<=$NRUN;i++))
do
  cmd_tmp="<(awk 'NR>1{print \$2}' res.hist$i)"
  cmd_awk="$cmd_awk $cmd_tmp"
done

TMP_OUTPUT1=${OUTPUT1}_tmp
eval paste $cmd_awk | awk '{sum1=0;sum2=0;for(i=1;i<=NF;i++){sum1+=$i;sum2+=$i*$i}sum1/=NF;sum2/=NF;print sum1,(sum2-sum1*sum1)*NF/(NF-1),sqrt((sum2-sum1*sum1)/NF),$0}' > $TMP_OUTPUT1
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
eval paste $cmd_awk | awk '{sum1=0;sum2=0;for(i=1;i<=NF;i++){sum1+=$i;sum2+=$i*$i}sum1/=NF;sum2/=NF;print sum1,(sum2-sum1*sum1)*NF/(NF-1),sqrt((sum2-sum1*sum1)/NF),$0}' > $TMP_OUTPUT2
paste <(awk 'NR>2{print $1}' res.Sk1) <(awk '{print $0}' $TMP_OUTPUT2) > $OUTPUT2
rm $TMP_OUTPUT2

# OUTPUT3 (n_stat: mean, var, prob of neg dens, and KL-divergence)

cmd_awk=""
for ((i=1;i<=$NRUN;i++))
do
  cmd_tmp="<(awk '{print \$2}' res.n_stat$i)" 
  cmd_awk="$cmd_awk $cmd_tmp"
done

TMP_OUTPUT3=${OUTPUT3}_tmp
eval paste $cmd_awk | awk '{sum1=0;sum2=0;for(i=1;i<=NF;i++){sum1+=$i;sum2+=$i*$i}sum1/=NF;sum2/=NF;print sum1,(sum2-sum1*sum1)*NF/(NF-1),sqrt((sum2-sum1*sum1)/NF),$0}' > $TMP_OUTPUT3
paste <(awk '{print $1}' res.n_stat1) <(awk '{print $0}' $TMP_OUTPUT3) > $OUTPUT3

awk '{printf("%s\t%g\t%g\n",$1,$2,$4)}' $OUTPUT3
rm $TMP_OUTPUT3

########
# plot #
########

PYSCR_PLOT=stat_plot.py
python ../$PYSCR_PLOT $NRUN
