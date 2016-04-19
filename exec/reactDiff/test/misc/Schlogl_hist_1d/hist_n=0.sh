#!/bin/bash

RUNNAME=AVG1_DT0.001_RUN
NRUN=16

PYSCR=hist_n=0.py
OUTPUT_STAT=res.hist_n=0_stat

################################################################################

for ((i=1;i<=NRUN;i++))
do
  dir=$RUNNAME$i
  cd $dir
  echo "***** $dir"
  python ../$PYSCR 
  cd ..
  echo 
done

for ((i=1;i<=NRUN;i++))
do
  dir=$RUNNAME$i
  cp $dir/res.hist_n=0 $RUNNAME/res.hist_n=0_$i
done

cmd_awk=""
for ((i=1;i<=NRUN;i++))
do
  dir=$RUNNAME$i
  cmd_tmp="<(awk '{print \$2}' $dir/res.hist_n=0)"
  cmd_awk="$cmd_awk $cmd_tmp"
done

# 1st column=quantity 2nd=mean 3rd=var (not normalized by NRUN) 4th=standard error (normalized by NRUN)
TMP_OUTPUT_STAT=$RUNNAME/${OUTPUT_STAT}_tmp
eval paste $cmd_awk | awk '{sum1=0;sum2=0;for(i=1;i<=NF;i++){sum1+=$i;sum2+=$i*$i}sum1/=NF;sum2/=NF;print sum1,(sum2-sum1*sum1)*NF/(NF-1),sqrt((sum2-sum1*sum1)/NF),$0}' > $TMP_OUTPUT_STAT
paste <(awk '{print $1}' $dir/res.hist_n=0) <(awk '{print $0}' $TMP_OUTPUT_STAT) > $RUNNAME/$OUTPUT_STAT
rm $TMP_OUTPUT_STAT

awk '{printf("%s\t%g\t%g\n",$1,$2,$4)}' $RUNNAME/$OUTPUT_STAT 
