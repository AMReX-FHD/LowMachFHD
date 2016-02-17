#!/bin/bash

RUNNAME=TEST
NRUN=16

NCELL=64
DV=10.
AVGTYPE=1

PYSCR=face_avg.py
OUTPUT_STAT=res.face_avg_stat

################################################################################

echo "NCELL=$NCELL"
echo "DV=$DV"
echo "AVGTYPE=$AVGTYPE"
echo

for ((i=1;i<=NRUN;i++))
do
  dir=$RUNNAME$i
  cd $dir
  echo "***** $dir"
  python ../$PYSCR $NCELL $DV $AVGTYPE | tee res.face_avg
  cd ..
  echo 
done

for ((i=1;i<=NRUN;i++))
do
  dir=$RUNNAME$i
  cp $dir/res.face_avg $RUNNAME/res.face_avg$i
done

cmd_awk=""
for ((i=1;i<=NRUN;i++))
do
  dir=$RUNNAME$i
  cmd_tmp="<(awk 'NR>2{print \$2}' $dir/res.face_avg)"
  cmd_awk="$cmd_awk $cmd_tmp"
done

# 1st column=quantity 2nd=mean 3rd=var (not normalized by NRUN) 4th=standard error (normalized by NRUN)
TMP_OUTPUT_STAT=$RUNNAME/${OUTPUT_STAT}_tmp
eval paste $cmd_awk | awk '{sum1=0;sum2=0;for(i=1;i<=NF;i++){sum1+=$i;sum2+=$i*$i}sum1/=NF;sum2/=NF;print sum1,(sum2-sum1*sum1)*NF/(NF-1),sqrt((sum2-sum1*sum1)/NF),$0}' > $TMP_OUTPUT_STAT
paste <(awk 'NR>2{print $1}' $dir/res.face_avg) <(awk '{print $0}' $TMP_OUTPUT_STAT) > $RUNNAME/$OUTPUT_STAT
rm $TMP_OUTPUT_STAT

awk '{printf("%s\t%g\t%g\n",$1,$2,$4)}' $RUNNAME/$OUTPUT_STAT 
