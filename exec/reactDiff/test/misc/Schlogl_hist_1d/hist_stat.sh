#!/bin/bash

RUNNAME=TEST
NRUN=8
OUTPUT=res.hist_stat

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
  cp $RUNNAME$i/res.hist $RUNNAME/res.hist$i
  cp $RUNNAME$i/screen_out $RUNNAME/screen_out$i
done

##############
# statistics #
##############

cd $RUNNAME

cmd_awk=""
for ((i=1;i<=$NRUN;i++))
do
  cmd_tmp="<(awk 'NR>1{print \$2}' res.hist$i)"
  cmd_awk="$cmd_awk $cmd_tmp"
done

TMP_OUTPUT=${OUTPUT}_tmp
eval paste $cmd_awk | awk '{sum1=0;sum2=0;for(i=1;i<=NF;i++){sum1+=$i;sum2+=$i*$i}sum1/=NF;sum2/=NF;print sum1,(sum2-sum1*sum1)*NF/(NF-1),$0}' > $TMP_OUTPUT
paste <(awk 'NR>1{print $1}' res.hist1) <(awk '{print $0}' $TMP_OUTPUT) > $OUTPUT 
rm $TMP_OUTPUT 

########
# plot #
########

#python ../hist_stat_plot.py 
python ../hist_stat_plot.py $NRUN
