#!/bin/bash

RUNNAME=TEST
NRUN=500

NC=64
NC21=$(($NC/2+1))
NC22=$(($NC/2+2))

if [ -d $RUNNAME ]
then
  echo "ERROR: $RUNNAME already exists"
  exit
fi

mkdir $RUNNAME

for ((i=1;i<=$NRUN;i++))
do
  RUNNAMENO=${RUNNAME}_RUN$i

  if [ ! -d $RUNNAMENO ]
  then
    echo "ERROR: $RUNNAMENO does not exist"
    exit
  fi

  echo $RUNNAMENO
  cd $RUNNAMENO
  awk -v nc21=$NC21 -v nc22=$NC22 '{if (NR>2 && NR!=nc22) printf "%e\t%e\n",$1,$nc21}' SU.S_k.pair=1.Re.dat > ../$RUNNAME/res.Sk$i
  cd ..

  mv $RUNNAMENO $RUNNAME
done

cd $RUNNAME

cmd_awk=""
for ((i=1;i<=$NRUN;i++))
do
  cmd_tmp="<(awk '{print \$2}' res.Sk$i)"
  cmd_awk="$cmd_awk $cmd_tmp"
done

TMP="tmp.Sk_stat"
eval paste $cmd_awk | awk '{sum1=0;sum2=0;for(i=1;i<=NF;i++){sum1+=$i;sum2+=$i*$i}sum1/=NF;sum2/=NF;printf "%e\t%e\n",sum1,sqrt((sum2-sum1*sum1)/NF)}' > $TMP
paste <(awk '{print $1}' res.Sk1) <(awk '{print $0}' $TMP) > res.Sk_stat
rm $TMP

