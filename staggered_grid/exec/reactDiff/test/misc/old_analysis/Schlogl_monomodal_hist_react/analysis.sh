#!/bin/bash

RUNNAME=prev_ImMidTau_dt1
NRUN=16
SKIPLINE=10000
NAV=1.
NCELL=4
DV=5.

if [ -d $RUNNAME ]
then
  echo "$RUNNAME already exists..."
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
  #####################
  # enter RUN directory
  #####################

  RUNDIR=${RUNNAME}_RUN$i
  echo $RUNDIR
  cd $RUNDIR

  cp scr_out ../$RUNNAME/scr_out$i

  #####################
  # discard SKIPLINE lines from the beginning of fort.9
  #####################

  LINE1=`wc fort.10 | awk '{print $1}'`
  LINE2=`wc fort.9  | awk '{print $1}'`
  DIV=`echo print "${LINE1}./($LINE2-$SKIPLINE)" | python`

  if [ $DIV != "$NCELL.0" ]
  then
    echo "ERROR: check SKIPLINE $LINE1 / ( $LINE2 - $SKIPLINE ) = $DIV"
    exit
  fi

  LINE3=`echo print "$LINE2-$SKIPLINE" | python`
  tail --lines $LINE3 fort.9 > fort.9_trunc

  #####################
  # for fort.10, run scripts hist_n.py and hist_n_near_zero.py and copy results
  #####################

  python ../hist_n_cell.py $NAV $DV no no
  mv cell.hist ../$RUNNAME/cell.hist$i
  mv cell.hist_poiss ../$RUNNAME/cell.hist_poiss   # same for different i 
  mv cell.hist_cont ../$RUNNAME/cell.hist_cont     # same for different i
  mv cell.n_stat ../$RUNNAME/cell.n_stat$i
  #mv cell_hist_linear.png ../$RUNNAME/cell_hist_linear$i.png
  #mv cell_hist_semilogy.png ../$RUNNAME/cell_hist_semilogy$i.png

  python ../hist_n_near_zero_cell.py
  mv cell.hist_near_zero ../$RUNNAME/cell.hist_near_zero$i

  #####################
  # for fort.9_trunc, run scripts hist_n.py and hist_n_near_zero.py and copy results
  #####################

  DVTOT=`echo print "$DV*$NCELL" | python`

  python ../hist_n_tot.py $NAV $DVTOT no no
  mv tot.hist ../$RUNNAME/tot.hist$i
  mv tot.hist_poiss ../$RUNNAME/tot.hist_poiss   # same for different i 
  mv tot.hist_cont ../$RUNNAME/tot.hist_cont     # same for different i
  mv tot.n_stat ../$RUNNAME/tot.n_stat$i
  #mv tot_hist_linear.png ../$RUNNAME/tot_hist_linear$i.png
  #mv tot_hist_semilogy.png ../$RUNNAME/tot_hist_semilogy$i.png

  python ../hist_n_near_zero_tot.py
  mv tot.hist_near_zero ../$RUNNAME/tot.hist_near_zero$i

  #####################
  # return to pwd 
  #####################

  echo
  cd ..

  mv $RUNDIR $RUNNAME
done

#########################################################
# statistics (sample mean and variance)                 #
# 1st col = quantity_name, k for S(k), or n for hist    #
# 2nd col = sample mean                                 #
# 3rd col = std error (normalized by NRUN               #
# additional NRUN columns from each sample run          #
#########################################################

cd $RUNNAME

#################
# stat for cell #
#################

## hist

cmd_awk=""
for ((i=1;i<=$NRUN;i++))
do
  cmd_tmp="<(awk 'NR>=2{print \$2}' cell.hist$i)"
  cmd_awk="$cmd_awk $cmd_tmp"
done

TMP="tmp.hist_stat"
eval paste $cmd_awk | awk '{sum1=0;sum2=0;for(i=1;i<=NF;i++){sum1+=$i;sum2+=$i*$i}sum1/=NF;sum2/=NF;printf "%e\t%e\t",sum1,sqrt((sum2-sum1*sum1)/NF);print $0}' > $TMP
paste <(awk 'NR>=2{print $1}' cell.hist1) <(awk '{print $0}' $TMP) > cell.hist_stat
rm $TMP

## n_stat

cmd_awk=""
for ((i=1;i<=$NRUN;i++))
do
  cmd_tmp="<(awk '{print \$2}' cell.n_stat$i)"
  cmd_awk="$cmd_awk $cmd_tmp"
done

TMP="tmp.n_stat_stat"
eval paste $cmd_awk | awk '{sum1=0;sum2=0;for(i=1;i<=NF;i++){sum1+=$i;sum2+=$i*$i}sum1/=NF;sum2/=NF;printf "%e\t%e\t",sum1,sqrt((sum2-sum1*sum1)/NF);print $0}' > $TMP
paste <(awk '{print $1}' cell.n_stat1) <(awk '{print $0}' $TMP) > cell.n_stat_stat
rm $TMP

## hist_near_zero

cmd_awk=""
for ((i=1;i<=$NRUN;i++))
do
  cmd_tmp="<(awk '{print \$2}' cell.hist_near_zero$i)"
  cmd_awk="$cmd_awk $cmd_tmp"
done

TMP="tmp.hist_near_zero_stat"
eval paste $cmd_awk | awk '{sum1=0;sum2=0;for(i=1;i<=NF;i++){sum1+=$i;sum2+=$i*$i}sum1/=NF;sum2/=NF;printf "%e\t%e\t",sum1,sqrt((sum2-sum1*sum1)/NF);print $0}' > $TMP
paste <(awk '{print $1}' cell.hist_near_zero1) <(awk '{print $0}' $TMP) > cell.hist_near_zero_stat
rm $TMP

# plots
python ../stat_plot_cell.py

################
# stat for tot #
################

## hist

cmd_awk=""
for ((i=1;i<=$NRUN;i++))
do
  cmd_tmp="<(awk 'NR>=2{print \$2}' tot.hist$i)"
  cmd_awk="$cmd_awk $cmd_tmp"
done

TMP="tmp.hist_stat"
eval paste $cmd_awk | awk '{sum1=0;sum2=0;for(i=1;i<=NF;i++){sum1+=$i;sum2+=$i*$i}sum1/=NF;sum2/=NF;printf "%e\t%e\t",sum1,sqrt((sum2-sum1*sum1)/NF);print $0}' > $TMP
paste <(awk 'NR>=2{print $1}' tot.hist1) <(awk '{print $0}' $TMP) > tot.hist_stat
rm $TMP

## n_stat

cmd_awk=""
for ((i=1;i<=$NRUN;i++))
do
  cmd_tmp="<(awk '{print \$2}' tot.n_stat$i)"
  cmd_awk="$cmd_awk $cmd_tmp"
done

TMP="tmp.n_stat_stat"
eval paste $cmd_awk | awk '{sum1=0;sum2=0;for(i=1;i<=NF;i++){sum1+=$i;sum2+=$i*$i}sum1/=NF;sum2/=NF;printf "%e\t%e\t",sum1,sqrt((sum2-sum1*sum1)/NF);print $0}' > $TMP
paste <(awk '{print $1}' tot.n_stat1) <(awk '{print $0}' $TMP) > tot.n_stat_stat
rm $TMP

## hist_near_zero

cmd_awk=""
for ((i=1;i<=$NRUN;i++))
do
  cmd_tmp="<(awk '{print \$2}' tot.hist_near_zero$i)"
  cmd_awk="$cmd_awk $cmd_tmp"
done

TMP="tmp.hist_near_zero_stat"
eval paste $cmd_awk | awk '{sum1=0;sum2=0;for(i=1;i<=NF;i++){sum1+=$i;sum2+=$i*$i}sum1/=NF;sum2/=NF;printf "%e\t%e\t",sum1,sqrt((sum2-sum1*sum1)/NF);print $0}' > $TMP
paste <(awk '{print $1}' tot.hist_near_zero1) <(awk '{print $0}' $TMP) > tot.hist_near_zero_stat
rm $TMP

# plots
python ../stat_plot_tot.py
