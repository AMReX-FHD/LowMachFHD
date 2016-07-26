RUNNAME=cs1_64_ExMidTau_0.05

FITDATDIR=fit_data
SCRDIR=run_scripts

if [ ! -d $FITDATDIR ]
then
  echo "ERROR: $FITDATDIR not found"
  exit
fi

if [ ! -d $SCRDIR ]
then
  echo "ERROR: $SCRDIR not found"
  exit
fi

if [ -d $RUNNAME ]
then
  echo "ERROR: $RUNNAME already exists..."
  exit
fi

mkdir $RUNNAME

cp ${RUNNAME}.fit $RUNNAME
mv ${RUNNAME}.fit $FITDATDIR

mv ${RUNNAME}_RUN* $RUNNAME

mv ${RUNNAME}_plot.sh $SCRDIR
mv ${RUNNAME}_analysis.sh $SCRDIR

