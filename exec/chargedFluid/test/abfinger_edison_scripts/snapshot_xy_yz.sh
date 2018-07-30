#!/bin/bash

RUNNAME=$1

if [ -z "$RUNNAME" ]
then
  echo "ERROR: need to set \$RUNNAME"
  exit
fi

if [ ! -d $RUNNAME ]
then
  echo "$RUNNAME does not exist"
  exit
fi

VISIT=${RUNNAME}.visit
ls $RUNNAME/plt*/Header > $VISIT

NLINE=`wc $VISIT | awk '{print $1}'`
echo "$NLINE frames"
LASTFRAME=$(($NLINE-1))

PNG=${RUNNAME}_png_xy_yz
if [ -d $PNG ]
then
  echo "$PNG already exists"
  exit
else
  mkdir $PNG
fi

visit -cli -nowin -l sbatch/srun -p debug -np 24 -t 00:10:00 -s snapshot_xy_yz.py $VISIT $PNG Na
mv visitlog.py $PNG/visitlog.py_Na

visit -cli -nowin -l sbatch/srun -p debug -np 24 -t 00:10:00 -s snapshot_xy_yz.py $VISIT $PNG Cl 
mv visitlog.py $PNG/visitlog.py_Cl

visit -cli -nowin -l sbatch/srun -p debug -np 24 -t 00:10:00 -s snapshot_xy_yz.py $VISIT $PNG H
mv visitlog.py $PNG/visitlog.py_H

visit -cli -nowin -l sbatch/srun -p debug -np 24 -t 00:10:00 -s snapshot_xy_yz.py $VISIT $PNG OH
mv visitlog.py $PNG/visitlog.py_OH

PNGDIR0=$PNG/0_orig
PNGDIR1=$PNG/1_crop
PNGDIR2=$PNG/2_resize

mkdir $PNGDIR0
mkdir $PNGDIR1
mkdir $PNGDIR2

cp snapshot_xy_yz_crop.sh $PNGDIR1
cp snapshot_xy_yz_resize.sh $PNGDIR2

mv $PNG/*.png $PNGDIR0
echo "snapshot_xy_yz_crop.sh"
cd $PNGDIR1
./snapshot_xy_yz_crop.sh $LASTFRAME
cd ../..
echo "snapshot_xy_yz_resize.sh"
cd $PNGDIR2
./snapshot_xy_yz_resize.sh $LASTFRAME
cd ../..

mv slurm-*.out $PNG
