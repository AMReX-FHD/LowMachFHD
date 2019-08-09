#!/bin/bash

LASTFRAME=$1

if [ -z "$LASTFRAME" ]
then
  echo "ERROR: need to set \$LASTFRAME"
  exit
fi

for i in $(seq -f "%04g" 0 $LASTFRAME)
do
  echo $i
  convert ../1_crop/Na${i}.png -resize 367x307 Na${i}.png
  convert ../1_crop/Cl${i}.png -resize 367x307 Cl${i}.png
  convert ../1_crop/H${i}.png  -resize 367x307 H${i}.png
  convert ../1_crop/OH${i}.png -resize 367x307 OH${i}.png
done

