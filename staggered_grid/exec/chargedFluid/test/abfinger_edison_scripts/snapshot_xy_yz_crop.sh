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
  convert ../0_orig/Na${i}.png -crop 1101x921+23+60 +repage Na${i}.png
  convert ../0_orig/Cl${i}.png -crop 1101x921+23+60 +repage Cl${i}.png
  convert ../0_orig/H${i}.png  -crop 1101x921+23+60 +repage H${i}.png
  convert ../0_orig/OH${i}.png -crop 1101x921+23+60 +repage OH${i}.png
done

