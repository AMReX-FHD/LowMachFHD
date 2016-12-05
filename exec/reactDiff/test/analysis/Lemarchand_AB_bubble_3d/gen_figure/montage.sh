#!/bin/bash

RUNNAME=TEST

convert ${RUNNAME}_0002.png -set caption 't = 50'   tmp1.png 
convert ${RUNNAME}_0004.png -set caption 't = 100'  tmp2.png 
convert ${RUNNAME}_0006.png -set caption 't = 150'  tmp3.png 
convert ${RUNNAME}_0008.png -set caption 't = 200'  tmp4.png 
montage -label '%[caption]' -pointsize 30 tmp1.png tmp2.png tmp3.png tmp4.png colorbar.png -tile x1 -mode Concatenate ${RUNNAME}.png

rm -rf tmp?.png
eog ${RUNNAME}.png &
