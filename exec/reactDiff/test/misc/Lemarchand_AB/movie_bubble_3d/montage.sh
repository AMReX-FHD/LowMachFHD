#!/bin/bash

mkdir tmp 

convert PNG0/frame0050.png -resize 512x500 -set caption 't = 50'  tmp/frame0050.png 
convert PNG0/frame0100.png -resize 512x500 -set caption 't = 100' tmp/frame0100.png
convert PNG0/frame0150.png -resize 512x500 -set caption 't = 150' tmp/frame0150.png
convert PNG0/frame0200.png -resize 512x500 -set caption 't = 200' tmp/frame0200.png
montage -label '%[caption]' -pointsize 30 tmp/frame0050.png tmp/frame0100.png tmp/frame0150.png tmp/frame0200.png colorbar.png -tile x1 -mode Concatenate tile.png

rm -rf tmp
eog tile.png &
