#!/bin/bash

SRCGIF=1.gif
DESTGIF=2.gif 

convert $SRCGIF -pointsize 20 \
  -annotate +93+100 "Diffusion Fluctuation: ON" \
  -annotate +391+100 "Diffusion Fluctuation: OFF" \
  -annotate 270x270+45+367 "Reaction Fluctuation: ON" \
  -annotate 270x270+45+675 "Reaction Fluctuation: OFF" \
  -pointsize 30 -annotate +145+45 "Number Density of Species U" \
  $DESTGIF
animate $DESTGIF &
