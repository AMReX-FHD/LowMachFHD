#!/bin/bash

SRCGIF=1.gif
DESTGIF=2.gif 

convert $SRCGIF -pointsize 10 \
  -annotate +20+35 "Particle 256x256 dt=0.05" \
  -annotate +168+35 "Particle 128x128 dt=0.2" \
  -annotate +20+183 "RDME 256x256 dt=0.05" \
  -annotate +168+183 "RDME 128x128 dt=0.2" \
  -annotate +20+331 "ExMid 256x256 dt=0.025" \
  -annotate +168+331 "ExMid 128x128 dt=0.1" \
  -annotate +20+479 "ImMid 256x256 dt=0.25" \
  -annotate +168+479 "ImMid 128x128 dt=0.5" \
  -annotate +20+627 "ImMid 256x256 dt=0.5" \
  -annotate +168+627 "ImMid 128x128 dt=1" \
  -pointsize 15 -annotate +50+15 "Number Density of Species W" \
  $DESTGIF
eog $DESTGIF &
