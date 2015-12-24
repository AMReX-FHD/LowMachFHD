#!/bin/bash

################################################################################
# input: change accordingly                                                    #
################################################################################
                                          # choose one and comment out the other
#JOB=ONLY_FIRST_FRAME                     #  generate only one png file for test  
JOB=ALL_FRAMES                            #  generate all png files 

SCRIPT=../src_python/movie_thumbnails.py  # main python script filename 
THUMBNAILS_INFO=input_thumbnails_info.py

OUTPNG_DIR=PNGS                           # for png files
OUTPNG_FILENAME=${OUTPNG_DIR}/frame

DELAY=20                                  # for the resulting animiated gif file
OUTGIF_FILENAME=1.gif

################################################################################

# delete previous results (if present)
if [ -d $OUTPNG_DIR ]; then rm -rf $OUTPNG_DIR; fi
if [ -f $OUTGIF_FILENAME ]; then rm $OUTGIF_FILENAME; fi

# create a directory for png files
mkdir ${OUTPNG_DIR}

# for the first option, generate a single png and display it
if [ "$JOB" = ONLY_FIRST_FRAME ]
then
  visit -nowin -cli -s $SCRIPT $THUMBNAILS_INFO ${OUTPNG_FILENAME} $JOB 
  display ${OUTPNG_FILENAME}0000.png &
fi

# for the second option, generate png files and convert them into a gif file.
# then, display it. 
if [ "$JOB" = ALL_FRAMES ]
then
  visit -nowin -cli -s $SCRIPT $THUMBNAILS_INFO ${OUTPNG_FILENAME} $JOB 
  convert -delay $DELAY ${OUTPNG_FILENAME}*.png ${OUTGIF_FILENAME}
  animate ${OUTGIF_FILENAME} &
fi
