#!/bin/bash

dirlist=`ls -d 128x128_* 256x256_*`

for dir in $dirlist
do
  cd $dir

  if [ -d plt00000000 ]
  then
    echo "ls `pwd`/plt*/Header > movie.visit"
    ls `pwd`/plt*/Header > movie.visit
  else
    echo "ls Turing2D.000*.vtk > movie.visit"
    ls *.vtk > movie.visit
  fi

  cd ..
done
