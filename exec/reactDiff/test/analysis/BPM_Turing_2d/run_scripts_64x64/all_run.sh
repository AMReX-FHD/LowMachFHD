#!/bin/bash

RUNNAME=cs1_64_ExMidTau_0.05
NRUN=16

PLOTSCR=${RUNNAME}_plot.sh
ANALSCR=${RUNNAME}_analysis.sh

if [ ! -f $PLOTSCR ]
then
  echo "ERROR: $PLOTSCR not found"
  exit
fi

if [ ! -f $ANALSCR ]
then
  echo "EROOR: $ANALSCR not found"
  exit
fi

./$PLOTSCR plot

for ((i=1;i<=NRUN;i++))
do
  ./$ANALSCR $i
done

./fit.sh $RUNNAME $NRUN
