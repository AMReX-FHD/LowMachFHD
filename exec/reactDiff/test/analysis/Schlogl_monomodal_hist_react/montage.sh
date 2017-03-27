#!/bin/bash

RUNNAME=prev_ImMidTau_dt1

cd $RUNNAME

montage -tile 2x3 -geometry 450x300 cell_hist_linear.png tot_hist_linear.png cell_hist_semilogy.png tot_hist_semilogy.png cell_hist_near_zero.png tot_hist_near_zero.png res.png

eog res.png
