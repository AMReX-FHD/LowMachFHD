# variables "datafile1" and "datafile2" will be provided by -e option when this script is called. (see anim_gif_comp.sh)

set term png
EXT=".png"

set linestyle 1 lt 1 lc rgb "red"
set linestyle 2 lt 1 lc rgb "blue"
set linestyle 3 lt 1 lc rgb "green"
set linestyle 4 lt 2 lc rgb "gray"

set output datafile.EXT

set key reverse Left
plot [0:256][0:2] 0.5 w l ls 4 t "",1.6 w l ls 4 t "",datafile1 u 1:3 w l ls 1 t "",datafile2 u 1:3 w l ls 2 t ""

