# variable "datafile" will be provided by -e option when this script is called. (see anim_gif_Schlogl.sh)

set term png
EXT=".png"

set output datafile.EXT

set key reverse Left
plot [0:128][0:2] datafile u 1:3 w l t ""

