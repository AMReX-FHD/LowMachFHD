# variable "datafile" will be provided by -e option when this script is called. (see anim_gif_AB.sh)

set term png
EXT=".png"

set output datafile.EXT

set key reverse Left
plot [:][-2:12] datafile u 1:3 w l t "n_A",'' u 1:4 w l t "n_B"

