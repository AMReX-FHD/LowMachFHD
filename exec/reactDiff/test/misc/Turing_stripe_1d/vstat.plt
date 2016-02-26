set term png
EXT=".png"

set output datafile.EXT

set key reverse Left
plot [:][-2:12] datafile u 1:3 w l t "n_1",'' u 1:4 w l t "n_2"

