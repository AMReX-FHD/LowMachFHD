penlen=sqrt(10)

res="ImMidTau_b5_a0.5/res.Sk_xi"
theo="theo_data/theo_ImMidTau_a0.5_b5.dat"

set term post enhanced color "Times, 20"
set output "Sk.eps"

set key reverse Left
set xlabel "~{k}{.6\\~} l"
set ylabel "Structure Factor"
set bmargin 4 

plot [0:2*penlen][0.8:2.2] res u (sin($1/2/penlen)*2*penlen):2 w l lw 2 t "numerics",'' u 1:3 w l lw 2 t "exact", theo u (sin($1/2/penlen)*2*penlen):2 w l lw 2 t "theory"
