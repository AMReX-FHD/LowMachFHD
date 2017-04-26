set term post enhanced color "Roman, 18"
set output "Sk.eps"

set style line 11 lt 2 lw 2 lc rgb "red" 
set style line 12 lw 2 lc rgb "red" pt 4 ps 0.8
set style line 21 lt 2 lw 2 lc rgb "green" 
set style line 22 lw 2 lc rgb "green" pt 4 ps 0.8
set style line 31 lt 2 lw 2 lc rgb "blue" 
set style line 32 lw 2 lc rgb "blue" pt 4 ps 0.8

set key left bottom reverse Left samplen 2 spacing 1.3
set bmargin 4
set log xy
set format y "%.0e"
set xlabel "~{k}{0.6\\~}"
set ylabel "{/Symbol=n} S_{ij}"

#nu=1000
#c11=0.4344465040e-11
#c12=0.1844617185e-10
#c22=0.8639076145e-10

#nu=100
#c11=0.4324414486e-10
#c12=0.1830737676e-9
#c22=0.8542476270e-9

#nu=10
#c11=0.4138039156e-9
#c12=0.1712915944e-8
#c22=0.7762870575e-8

nu=1
c11=0.2675864886e-8
c12=0.1018456594e-7
c22=0.4076878832e-7

plot [0.08:3][0.5e-10:0.2e-3] c11/x**4*nu ls 11 t "S_{11}, theory", "diffbarrier_2d_vstat.S_k.pair=2.Re.dat" u (sin($1/2)/(1./2)):($2*nu) w p ls 12 t "S_{11}", c22/x**4*nu ls 21 t "S_{22}, theory", "diffbarrier_2d_vstat.S_k.pair=3.Re.dat" u (sin($1/2)/(1./2)):($2*nu) w p ls 22 t "S_{22}", c12/x**4*nu ls 31 t "S_{12}, theory", "diffbarrier_2d_vstat.S_k.pair=4.Re.dat" u (sin($1/2)/(1./2)):($2*nu) w p ls 32 t "S_{12}"
