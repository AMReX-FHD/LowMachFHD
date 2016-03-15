#!/usr/bin/gnuplot

set terminal post enhanced color "Roman, 18"
set output "Sk.eps"

set border 31 front linetype -1 linewidth 2.000

set linestyle 1 lt 1 lc rgb "red" lw 4
set linestyle 2 lt 1 lc rgb "blued" lw 4
set linestyle 3 lt 2 lc rgb "black" lw 4

set key center bottom Left reverse spacing 1.3

set xlabel "x = l k"
set ylabel "Static structure factor"

set yrange [0.5:2.5]
set ytics 0.5,0.5,2.5
set mytics 5

## Schlogl model with k1 = 1.5, k2 = 0.6, k3 = 2.05, k4 = 2.95
nb = 1.           ## equilibrium number density d(nb) = 0 
r = 1.75          # r = -d'(nb)
D0 = 3.55         # D0 = D(nb)
A = D0/nb/r

dt = 0.01
dx = 0.25
rm = 5.71429      # rate multiplier
chi = 1.5625      # diffusion coefficient

alpha = rm*r*dt
beta = chi*dt/dx**2   
l = sqrt(chi/(rm*r))   # length scale

xmax = sqrt(beta/alpha)*3.141592
Sk0(x) = nb*(A+x**2)/(1+x**2)

plot [-xmax:xmax] "Schlogl.S_k.pair=1.Re.dat" u (l*$1):2 w l ls 1 t "Numerical result", "res_theo" u 1:2 w l ls 2 t "Theoretical result", Sk0(x) w l ls 3 t "Theoretical result ({/Symbol=D}x{/Symbol=\256}0, {/Symbol=D}t{/Symbol=\256}0)"

