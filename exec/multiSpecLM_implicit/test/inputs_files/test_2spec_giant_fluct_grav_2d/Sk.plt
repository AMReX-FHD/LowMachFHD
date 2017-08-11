set term post enhanced color "Roman, 18"
set output "Sk.eps"

set style line 1 lt 2 lw 2 lc rgb "black" 
set style line 2 lt 1 lw 2 lc rgb "black" 
set style line 3 lt 1 lw 2 lc rgb "red" pt 7 ps 0.8

set key right top reverse Left samplen 2 spacing 1.3
set bmargin 4
set log xy
set format y "%.0e"
set xlabel "~{k}{0.6\\~}"
set ylabel "S_{11}"

kB = 1.; T = 1.; nu = 1.e3; D = 1.

#Nx = 64; Ny = 32; c1 = 0.009; c2 = 0.011; dx = 1; dy = 1
Nx = 128; Ny = 64; c1 = 0.008; c2 = 0.012; dx = 1; dy = 1
#Nx = 128; Ny = 128; c1 = 0.006; c2 = 0.014; dx = 1; dy = 1

g = -1.e6; rhobar1 = 2.; rhobar2 = 1.

kmax = 2/dx
kmin = 2/dx*sin(3.141592/Nx)

c0 = 0.5*(c1+c2)
h = (c1-c2)/(Ny*dy)

rho = 1/(c0/rhobar1+(1-c0)/rhobar2)
beta = rho*(1/rhobar2-1/rhobar1)

plot [kmin:kmax] kB*T/(nu*D*x**4)*h**2 ls 1 t "theory (no gravity)", kB*T/(nu*D*x**4+g*rho*beta*h)*h**2 ls 2 t "theory (gravity)", "diffbarrier_2d_vstat.S_k.pair=2.Re.dat" u (sin($1/2)/(1./2)):2 w lp ls 3 t ""
