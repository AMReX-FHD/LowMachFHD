# mole fraction based LMA
gnuplot> set format y "%e"; dx = 1.e-6; k2 = 3.3333333e5; c0 = 0.5; r = (2-c0)/c0*k2; D = 1.e-5; d = sqrt(D/r); n0 = 3.3333333e22; plot [-2/dx:2/dx] "res.Sk" u (sin(0.5*dx*$1)/(0.5*dx)):2 w l, 3/(8*n0) 

# number density based LMA
gnuplot> set format y "%e"; dx = 1.e-6; k2 = 3.3333333e5; c0 = 0.5; r = (2-c0)/c0*k2; D = 1.e-5; d = sqrt(D/r); n0 = 3.3333333e22; plot [-2/dx:2/dx] "res.Sk" u (sin(0.5*dx*$1)/(0.5*dx)):2 w l, 3/(8*n0)*(8./9+(x*d)**2)/(1+(x*d)**2) 
