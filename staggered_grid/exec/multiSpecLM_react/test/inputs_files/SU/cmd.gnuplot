gnuplot> set format y "%e"; m=3e-23; nbar=3.3333333e20; r=1e6; Gamma=6.6666667e26; D=1e-5; dx=1e-6; plot [-2/dx:2/dx] "res.Sk" u (sin(0.5*dx*$1)/(0.5*dx)):2 w l,(D*nbar*x**2+Gamma)/(D*x**2+r)*m**2
