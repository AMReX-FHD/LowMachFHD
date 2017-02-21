gnuplot> dx=1e-6;dz=1e-5; plot "res.n_avg" u 0:($1*dx*dx*dz) w l
gnuplot> dx=1e-6;dz=1e-5; plot "res.n_avg" u 0:($2*dx*dx*dz) w l

gnuplot> nu=1.e19;k2k1ratio=2.;pendxratiosq=10.;dx=1.e-6;dk=(2*3.14159265/64/dx);plot "SU.S_k.pair=1.Re.dat" u (sin($1*dx/2)/(dx/2)):($33) w l,nu*(1+1./(k2k1ratio-1.)/(1+pendxratiosq*(dx*x)**2)),'' u (sin($1*dx/2)/(dx/2)):($29) w l,nu*(1+1./(k2k1ratio-1.)/(1+pendxratiosq*(dx*x)**2+pendxratiosq*(dx*dk*4)**2))
