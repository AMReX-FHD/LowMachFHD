** 64x64, dx=1e-6
gnuplot> D=1e-6;rm=1;pendepth=sqrt(D/(1e5*rm));alpha=1e-2;S0=1.e23*alpha;dx=1e-6;k2k1ratio=2;set format y "%.1e";set samples 1000;plot S0*(1+1./(k2k1ratio-1.)/(1+(pendepth*x)**2)),"SU.S_k.pair=1.Re.dat" u (sin($1*dx/2)/(dx/2)):($33) w l
awk '{if (NR>2 && NR!=34) printf "%e\t%e\n",$1,$33}' SU.S_k.pair=1.Re.dat

** 128x128, dx=0.5e-6
gnuplot> D=1e-6;rm=1;pendepth=sqrt(D/(1e5*rm));alpha=1e-2;S0=1.e23*alpha;dx=0.5e-6;k2k1ratio=2;set format y "%.1e";set samples 1000;plot S0*(1+1./(k2k1ratio-1.)/(1+(pendepth*x)**2)),"SU.S_k.pair=1.Re.dat" u (sin($1*dx/2)/(dx/2)):($65) w l
awk '{if (NR>2 && NR!=66) printf "%e\t%e\n",$1,$65}' SU.S_k.pair=1.Re.dat

gnuplot> D=1e-6;rm=1;pendepth=sqrt(D/(1e5*rm));alpha=1e-2;S0=1.e23*alpha;dx=1e-6;k2k1ratio=2;set format y "%.1e";set samples 1000;plot S0*(1+1./(k2k1ratio-1.)/(1+(pendepth*x)**2)),"RUN_TEST/SU.S_k.pair=1.Re.dat" u (sin($1*dx/2)/(dx/2)):($33) w l,"RUN_0.5dx_4dz/SU.S_k.pair=1.Re.dat" u (sin($1*0.5*dx/2)/(0.5*dx/2)):($65) w l
