# cgs units
# 1 = HCl, 2 = NaOH, 3 = NaCl, 4 = H2O

m1 = 6.06e-23 
m2 = 6.64e-23
m3 = 9.71e-23
m4 = 2.99e-23 

R = 6.02e23
print "molmass (HCl,NaOH,NaCl,H2O) = (%.2e,%.2e,%.2e,%.2e)" % (m1,m2,m3,m4)
print "mol weight (HCl,NaOH,NaCl,H2O) = (%.2f,%.2f,%.2f,%.2f)" % (R*m1,R*m2,R*m3,R*m4)
print "m1+m2-m3-m4 = %e" % (m1+m2-m3-m4)

# rho = rho0*(1+a1*c1+a2*c2+a3*c3) [a]=cm3/mol [c]=mol/cm3
# from de Wit's paper
a1 = 18.
a2 = 44.
a3 = 41.

# low-Mach EOS: 1 = rho1/rhobar1 + rho2/rhobar2 + rho3/rhobar3 + rho4/rhobar4
# -> rhobar4 = rho0
# -> rhobari = rho0/(1-rho0*ai/R/mi)
rho0 = 1
rhobar1 = rho0/(1-rho0*a1/R/m1) 
rhobar2 = rho0/(1-rho0*a2/R/m2) 
rhobar3 = rho0/(1-rho0*a3/R/m3)
rhobar4 = rho0
print "rhobar = (%g,%g,%g,%g)" % (rhobar1,rhobar2,rhobar3,rhobar4) 

# from de Wit's paper
# 1M HCl rho = 1.0171
w1 = 36.46/1017.1 
w2 = 0.   
w3 = 0.
w4 = 1-w1 

tmp = w1/rhobar1+w2/rhobar2+w3/rhobar3+w4/rhobar4
print "(w1,w2,w3,w4) = (%f,%f,%f,%f), rho = %f" % (w1,w2,w3,w4,1/tmp)

# from de Wit's paper
# 0.4M NaOH rho = 1.0167 (for double diffusion)
# 0.5M NaOH rho = 1.0208 (for Rayleigh-Taylor)
w1 = 0.
w2 = 40.00*0.4/1016.7  
#w2 = 40.00*0.5/1020.8   
w3 = 0.
w4 = 1-w2 

tmp = w1/rhobar1+w2/rhobar2+w3/rhobar3+w4/rhobar4
print "(w1,w2,w3,w4) = (%f,%f,%f,%f), rho = %f" % (w1,w2,w3,w4,1/tmp)
