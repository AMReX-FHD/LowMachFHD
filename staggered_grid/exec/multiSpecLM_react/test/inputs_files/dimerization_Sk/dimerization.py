import numpy as np
import math

# mass density and number density

m = 3.e-23
m1 = m
m2 = 2*m
print "m1 = %e, m2 = %e" % (m1,m2) 

rho = 1.
c0 = 0.5
c1 = c0 
c2 = 1.-c0
print "rhobar1 = rhobar2 = %e" % rho
print "w1 = %e, w2 = %e" % (c1,c2)

rho1 = rho*c1
rho2 = rho*c2
print "rho = %e" % rho
print "rho1 = %e, rho2 = %e" % (rho1,rho2)

n1 = rho1/m1
n2 = rho2/m2
n0 = n1+2*n2
print "n1 = %e, n2 = %e" % (n1,n2)
print "n0 = %e" % n0

# rate const

D = 1.e-5
dx = 1.e-6
dz = 1.
dv = dx**2*dz
print "D = %e" % D
print "dx = %e, dz = %e" % (dx,dz)

ratio1 = math.sqrt(10.)
k2 = D*c0/(2-c0)/(ratio1*dx)**2 
k1 = k2/n0*(1-c0)/2/c0**2 
r = (2-c0)/c0*k2
print "k1 = %e, k2 = %e" % (k1,k2)
print "r = %e, Gamma = " % (r)

pen_dep = math.sqrt(D/r)
print "pen_dep = %e" % pen_dep

# misc

print "N1 = %e, N2 = %e" % (n1*dv,n2*dv)

kB = 1.38e-16
T = 300.
print "kB = %e, T = %e" % (kB,T)

urand = math.sqrt(kB*T/rho/dz)/dx
print "urand = %e" % urand
print "cell Pe = %e" % (urand*dx/D) 

dt = 1.e-8
nu = 1.e-2
print "dt = %e" % dt
print "nu = %e" % nu

print "D*dt/dx^2 = %e" % (D*dt/dx**2)
print "nu*dt/dx^2 = %e" % (nu*dt/dx**2)

###########################
# mole fraction based LMA #
###########################

x1 = 2*c0/(1+c0)
x2 = 1-x1
print "x1 = %e, x2 = %e" % (x1,x2)

ratio1 = math.sqrt(10.)
kx2 = 0.25*rho/m*c0*(1+c0)**2*D/(ratio1*dx)**2
kx1 = 0.25*(1-c0**2)/c0**2*kx2

r = 4*m/rho/c0/(1+c0)**2*kx2
pen_dep = math.sqrt(D/r)

print "kx1 = %e, kx2 = %e" % (kx1,kx2)
print "r = %e, pen_dep = %e" % (r,pen_dep)
print "kx1*x1^2 = %e" % (kx1*x1**2)
print "kx2*x2 = %e" % (kx2*x2)

######################
# distribution of x1 #
######################

Mbar = n1*dv
Nbar = n2*dv
x0 = Mbar/(Mbar+Nbar)
print "dz=%e, Mbar=%e, Nbar=%e, x0=%f" % (dz,Mbar,Nbar,x0)

# dz = 1
xmin = 0.66665
xmax = 0.66668
# dz = 1e-1
#xmin = 0.66662
#xmax = 0.66672
# dz = 1e-2
#xmin = 0.6665
#xmax = 0.6668
# dz = 1e-3
#xmin = 0.666
#xmax = 0.6672
# dz = 1e-4
#xmin = 0.665
#xmax = 0.668

Nbin = 101
output = "theo.hist"

xx = []
pp = []
for i in range(Nbin):
   x = xmin+(xmax-xmin)/(Nbin-1)*i
   A = (1-x)**2/(1-x0)+x**2/x0
   xx.append(x)
   pp.append(math.exp(-0.5*(Mbar+Nbar)*(1-1/A))/math.sqrt(2*math.pi*Mbar*Nbar)*math.sqrt((Mbar+Nbar)/A)**3)

np.savetxt(output,np.c_[xx,pp])
print "** %s generated" % output
