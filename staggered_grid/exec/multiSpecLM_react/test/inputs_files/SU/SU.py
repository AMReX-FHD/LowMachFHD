import math

# mass density and number density

rho = 1.
m = 3.e-23
print "rhobar1 = rhobar2 = %e" % rho
print "m1 = m2 = %e" % m 

c1 = 0.01
c2 = 1.-c1
print "w1 = %e, w2 = %e" % (c1,c2)

rho1 = rho*c1
rho2 = rho*c2
print "rho = %e" % rho
print "rho1 = %e, rho2 = %e" % (rho1,rho2)

n1 = rho1/m
n2 = rho2/m
print "n1 = %e, n2 = %e" % (n1,n2)

# rate const

D = 1.e-5
dx = 1.e-6
dz = 1.e-5
dv = dx**2*dz
print "D = %e" % D
print "dx = %e, dz = %e" % (dx,dz)

ratio1 = 2.             # k2/k1
ratio2 = math.sqrt(10.) # pen_dep/dx

k1 = D/(ratio1-1)/ratio2**2/dx**2
k2 = ratio1*k1
k3 = n1*(k2-k1)
r = k2-k1
Gamma = k3/2*(1+(k1+k2)/(k2-k1))
print "k1 = %e, k2 = %e, k3 = %e" % (k1,k2,k3)
print "r = %e, Gamma = %e" % (r,Gamma)

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
