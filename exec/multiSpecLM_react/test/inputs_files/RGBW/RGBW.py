import math
import numpy as np
import numpy.linalg

# assume cgs units

# cell size
dx = 1.e-6
dz = 1.e-2 
dV = dx**2*dz
print "dx = dy = %e, dz = %e, dV = %e" % (dx,dz,dV)

# assume rhobar = rho for all species (Boussinesq approximation)
rho = 1.
print "rhobarR = rhobarG = rhobarB = rhobarW = %e, rhoTot = %e" % (rho,rho)

# molmass
mR = 3.e-23
mG = 6.e-23
mB = 3.e-23
mW = 3.e-23
print "mR = %e, mG = %e, mB = %e, mW = %e" % (mR,mG,mB,mW)

# rough estimate for NTot
Napprox = rho/mW*dV
print "total number of molecules per cell: %e (approx)" % Napprox 

# physical parameters
kT = 4.e-14     # approx 300K
D = 1.e-5
nu = 0.01
print "kT = %e" % kT
print "diff coeff: DR = DG = DB = %e" % D
print "nu = %e" % nu

# random advection
uRand = math.sqrt(kT/mW/Napprox)
print "magnitude of random advection: %e" % uRand

# cell Peclet number 
print "cell Peclet number: %e" % (uRand*dx/D)

# ratio1 = wR+wG+wB, ratio2 = wR/(wR+wG+wB)
ratio1 = 3.e-7
ratio2 = 0.1 
print "ratio1 = %e, ratio2 = %e" % (ratio1,ratio2)

wR = ratio1*ratio2
wG = ratio1*(1-2*ratio2)
wB = ratio1*ratio2
wW = 1-ratio1

nR = rho*wR/mR
nG = rho*wG/mG
nB = rho*wB/mB
nW = rho*wW/mW

print "wR: %e\tnR: %e\tNR: %e" % (wR,nR,nR*dV)
print "wG: %e\tnG: %e\tNG: %e" % (wG,nG,nG*dV)
print "wB: %e\tnB: %e\tNB: %e" % (wB,nB,nB*dV)
print "wW: %e\tnW: %e\tNW: %e" % (wW,nW,nW*dV)
print "wTot: %e\tnTot:%e\tNTot: %e" % (wR+wG+wB+wW,nR+nG+nB+nW,(nR+nG+nB+nW)*dV)

# rate constant
kd = 1.e5
ka = kd*nG/nR/nB
print "rate const: kd = %e, ka = %e" % (kd,ka)

# linearized reaction rate and penetration length
r = ka*(nR+nB)+kd
print "linearized reaction rate = %e" % r
print "penetration length / dx = %e" % (math.sqrt(D/r)/dx)

# reaction counts
dt = 1.e-8
print "dt = %e" % dt
print "fwd reactions per time = %e" % (ka*nR*nB*dV*dt)
print "bwd reactions per time = %e" % (kd*nG*dV*dt)

# cfl numbers
print "D*dt/dx^2 = %e" % (D*dt/dx**2)
print "nu*dt/dx^2 = %e" % (nu*dt/dx**2)

################################################################################
# structure factor
################################################################################

ww = np.array([wR,wG,wB,wW])            # mass fraction vector
mm = np.array([mR,mG,mB,mW])            # molecular mass vector
mbar = 1./sum(ww/mm)                    # mean molecular mass
xx = mbar*(ww/mm)                       # mole fraction vector

W = np.diag(ww)
X = np.diag(xx)

wwww = np.outer(ww,ww)
xxxx = np.outer(xx,xx)

one = np.ones(4)
oneone = np.outer(one,one)

tmp = np.linalg.inv(X-xxxx+oneone)
Sw = mbar/rho*np.dot(np.dot(W-wwww,tmp),W-wwww)

print rho*rho*Sw
