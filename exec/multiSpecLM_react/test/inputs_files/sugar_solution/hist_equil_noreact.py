import math
import numpy as np
import numpy.linalg

# species 1: glucose  (solute; R)
# species 2: sucrose  (solute; G)
# species 3: fructose (solute; B)
# species 4: water    (solvent; W)
#
# values are given in cgs units

print
print "#########################"
print "# determine composition #"
print "#########################"

# cell size
dx = 1.e-5
dz = 1.e-5
dV = dx**2*dz
print "dx = dy = %e, dz = %e, dV = %e" % (dx,dz,dV)

# Avogadro number
NA = 6.022141e23

# mass of a single molecule
mR = 180.1559/NA 
mG = 342.2965/NA 
mB = 180.1559/NA
mW = 18.01528/NA
print "mR = %e, mG = %e, mB = %e, mW = %e" % (mR,mG,mB,mW)

# mean number of molecules per cell
NR = 10 
NG = 5
NB = 0

# number density
nR = NR/dV
nG = NG/dV
nB = NB/dV

# mass density
rhoR = mR*nR
rhoG = mG*nG
rhoB = mB*nB

print "NR = %e, NG = %e, NB = %e" % (NR,NG,NB)
print "nR = %e, nG = %e, nB = %e" % (nR,nG,nB)
print "cR = %e, cG = %e, cB = %e" % (1e3*nR/NA,1e3*nG/NA,1e3*nB/NA)
print "rhoR = %e, rhoG = %e, rhoB = %e" % (rhoR,rhoG,rhoB)

# assume rhobar = rho for all species (Boussinesq approximation)
rho = 1.
print "rhobarR = rhobarG = rhobarB = rhobarW = %e, rhoTot = %e" % (rho,rho)

# W molecules
rhoW = rho-rhoR-rhoG-rhoB
nW = rhoW/mW
NW = nW*dV
print "rhoW = %.16e, nW = %e, NW = %e" % (rhoW,nW,NW)

# mass fraction
wR = rhoR/rho
wG = rhoG/rho
wB = rhoB/rho
wW = rhoW/rho
print "wR = %e, wG = %e, wB = %e, wW = %.16e" % (wR,wG,wB,wW)

print "mR/dV = %e, mG/dV = %e, mB/dV = %e" % (mR/dV,mG/dV,mB/dV)

print
print "##################"
print "# determine Dbar #"
print "##################"

# self-diffusion coefficient of water 
#  from Krynicki et al. Faraday Discuss. Chem. Soc. 66, 199 (1978)
DW = 2.30e-5
# tracer-diffusion coefficient of sucrose in water
#  from Tilley et al. J. Phys. Chem. 71, 2756 (1967)
DG = 5.25e-6
# tracer-diffusion coefficient of glucose in water
#  from Venancio et al. Biotechnology Techniques 11, 183 (1997)
DR = 7.0e-6
# tracer-diffusion coefficient of fructose in water
#  from Venancio et al. Biotechnol. Tech. 11, 183 (1997)
DB = 6.84e-6

D12 = DR*DG/DW
D13 = DR*DB/DW
D23 = DG*DB/DW 
D14 = DR
D24 = DG
D34 = DB
print "diff coeff: DR = %e, DG = %e, DB = %e, DW = %e" % (DR,DG,DB,DW)
print "D12 = %e, D13 = %e, D23 = %e, D14 = %e, D24 = %e, D34 = %e" % (D12,D13,D23,D14,D24,D34)

print
print "######################"
print "# cell Peclet number #"
print "######################"

# rough estimate for NTot
Napprox = rho/mW*dV
print "total number of molecules per cell: %e (approx)" % Napprox 

# physical parameters
kT = 1.381e-16*300. 
# kinematic viscosity from http://www.viscopedia.com/viscosity-tables/substances/water/
nu = 0.0085
print "kT = %e" % kT
print "nu = %e" % nu

# random advection
uRand = math.sqrt(kT/mW/Napprox)
print "magnitude of random advection: %e" % uRand

# cell Peclet number 
print "cell Peclet number: %e" % (uRand*dx/min(DR,DG,DB))

print
print "###############"
print "# CFL numbers #"
print "###############"

# reaction counts
dt = 5e-7
print "dt = %e" % dt

# cfl numbers
print "D*dt/dx^2 = %e" % (max(DR,DG,DB)*dt/dx**2)
print "nu*dt/dx^2 = %e" % (nu*dt/dx**2)
print "u*dt/dx = %e" % (uRand*dt/dx)

print
print "####################"
print "# structure factor #"
print "####################"

print "wR = %e, wG = %e, wB = %e, wW = %e" % (wR,wG,wB,wW)

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
