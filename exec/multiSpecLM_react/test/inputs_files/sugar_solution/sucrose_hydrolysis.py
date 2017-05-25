import math
import numpy as np
import numpy.linalg

# assume cgs units

# C12H22O11 (sucrose) + H2O (water) -> C6H12O6 (glucose) + C6H12O6 (fructose)
# R = glucose
# G = sucrose
# B = fructose
# W = water

# Avogadro number
NA = 6.022e23

# mass of a single molecule
mR = 180.16/NA 
mG = 342.30/NA 
mB = 180.16/NA
mW =  18.02/NA
print "mR = %e, mG = %e, mB = %e, mW = %e" % (mR,mG,mB,mW)

# equilibrium constant (w.r.t. concentration)
#  from Goldberg et al. J. Biol. Chem. 264, 9901 (1989)
K0 = 4.44e4

# equilibrium constant (w.r.t. number density)
K = NA/1e3*K0
print "equilibrium constant: K0 = %e (conc), K = %e (num dens)" % (K0,K)

# hydrolysis rate constant (HCl-catalyzed)
#  from Tombari et al. J. Phys. Chem. B 111, 496 (2007)
kd = 1.805e5

# dissociation rate constant (w.r.t. number density)
ka = kd/K
print "rate const: kd = %e, ka = %e" % (kd,ka)

# cell size
dx = 2.5e-6
dz = 0.16
dV = dx**2*dz
print "dx = dy = %e, dz = %e, dV = %e" % (dx,dz,dV)

# number densities of R/G/B molecules
NG = 10 
nG = NG/dV
nR = math.sqrt(K*nG)
nB = math.sqrt(K*nG)
NR = nR*dV
NB = nB*dV

print "NR = %e, NG = %e, NB = %e" % (NR,NG,NB)
print "nR = %e, nG = %e, nB = %e" % (nR,nG,nB)
print "cR = %e, cG = %e, cB = %e" % (1e3*nR/NA,1e3*nG/NA,1e3*nB/NA)

# mass densities of R/G/B molecules
rhoR = nR*mR
rhoG = nG*mG
rhoB = nB*mB
print "rhoR = %e, rhoG = %e, rhoB = %e" % (rhoR,rhoG,rhoB)

# assume rhobar = rho for all species (Boussinesq approximation)
rho = 1.
print "rhobarR = rhobarG = rhobarB = rhobarW = %e, rhoTot = %e" % (rho,rho)

# W molecules
rhoW = rho-rhoR-rhoG-rhoB
nW = rhoW/mW
NW = nW*dV
print "rhoW = %e, nW = %e, NW = %e" % (rhoW,nW,NW)

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
D14 = DR
D23 = DG*DB/DW 
D24 = DG
D34 = DB
print "diff coeff: DR = %e, DG = %e, DB = %e, DW = %e" % (DR,DG,DB,DW)
print "D12 = %e, D13 = %e, D14 = %e, D23 = %e, D24 = %e, D34 = %e" % (D12,D13,D14,D23,D24,D34)

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

# linearized reaction rate and penetration length
r = ka*(nR+nB)+kd
print "linearized reaction rate = %e" % r
print "penetration length / dx = %e" % (math.sqrt(min(DR,DG,DB)/r)/dx)

# reaction counts
dt = 2.5e-8
print "dt = %e" % dt
print "fwd reactions per time = %e" % (ka*nR*nB*dV*dt)
print "bwd reactions per time = %e" % (kd*nG*dV*dt)

# cfl numbers
print "D*dt/dx^2 = %e" % (max(DR,DG,DB,DW)*dt/dx**2)
print "nu*dt/dx^2 = %e" % (nu*dt/dx**2)
print "u*dt/dx = %e" % (uRand*dt/dx)

################################################################################
# structure factor
################################################################################

wR = rhoR/rho
wG = rhoG/rho
wB = rhoB/rho
wW = rhoW/rho

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
