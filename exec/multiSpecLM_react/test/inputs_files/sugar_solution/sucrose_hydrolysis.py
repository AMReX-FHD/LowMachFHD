import math
import numpy as np

################################################################################
# hydrolysis of sucrose:
# sucrose + water <-> glucose + fructose
#    G    +   W   <->    R    +    B
# C12H22O11 + H2O <-> C6C12O6 + C6O12O6
################################################################################

print
print "#########################"
print "# Hydrolysis of sucrose #"
print "#########################"
print
print "* in cgs units"
print 

################################################################################
# physical parameters
################################################################################

# basic constants
NA = 6.022e23           # Avogadro's number
kB = 1.381e-16          # Bolzmann's constant

# parameters for aqueous solution
T = 293.                # temperature
rho0 = 1.               # mass density
eta = 1.e-2             # shear viscosity

# molecular mass of each species (per single molecule)
mR = 180.16/NA
mG = 342.30/NA
mB = 180.16/NA
mW =  18.02/NA

# equilibrium constant
#  from Goldberg et al. J. Biol. Chem. 264, 9901 (1989)
Kc = 4.44e4             # for concentration (mol/L), c_i = n_i/NA*1e3
K = NA/1e3*Kc           # for number density

# rate constant (for dissociation)
#  from Tombari et al. J. Phys. Chem. B 111, 496 (2007)
#kd = 1.805e-5
kd = 10. 

# trace diffusion coefficient of glucose in water
#  from Venancio et al. Biotechnology Techniques 11, 183 (1997)
DR = 7.0e-6

# trace diffusion coefficient of sucrose in water
#  from Tilley et al. J. Phys. Chem. 71, 2756 (1967)
DG = 5.25e-6

# trace diffusion coefficient of fructose in water
#  from Venancio et al. Biotechnol. Tech. 11, 183 (1997)
DB = 6.84e-6

# self-diffusion coefficient of water 
#  from Krynicki et al. Faraday Discuss. Chem. Soc. 66, 199 (1978)
DW = 2.30e-5

################################################################################
# numerical parameters
################################################################################

# number of sucrose molecules
NG = 10

# cell dimensions
dx = 1.e-4
dz = 1.e-4

# time step size
dt = 5.e-5

################################################################################
# equilibrium composition
################################################################################

dv = dx*dx*dz

nG = NG/dv
nR = math.sqrt(K*nG) 
nB = nR

rhoR = mR*nR
rhoG = mG*nG
rhoB = mB*nB
rhoW = rho0 - rhoR - rhoG - rhoB

nW = rhoW/mW

NR = nR*dv 
NB = nB*dv 
NW = nW*dv 

cR = 1e3/NA*nR
cG = 1e3/NA*nG
cB = 1e3/NA*nB
cW = 1e3/NA*nW

print "dx= %e\tdz= %e\tdv= %e"  % (dx,dz,dv)
print 

print "mR=   %e\tmG=   %e\tmB=   %e\tmW=   %e" % (mR,mG,mB,mW)
print "nR=   %e\tnG=   %e\tnB=   %e\tnW=   %e" % (nR,nG,nB,nW)
print "rhoR= %e\trhoG= %e\trhoB= %e\trhoW= %e" % (rhoR,rhoG,rhoB,rhoW)
print "NR=   %e\tNG=   %e\tNB=   %e\tNW=   %e" % (NR,NG,NB,NW)
print "cR=   %e\tcG=   %e\tcB=   %e\tcW=   %e" % (cR,cG,cB,cW)
print

#print "check:"
#print "cR*cB/cG %e" % (cR*cB/cG)
#print "Kc       %e" % (Kc)
#print 

if (rhoW/rho0<0.99):
  print "Warning: rhoW/rho0= %e" % (rhoW/rho)

################################################################################
# reaction
################################################################################

ka = kd/K
r = kd+ka*(nR+nB)

print "K=  %e\tKc= %e" % (K,Kc)
print "kd= %e\tka= %e" % (kd,ka)
print "r= %e" % (r)
print

################################################################################
# diffusion
################################################################################

D12 = DR*DG/DW
D13 = DR*DB/DW
D14 = DR
D23 = DG*DB/DW
D24 = DG
D34 = DB

pendep = math.sqrt(min(DR,DG,DB)/r)

print "D12= %e\tD13= %e\tD23= %e\tD14= %e\tD24= %e\tD34= %e" % (D12,D13,D23,D14,D24,D34)
print "pendep= %e\tdx= %e\tpendep/dx= %e" % (pendep,dx,pendep/dx)
print

################################################################################
# advection
################################################################################

urand = math.sqrt(kB*T/rho0/dv)
cellPe = urand*dx/min(DR,DG,DB)

print "urand= %e" % (urand)
print "cell Pe= %e" % (cellPe)
print

################################################################################
# CFL numbers 
################################################################################

advCFL = urand*dt/dx
massCFL = max(DR,DG,DB)*dt/dx**2
momCFL = (eta/rho0)*dt/dx**2

print "dt= %e" % (dt)
print "adv_CFL= %e" % (advCFL)
print "mass_diff_CFL= %e" % (massCFL)
print "mom_diff_CFL= %e" % (momCFL)
print

################################################################################
# reaction occurrences 
################################################################################

print "forward reactions/cell/timestep= %e" % (kd*nG*dv*dt)
print "backard raactions/cell/timestep= %e" % (ka*nR*nB*dv*dt)
print

################################################################################
# structure factor
################################################################################

wR = rhoR/rho0
wG = rhoG/rho0
wB = rhoB/rho0
wW = rhoW/rho0

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
Sw = mbar/rho0*np.dot(np.dot(W-wwww,tmp),W-wwww)

print rho0*rho0*Sw

