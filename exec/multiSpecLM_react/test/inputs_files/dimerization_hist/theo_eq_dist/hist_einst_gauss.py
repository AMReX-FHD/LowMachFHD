import numpy as np
import math

N0 = 40
kB = 1.
dx = 0.1

N2eq = N0/4.
Neq = (N0-2.*N2eq)+N2eq
x2eq = N2eq/Neq

N2 = dx*np.arange(-1/dx,(N0/2+1)/dx+1)
N = (N0-2.*N2)+N2
x2 = N2/N

dS2 = -kB*N0**2/N2eq/(N0-N2eq)/(N0-2*N2eq)

PGauss = np.exp(0.5*dS2/kB*(N2-N2eq)**2)
PGauss = PGauss/sum(PGauss)/dx

np.savetxt("res.hist_einst_gauss",np.transpose([N2,PGauss]))
