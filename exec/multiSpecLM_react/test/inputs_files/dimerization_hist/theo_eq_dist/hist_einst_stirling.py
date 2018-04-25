import numpy as np
import math

N0 = 40
kB = 1.

N2eq = N0/4.
Neq = (N0-2.*N2eq)+N2eq
x2eq = N2eq/Neq

N2 = np.arange(N0/2+1)
N2 = N2[1:-1]                   # omit N1=0 and N2=0
N = (N0-2.*N2)+N2
x2 = N2/N

dS = -kB*N*((1-x2)*np.log((1-x2)/(1-x2eq))+x2*np.log(x2/x2eq))

P = np.exp(dS/kB)
P = P/sum(P)

np.savetxt("res.hist_einst_stirling",np.transpose([N2,P]))
