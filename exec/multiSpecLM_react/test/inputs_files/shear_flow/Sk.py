import math
import numpy as np
import numpy.linalg

################################################################################
# structure factor
################################################################################

rho0 = 1.

w1 = 0.01 
w2 = 0.02
w3 = 1-w1-w2 
print "w1 = %e, w2 = %e, w3 = %e" % (w1,w2,w3)

m1 = 1. 
m2 = 2.
m3 = 1. 
print "m1 = %e, m2 = %e, m3 = %e" % (m1,m2,m3)

ww = np.array([w1,w2,w3])            # mass fraction vector
mm = np.array([m1,m2,m3])            # molecular mass vector
mbar = 1./sum(ww/mm)                 # mean molecular mass
xx = mbar*(ww/mm)                    # mole fraction vector

W = np.diag(ww)
X = np.diag(xx)

wwww = np.outer(ww,ww)
xxxx = np.outer(xx,xx)

one = np.ones(3)
oneone = np.outer(one,one)

tmp = np.linalg.inv(X-xxxx+oneone)
Sw = mbar/rho0*np.dot(np.dot(W-wwww,tmp),W-wwww)

print rho0*rho0*Sw
