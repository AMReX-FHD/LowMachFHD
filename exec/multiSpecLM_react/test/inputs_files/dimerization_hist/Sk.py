import numpy as np

rho0 = 40.

w1 = 0.5
w2 = 0.5
wvec = np.array([w1,w2])

m1 = 1.
m2 = 2.
mvec = np.array([m1,m2])

mbar = 1./sum(wvec/mvec)
xvec = mbar*(wvec/mvec)

W = np.diag(wvec)
X = np.diag(xvec)

ww = np.outer(wvec,wvec)
xx = np.outer(xvec,xvec)

onevec = np.ones(2)
oneone = np.outer(onevec,onevec)

tmp = np.linalg.inv(X-xx+oneone)
Sw = mbar/rho0*np.dot(np.dot(W-ww,tmp),W-ww)

print rho0*rho0*Sw
