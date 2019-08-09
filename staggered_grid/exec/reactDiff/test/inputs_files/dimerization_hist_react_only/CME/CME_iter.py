import numpy as np
import math
import sys

#####

N0 = 40
N1bar = 20 
eps = 1e-3

if (N0%2!=0):
  print "ERROR: N0 should be even"
  sys.exit(0)

#####

def x1(N1):
  return 2*N1/float(N0+N1)

def x2(N1):
  return (N0-N1)/float(N0+N1)

#####

def avgN1(k2k1):
  k1 = 1.
  k2 = k2k1*k1
  print "k2/k1 = %e, N0 = %d, N1bar = %d" % (k2k1,N0,N1bar)

  c = np.zeros(N0/2)
  c[0] = k2*x2(0)/(k1*x1(2)**2)
  for i in range(1,N0/2):
    c[i] = (k1*x1(2*i)**2+k2*x2(2*i)-k2*x2(2*i-2)/c[i-1])/(k1*x1(2*i+2)**2)
    if (c[i]<0):
      print "** Warning: P(%d)/P(%d) is negative" % (2*i+2,2*i)

  logP = np.zeros(N0/2+1)
  logP[0] = 0.
  maxlogP = logP[0]
  signP = np.ones(N0/2+1)
  for i in range(N0/2):
    logP[i+1] = math.log(abs(c[i]))+logP[i]
    if (maxlogP<logP[i+1]):
      maxlogP = logP[i+1]
    if (c[i]<0):
      signP[i+1] = -signP[i]

  sumP = 0.
  for i in range(N0/2+1):
    sumP += signP[i]*math.exp(logP[i]-maxlogP)

  P = np.zeros(N0/2+1)
  for i in range(N0/2+1):
    P[i] = signP[i]*math.exp(logP[i]-maxlogP)/sumP

  N1vec = np.arange(0,N0+1,2)
  x1vec = 2.*N1vec/(N0+N1vec)
  avN1 = sum(P*N1vec)
  avx1 = sum(x1vec*P)
  print "<N1> = %f, <x1> = %f" % (avN1,avx1)
  print

  return avN1 
  
#####

# initial guess
x1bar = x1(N1bar)
ainit = x1bar**2/(1-x1bar)
a0 = ainit 
f0 = avgN1(a0)

# Newton's method
for trial in range(3):
  print "----- TRIAL %d" % (trial+1)

  # estimate derivative 
  a1 = (1-eps)*a0
  f1 = avgN1(a1)
  a2 = (1+eps)*a0
  f2 = avgN1(a2)
  df0 = (f2-f1)/(a2-a1)

  # correction
  da = (N1bar-f0)/df0 
  acorr = a0+da
  fcorr = avgN1(acorr)
  print "** k2/k1 = %f" % (acorr)
  print "** factor = %f\n" % (acorr/ainit)

  # for next step
  a0 = acorr
  f0 = fcorr
