import numpy as np
import math
import sys

#####

#N0 = 20;      N1bar = 10;  factor = 0.967923; output = "res.hist_CME_20"
N0 = 40;      N1bar = 20;  factor = 0.983646; output = "res.hist_CME_40"
#N0 = 50;      N1bar = 25;  factor = 0.986867; output = "res.hist_CME_50"
#N0 = 60;      N1bar = 30;  factor = 0.989028; output = "res.hist_CME_60"
#N0 = 80;      N1bar = 40;  factor = 0.991745; output = "res.hist_CME_80"
#N0 = 100;     N1bar = 50;  factor = 0.993383; output = "res.hist_CME_100"

#N0 = 1000;    N1bar = 5e2; factor = 0.999334; output = "res.hist_CME_1e3"
#N0 = 10000;   N1bar = 5e3; factor = 0.999933; output = "res.hist_CME_1e4"
#N0 = 100000;  N1bar = 5e4; factor = 0.999993; output = "res.hist_CME_1e5"
#N0 = 1000000; N1bar = 5e5; factor = 0.999999; output = "res.hist_CME_1e6"

#####

if (N0%2!=0):
  print "ERROR: N0 should be even"
  sys.exit(0)

print "N0 = %d, N1bar = %d" % (N0,N1bar)

def x1(N1):
  return 2*N1/float(N0+N1)

def x2(N1):
  return (N0-N1)/float(N0+N1)

k1 = 1.
k2 = factor*k1*x1(N1bar)**2/x2(N1bar)
print "k2/k1 = %f, factor = %f" % (k2/k1,factor)

#####

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
print "<N1> = %f, <x1> = %f" % (sum(P*N1vec),sum(x1vec*P))

#### file output

np.savetxt(output,np.c_[N1vec,P])
print "%s generated" % output
