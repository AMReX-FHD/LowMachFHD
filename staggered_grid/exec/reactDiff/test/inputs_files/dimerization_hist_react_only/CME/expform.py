import numpy as np
import math
import sys

#####

#N0 = 20;      N1bar = 10;     output = "res.hist_expform_20"
N0 = 40;      N1bar = 20;     output = "res.hist_expform_40"
#N0 = 60;      N1bar = 30;     output = "res.hist_expform_60"

#####

if (N0%2!=0):
  print "ERROR: N0 should be even"
  sys.exit(0)

print "N0 = %d, N1bar = %d" % (N0,N1bar)

def x1(N1):
  return 2*N1/float(N0+N1)

def x2(N1):
  return (N0-N1)/float(N0+N1)

##### expression 1

x1eq = x1(N1bar)
x2eq = 1-x1eq 

logP1 = np.zeros(N0/2+1)
maxlogP1 = 0.
for i in range(1,N0/2):
  N1 = 2*i
  N2 = (N0-N1)/2
  logP1[i] = -N1*math.log(x1(N1)/x1eq)-N2*math.log(x2(N1)/x2eq)
  if (maxlogP1<logP1[i]):
    maxlogP1 = logP1[i]

sumP1 = 0.
for i in range(1,N0/2):
  sumP1 += math.exp(logP1[i]-maxlogP1)

P1 = np.zeros(N0/2+1)
for i in range(1,N0/2):
  P1[i] = math.exp(logP1[i]-maxlogP1)/sumP1
P1[0] = 0.
P1[N0/2] = 0.

N1vec = np.arange(0,N0+1,2)
x1vec = 2.*N1vec/(N0+N1vec)
print "<N1> = %f, <x1> = %f" % (sum(P1*N1vec),sum(x1vec*P1))

##### expression 1, Gaussian approx

mean = N1bar
var = N1bar*(1-(N1bar/float(N0))**2)

P1G = np.zeros(N0/2+1)
for i in range(N0/2+1):
  N1 = 2*i
  P1G[i] = math.exp(-(N1-mean)**2/2/var)/math.sqrt(2*math.pi*var)

P1G *= 2

#### file output

np.savetxt(output,np.c_[N1vec,P1,P1G])
print "%s generated" % output


