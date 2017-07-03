import numpy as np

datafile = "fort.10"
nspecies = 2
massratio = [1.,2.]
spec_analysis = 0
output = "res.hist"

# dz=1
xmin = 0.66665
xmax = 0.66668
# dz=1e-1
#xmin = 0.66662
#xmax = 0.66672
# dz=1e-2
#xmin = 0.6665
#xmax = 0.6668
# dz=1e-3
#xmin = 0.666
#xmax = 0.6672
# dz = 1e-4
#xmin = 0.665
#xmax = 0.668

Nbin = 51

# read data
wvec = np.loadtxt(datafile,usecols=range(1,nspecies+1),unpack=True)

# calculate mole fraction
nsumvec = np.zeros(len(wvec[0]))
for spec in range(nspecies):
   if (spec==spec_analysis):
      nsavvec = wvec[spec]/massratio[spec]
      nsumvec += nsavvec
   else:
      nsumvec += wvec[spec]/massratio[spec]
xvec = nsavvec/nsumvec

# compute histogram

if (xmin>min(xvec) or xmax<max(xvec)):
  print "x_min = %e, x_max = %e" % (min(xvec),max(xvec))

dx = (xmax-xmin)/(Nbin-1)
[hist,bin_edges] = np.histogram(xvec,Nbin,[xmin-dx/2,xmax+dx/2],density=True)
bins = bin_edges[0:Nbin]+dx/2

# write histogram
np.savetxt(output,np.c_[bins,hist])
