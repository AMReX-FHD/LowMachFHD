import sys
import numpy as np
import math

##############
# parameters #
##############

# datafile contains number densities of the cells at different times.
datafile="fort.10"

# inputs_hist_1d assumes ncell*1 cells in 2d.
# each cell has dx*dx in dimensions with cross_section=dz.
ncell = 64
dv = 10.
avg_type = 1                   # 1=arith 2=geom 3=harmo

# if additional arguments are given
if (len(sys.argv)!=1):   
  if (len(sys.argv)==4):     # for three additional arguments
    ncell = int(sys.argv[1])
    dv = float(sys.argv[2])
    avg_type = int(sys.argv[3])
  else:                        # otherwise, generate error
    print "Error: only one or four additional arguments can be given."

    sys.exit()

# print parameter values
print "** ncell=%d, dv=%g" %(ncell,dv)

####################
# datafile=fort.10 #
####################

lcnt = 0        # line count
tval_prev = 0.  # t value at the previous line
eps = 1.e-10    # small number for checking whether tval==tval_prev or not

sum1 = 0.       # sum of values of n (over all lines)
sum2 = 0.       # sum of values of n*n (over all lines)
sum3 = 0.       # sum of values of n*n_neigh (n_neigh = nleft and nright) (over all lines)
sum41 = 0.      # sum of arithmetic-mean face average values of (n,n_neigh) (over all lines)
sum42 = 0.      # sum of geometric-nean face average values of (n,n_neigh) (over all lines)
sum43 = 0.      # sum of harmonic-mean face average values of (n,n_neigh) (over all lines)

nvals = np.zeros(ncell)

def face_avg_arith(n1,n2):
  if (n1<0. or n2<0):
    return 0.
  else:
    tmp = dv*min(n1,n2)
    if (tmp<1.5):
      tmp = 1./(1.+math.exp(-12.*(tmp-0.5)))    # smoothed Headviside function
      return tmp*(n1+n2)/2.
    else:
      return (n1+n2)/2.

def face_avg_geom(n1,n2):
  if (n1<0. or n2<0.):
    return 0.
  else:
    return math.sqrt(n1*n2)

def face_avg_harmo(n1,n2):
  if (n1<1.e-10 or n2<1.e-10):
    return 0.
  else:
    return 2./(1./n1+1./n2) 

with open(datafile) as inf:
  for line in inf:

    # read t and n values
    lcnt += 1
    parts = line.split()
    tval = float(parts[0])
    nval = float(parts[1])

    nvals[(lcnt-1)%ncell] = nval

    # confirm that the values are consistent with ncell
    if ((lcnt-1)%ncell==0):                      # at the beginning of new block
      if (lcnt!=1 and abs(tval-tval_prev)<eps):  # if tval==tval_prev,
        print "Error: check ncell"               # generate error
        sys.exit()
    else:                                        # at other lines in the block
      if (abs(tval-tval_prev)>eps):              # if tval!=tval_prev,
        print "Error: check ncell"               # generate error
        sys.exit() 

    # calculate quantities at the end of the block
    if (lcnt%ncell==0):
      for i in range(ncell):
        nn = nvals[i]
        nleft = nvals[i-1]
        nright = nvals[(i+1)%ncell]
        sum1 += nn 
        sum2 += nn*nn 
        sum3 += nn*(nleft+nright)
        sum41 += face_avg_arith(nn,nleft)+face_avg_arith(nn,nright)
        sum42 += face_avg_geom(nn,nleft)+face_avg_geom(nn,nright)
        sum43 += face_avg_harmo(nn,nleft)+face_avg_harmo(nn,nright)

    # before reading the next line
    tval_prev = tval 

nblock = lcnt/ncell
if (nblock*ncell==lcnt):
  print "** %s: ncell=%d, nblock=%d, lcnt=%d" % (datafile,ncell,nblock,lcnt)
else:
  print "Error: nblock*ncell != lcnt"
  sys.exit()

sum1 /= lcnt     # first moment of n
sum2 /= lcnt     # second moment of n
sum3 /= 2*lcnt   # mean of (n*nleft+n*nright)/2

face_avg = np.zeros(3)
face_avg[0] = sum41/(2*lcnt)  # mean of (face_avg_arith(n,nleft)+face_avg_arith(n,nright))/2
face_avg[1] = sum42/(2*lcnt)  # mean of (face_avg_geom(n,nleft)+face_avg_geom(n,nright))/2
face_avg[2] = sum43/(2*lcnt)  # mean of (face_avg_harmo(n,nleft)+face_avg_harmo(n,nright))/2

var = sum2-sum1*sum1
cov = sum3-sum1*sum1
corr = cov/var

print "<n>= %.10f" % (sum1)
print "<n^2>= %.10f" % (sum2)
print "Var[n]= %.10f" % (var)

print "<n*n_neigh>= %.10f" % (sum3)
print "Cov[n,n_neigh]= %.10f" % (cov)
print "Corr[n,n_neigh]= %.10f" % (corr)

print "face_avg_arith= %.10f" % (face_avg[0]) 
print "face_avg_geom= %.10f" % (face_avg[1]) 
print "face_avg_harmo= %.10f" % (face_avg[2]) 

print "-1/(Nc-1)= %.10f" % (-1./(ncell-1))
print "(Nc-1)/Nc*<face_avg>/dv= %.10f" % (float(ncell-1)/ncell*face_avg[avg_type-1]/dv)
print "-1/Nc*<face_avg>/dv= %.10f" % (-1./ncell*face_avg[avg_type-1]/dv)
