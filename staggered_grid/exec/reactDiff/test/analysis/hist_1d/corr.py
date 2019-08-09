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

# if additional arguments are given
if (len(sys.argv)!=1):   
  if (len(sys.argv)==2):       # for one additional arguments
    ncell = int(sys.argv[1])
  else:                        # otherwise, generate error
    print "Error: only one additional argument can be given."
    sys.exit()

# print parameter values
print "** ncell=%d" %(ncell)

####################
# datafile=fort.10 #
####################

lcnt = 0        # line count
tval_prev = 0.  # t value at the previous line
eps = 1.e-10    # small number for checking whether tval==tval_prev or not

sum1 = 0.       # sum of values of n (over all lines)
sum2 = 0.       # sum of values of n*n (over all lines)
sum3 = 0.       # sum of values of n*n_neigh (n_neigh = nleft and nright) (over all lines)

nvals = np.zeros(ncell)

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

var = sum2-sum1*sum1
cov = sum3-sum1*sum1
corr = cov/var

print "<n>= %.10f" % (sum1)
print "<n^2>= %.10f" % (sum2)
print "Var[n]= %.10f" % (var)

print "<n*n_neigh>= %.10f" % (sum3)
print "Cov[n,n_neigh]= %.10f" % (cov)
print "Corr[n,n_neigh]= %.10f" % (corr)
