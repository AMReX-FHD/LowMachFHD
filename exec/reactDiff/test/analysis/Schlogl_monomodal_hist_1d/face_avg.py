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
dv = 5.
avg_type = 1   # 1=arith (C0-smoothed Heaviside)
               #   10=discontinous Heaviside
               #   11=C1-smoothed Heaviside
               #   12=C2-smoothed Heaviside
               # 2=geom
               # 3=harmo

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
print "** ncell=%d, dv=%g, avg_type=%d" %(ncell,dv,avg_type)

####################
# datafile=fort.10 #
####################

lcnt = 0        # line count
tval_prev = 0.  # t value at the previous line
eps = 1.e-10    # small number for checking whether tval==tval_prev or not

sum1 = 0.       # sum of values of n (over all lines)
sum2 = 0.       # sum of values of n*n (over all lines)
sum3 = 0.       # sum of values of n*n_neigh (n_neigh = nleft and nright) (over all lines)
sum41 = 0.      # sum of avg_type=1  face average values of (n,n_neigh) (over all lines)
sum410 = 0.     # sum of avg_type=10 face average values of (n,n_neigh) (over all lines)
sum411 = 0.     # sum of avg_type=11 face average values of (n,n_neigh) (over all lines)
sum412 = 0.     # sum of avg_type=12 face average values of (n,n_neigh) (over all lines)
sum42 = 0.      # sum of avg_type=2  face average values of (n,n_neigh) (over all lines)
sum43 = 0.      # sum of avg_type=3  face average values of (n,n_neigh) (over all lines)

nvals = np.zeros(ncell)

def face_avg_1(n1,n2):
  if (n1<0. or n2<0.):
    return 0.
  else:
    tmp1 = min(dv*n1,1.)
    tmp2 = min(dv*n2,1.)
    return (n1+n2)/2.*tmp1*tmp2

def face_avg_10(n1,n2):
  if (n1<0. or n2<0.):
    return 0.
  else:
    return (n1+n2)/2.

def face_avg_11(n1,n2):
  if (n1<0. or n2<0.):
    return 0.
  else:
    tmp1 = dv*n1
    if (tmp1<1.):
      tmp1 = (3.-2.*tmp1)*tmp1**2
    else:
      tmp1 = 1.
    tmp2 = dv*n2
    if (tmp2<1.):
      tmp2 = (3.-2.*tmp2)*tmp2**2
    else:
      tmp2 = 1.
    return (n1+n2)/2.*tmp1*tmp2
      
def face_avg_12(n1,n2):
  if (n1<0. or n2<0.):
    return 0.
  else:
    tmp1 = dv*n1
    if (tmp1<1.):
      tmp1 = (10-15*tmp1+6*tmp1**2)*tmp1**3
    else:
      tmp1 = 1.
    tmp2 = dv*n2
    if (tmp2<1.):
      tmp2 = (10-15*tmp2+6*tmp2**2)*tmp2**3
    else:
      tmp2 = 1.
    return (n1+n2)/2.*tmp1*tmp2  

def face_avg_2(n1,n2):
  if (n1<0. or n2<0.):
    return 0.
  else:
    return math.sqrt(n1*n2)

def face_avg_3(n1,n2):
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
        sum41  += face_avg_1(nn,nleft)  + face_avg_1(nn,nright)
        sum410 += face_avg_10(nn,nleft) + face_avg_10(nn,nright)
        sum411 += face_avg_11(nn,nleft) + face_avg_11(nn,nright)
        sum412 += face_avg_12(nn,nleft) + face_avg_12(nn,nright)
        sum42  += face_avg_2(nn,nleft)  + face_avg_2(nn,nright)
        sum43  += face_avg_3(nn,nleft)  + face_avg_3(nn,nright)

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

# mean of (face_avg(n,nleft)+face_avg(n,nright))/2
sum41  /= 2*lcnt  
sum410 /= 2*lcnt  
sum411 /= 2*lcnt  
sum412 /= 2*lcnt  
sum42  /= 2*lcnt  
sum43  /= 2*lcnt  

print "face_avg_arith1=  %.10f" % sum41 
print "face_avg_arith10= %.10f" % sum410 
print "face_avg_arith11= %.10f" % sum411 
print "face_avg_arith12= %.10f" % sum412 
print "face_avg_geom=    %.10f" % sum42 
print "face_avg_harmo=   %.10f" % sum43 

if (avg_type==1):
  face_avg = sum41
elif (avg_type==10):
  face_avg = sum410
elif (avg_type==11):
  face_avg = sum411
elif (avg_type==12):
  face_avg = sum412
elif (avg_type==2):
  face_avg = sum42
elif (avg_type==3):
  face_avg = sum43

print "face_avg= %.10f" % face_avg
print "-1/(Nc-1)= %.10f" % (-1./(ncell-1))
print "(Nc-1)/Nc*<face_avg>/dv= %.10f" % (float(ncell-1)/ncell*face_avg/dv)
print "-1/Nc*<face_avg>/dv= %.10f" % (-1./ncell*face_avg/dv)
