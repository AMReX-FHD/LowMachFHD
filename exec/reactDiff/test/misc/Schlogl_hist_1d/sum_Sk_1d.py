import sys

# this script analyzes the HydroGrid output for 1d system.
# it calculates:
# - first and second moments of n and nbar(=spatial average of n)
# - S(0) and sum of S(k)
# note that these two kinds of things can be related through Parseval's theorem.

##############
# parameters #
##############

# datafile1 contains number densities of the cells at different times.
# datafile2 contains the structure factor S(k) for nonzero k. 
datafile1="fort.10"
datafile2="Hist1D.S_k.pair=1.Re.dat"

# inputs_hist_1d assumes ncell*1 cells in 2d.
# each cell has dx*dx in dimensions with cross_section=dz.
ncell = 64
dx = 1.
dz = 10.

# if additional arguments are given
if (len(sys.argv)!=1):   
  if (len(sys.argv)==3):       # for two additional arguments
    datafile1 = sys.argv[1]  
    datafile2 = sys.argv[2] 
  elif (len(sys.argv)==6):     # for five additional arguments
    datafile1 = sys.argv[1]
    datafile2 = sys.argv[2]
    ncell = int(sys.argv[3])
    dx = float(sys.argv[4])
    dz = float(sys.argv[5])
  else:                        # otherwise, generate error
    print "Error: only two or five additional arguments can be given."
    sys.exit()

# print parameter values
print "** ncell=%d, dx=%g, dz=%g" %(ncell,dx,dz)

#####################
# datafile1=fort.10 #
#####################

lcnt = 0        # line count
tval_prev = 0.  # t value at the previous line
eps = 1.e-10    # small number for checking whether tval==tval_prev or not

sum1 = 0.       # sum of values of n (over all lines)
sum2 = 0.       # sum of values of n*n (over all lines)
sum3 = 0.       # sum of values of n (in a block)
sum4 = 0.       # sum of values of nbar*nbar (where nbar=sum3/ncell, spatial average of n)

with open(datafile1) as inf:
  for line in inf:

    # read t and n values
    lcnt += 1
    parts = line.split()
    tval = float(parts[0])
    nval = float(parts[1])

    # confirm that the values are consistent with ncell
    if ((lcnt-1)%ncell==0):                      # at the beginning of new block
      if (lcnt!=1 and abs(tval-tval_prev)<eps):  # if tval==tval_prev,
        print "Error: check ncell"               # generate error
        sys.exit()
    else:                                        # at other lines in the block
      if (abs(tval-tval_prev)>eps):              # if tval!=tval_prev,
        print "Error: check ncell"               # generate error
        sys.exit() 
  
    # calculate sum1, sum2, sum3
    sum1 += nval
    sum2 += nval*nval
    sum3 += nval

    if (lcnt%ncell==0):                          # at the end of the block
      nbar = sum3/ncell     
      sum4 += nbar*nbar
      sum3 = 0.

    # before reading the next line
    tval_prev = tval 

nblock = lcnt/ncell
if (nblock*ncell==lcnt):
  print "** %s: ncell=%d, nblock=%d, lcnt=%d" % (datafile1,ncell,nblock,lcnt)
else:
  print "Error: nblock*ncell != lcnt"
  sys.exit()

sum1 /= lcnt    # first moment of n
sum2 /= lcnt    # second moment of n
sum4 /= nblock  # second moment of nbar

print "<n>= %.5f" % (sum1)
print "<n^2>= %.5f" % (sum2)
print "Var[n]= %.5f" % (sum2-sum1*sum1)
print "<nbar*nbar>=%.5f" % (sum4)

second_moment_n = sum2 
second_moment_nbar = sum4

########################################
# datafile2="Hist1D.S_k.pair=1.Re.dat" #
########################################

lcnt = 0   # line count
sum1 = 0.  # sum of S(k) (note: S(0) is missing in the file)

with open(datafile2) as inf:
  for line in inf:

    lcnt += 1
    if (lcnt<=2):                   # skip the first two lines
      continue

    Skval = float(line.split()[1])  # read the second number as S(k)
    sum1 += Skval

nSkval = lcnt-2                     # note: there are two comment lines
if (nSkval==ncell-1):
  print "** %s: %d values of S(k) (avg=%g)" % (datafile2,nSkval,sum1/lcnt)
else:
  print "Error: ncell-1 != nSkval"
  sys.exit()

V = ncell*dx*dx*dz
S0 = V*second_moment_nbar

print "sum of S(k) (k!=0) / V = %g  (cf. <n^2> - <nbar^2> = %g)" % (sum1/V,second_moment_n-second_moment_nbar)
print "sum of S(k) / V = %g (cf. <n^2> = %g)" % ((S0+sum1)/V,second_moment_n)
