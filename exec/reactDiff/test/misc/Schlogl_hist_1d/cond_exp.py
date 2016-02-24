import sys
import numpy as np
import math

# this script is for the calculation of conditional expectation < n_{1+j} | n_1 >.
# note that this is a function of j as well as the value of n_1.
# cell_dist contains the list of j to be calculated.
# (nbin_val, nbin_cnt) is the histogram of n_1.
# in the output file, there are 2+len(cell_dist) columns.
# 1st col: nbin_val
# 2nd col: nbin_cnt
# j+3th col: sum of n_{1+cell_dist[j]}
# Hence, by calculating (j+3th col)/(2nd col), < n_{1+j} | n_1 > can be estimated.

# run this script on the directory containing fort.10
# (e.g. cd TEST1; python ../cond_exp.py)

# input parameters

datafile="fort.10"
ncell = 8

nbin_min = -0.5	 
nbin_max = 2.5       
nbin_num = 30 

cell_dist = [1,-1,2,-2,3,-3,4]

output = "res.cond_exp"

# auxiliary variables

nbin_dn = float(nbin_max-nbin_min)/nbin_num
nbin_val = nbin_min+nbin_dn*(np.arange(nbin_num)+0.5)
nbin_cnt = np.zeros(nbin_num)

one_frame_data = np.zeros(ncell)

A = np.zeros((nbin_num,len(cell_dist)))

# read file
lcnt = 0        # line count
tval_prev = 0.  # t value at the previous line
eps = 1.e-10    # small number for checking whether tval==tval_prev or not

with open(datafile) as inf:
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

    # save n value
    one_frame_data[(lcnt-1)%ncell] = nval
  
    if (lcnt%ncell==0):                          # at the end of the block
      # main calculation 
      for i in range(ncell):
        nval = one_frame_data[i]
        if (nval<nbin_min or nval>=nbin_max):
          print "nval=%f is out of [%f,%f]" % (nval,nbin_min,nbin_max)
        else:
          pos = int(math.floor((nval-nbin_min)/nbin_dn))
          nbin_cnt[pos] += 1
          for j in range(len(cell_dist)):
            A[pos,j] += one_frame_data[(i+cell_dist[j])%ncell]

    # before reading the next line
    tval_prev = tval 

# final check
nblock = lcnt/ncell
if (nblock*ncell==lcnt):
  print "** %s: ncell=%d, nblock=%d, lcnt=%d" % (datafile,ncell,nblock,lcnt)
else:
  print "Error: nblock*ncell != lcnt"
  sys.exit()

# output
out = open(output,"w")

for i in range(nbin_num):
  out.write("%f\t%e\t" % (nbin_val[i],nbin_cnt[i]))
  for j in range(len(cell_dist)):
    out.write("%e\t" % A[i,j])
  out.write("\n")

out.close()

print "%s generated" % output
