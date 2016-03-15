import sys
import numpy as np
import math

# get RUNNAME and NRUN from command line
if (len(sys.argv)==3):
  RUNNAME = sys.argv[1] 
  NRUN = int(sys.argv[2])
else:
  print "usage: python cond_exp_stat.py RUNNAME NRUN"
  sys.exit(0) 

# read the first datafile to know the numbers of columns and lines
datafile = "%s%d/res.cond_exp" % (RUNNAME,1)
with open(datafile) as inf:
  linecnt = 0
  for line in inf:
    linecnt += 1
nline = linecnt
ncol = len(line.split())
print "** %d lines with %d columns **" % (nline,ncol)

# main part
nn = np.zeros(nline)
nzero = np.zeros(nline)
ncnt = np.zeros(nline)
nsum = np.zeros((nline,ncol-2))
n1nj1 = np.zeros((nline,ncol-2))
n1nj2 = np.zeros((nline,ncol-2))

for i in range(NRUN):
  datafile = "%s%d/res.cond_exp" % (RUNNAME,i+1)
  print datafile
 
  with open(datafile) as inf:
    linecnt = 0
    for line in inf:
      parts = line.split()
      vals = np.array(map(float,parts))

      if (i==0):
        nn[linecnt] = vals[0]

      if (vals[1]==0.):
        nzero[linecnt] += 1
      else:
        ncnt[linecnt] += vals[1]
        nsum[linecnt] += vals[2:ncol]
        tmp = vals[2:ncol]/vals[1] 
        n1nj1[linecnt] += tmp 
        n1nj2[linecnt] += tmp*tmp 

      linecnt += 1

# output
output = "%s/res.cond_exp_stat" % RUNNAME
out = open(output,"w")
for k in range(nline):
  if (nzero[k]==0):
    out.write("%f" % nn[k])
    for i in range(ncol-2):
      tmp1 = n1nj1[k,i]/NRUN
      tmp2 = n1nj2[k,i]/NRUN
      std = math.sqrt((tmp2-tmp1*tmp1)/(NRUN-1))
      out.write("\t%g\t%g\t%g" % (nsum[k,i]/ncnt[k],tmp1,std))
    out.write("\n")
  elif (nzero[k]==NRUN):
    out.write("#\t%f\t%d\n" % (nn[k],nzero[k]))
  else:
    out.write("#\t%f\t%d\t%d" % (nn[k],nzero[k],ncnt[k]))
    for i in range(ncol-2):
      out.write("\t%g" % (nsum[k,i]/ncnt[k]))
    out.write("\n")
out.close()

print "res.cond_exp_stat generated"
