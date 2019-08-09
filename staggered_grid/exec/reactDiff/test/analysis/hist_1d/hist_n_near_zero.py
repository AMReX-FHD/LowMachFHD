import sys
import math
import numpy as np

##############
# parameters #
##############

datafile = "fort.10"


n_min = -0.199           # center value of left-most bin (=smallest value)
n_max = 0.399            # center value of right-most bin (=largest value)
nbin = 300 
#n_min = -0.099            # center value of left-most bin (=smallest value)
#n_max = 0.199             # center value of right-most bin (=largest value)
#nbin = 150 

output_hist = "res.hist_near_zero"

#####################
# histogram setting #
#####################

dn = (n_max-n_min)/(nbin-1)        # width of each bin
bins = n_min+dn*np.arange(nbin)    # center values of bins 

bin_edges = list(bins-0.5*dn)      # note: len(bins) = nbin 
bin_edges.append(bins[-1]+0.5*dn)  #       len(bin_edges) = nbin + 1

#############
# read data #
#############

print "** reading %s..." % datafile

f_data = []
with open(datafile) as inf:
  for line in inf:
    nval = float(line.split()[1])          # read the second number
    f_data.append(nval)                   
f_data = np.array(f_data)

#######################
# calculate histogram #
#######################

# actually, we calculate the probability density function.
n_hist = np.histogram(f_data,bin_edges)[0].astype('float')
n_hist /= len(f_data)*dn

#############################
# output file for histogram # 
#############################

out = open(output_hist,'w')
for i in range(nbin):
  out.write("%g\t%g\n" % (bins[i],n_hist[i]))
out.close() 
print "** %s generated." % output_hist
