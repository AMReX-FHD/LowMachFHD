import sys
import math
import numpy as np
import matplotlib.pyplot as plt

# this script computes the histogram for the number of molecules N of a chosen 
#  species in a cell near zero.
# screen output as well as file output are generated.
#
# usage:
#  python hist_N_near_zero.py
#  python hist_N_near_zero.py input_column norm_factor nbin Nmin Nmax
#    show_plot(yes/no) save_plot(yes/no)

##############
# parameters #
##############

# datafile contains (number or mass) densities of each cell at different times.
datafile = "fort.10"

# output filename for histogram
# if "none" is given, histogram is not saved.
output_hist = "res.hist_near_zero"
# output filename for figure
output_fig_hist = "hist_near_zero.png"

# input_column indicates which column of densities is to be analyzed.
# note: c convention is used. e.g. first column = 0, second column = 1, ...
input_column = 1

# nor_factor is required to convert densities into the number of molecules.
# hence, this quantity corresponds to the value of density representing a single molecule.
norm_factor = 3.e-9

# nbin is the number of bins
nbin = 101 

# Nmin/Nmax: min/max values of histogram bins
Nmin = -2.
Nmax = 4.

# flag for showing the plot on the screen
show_plot = True
# flag for saving the figure on a file
save_plot = True

# if additional arguments are given
if (len(sys.argv)!=1):
  if (len(sys.argv)==8):                # for seven additional arguments
    input_column = int(sys.argv[1])
    norm_factor = float(sys.argv[2])
    nbin = int(sys.argv[3])
    Nmin = float(sys.argv[4])
    Nmax = float(sys.argv[5])
    if (sys.argv[6]=="no"):
      show_plot = False
    if (sys.argv[7]=="no"):
      save_plot = False
  else:                                 # otherwise, abort
    print "Error: only seven additional arguments can be given."
    print "-> input_column norm_factor nbin Nmin Nmax show_plot(yes/no) save_plot(yes/no)"
    sys.exit(0)

#####################
# histogram setting #
#####################

dN = (Nmax-Nmin)/(nbin-1)               # width of each bin
bins = Nmin+dN*np.arange(nbin)          # center values of bins 

bin_edges = list(bins-0.5*dN)           # note: len(bins) = nbin 
bin_edges.append(bins[-1]+0.5*dN)       #       len(bin_edges) = nbin + 1

#############
# read data #
#############

print "** reading %s..." % datafile

f_data = []
with open(datafile) as inf:
  for line in inf:
    # read the number in input_column and then convert it into the number of molecules
    val = float(line.split()[input_column])/norm_factor
    f_data.append(val)                   
f_data = np.array(f_data)

#######################
# calculate histogram #
#######################

# actually, we calculate the probability density function.
Nhist = np.histogram(f_data,bin_edges)[0].astype('float')
Nhist /= len(f_data)*dN


############################
# show plots on the screen #
############################

fig, a = plt.subplots()

a.plot(bins,Nhist,'-r')

plt.xlabel("number of molecules N")
plt.ylabel("probability density")

if (show_plot):
  plt.show()
if (save_plot):
  fig.savefig(output_fig_hist)
  print "** %s generated." % output_fig_hist

#############################
# output file for histogram # 
#############################

out = open(output_hist,'w')
for i in range(nbin):
  out.write("%g\t%g\n" % (bins[i],Nhist[i]))
out.close() 
print "** %s generated." % output_hist
