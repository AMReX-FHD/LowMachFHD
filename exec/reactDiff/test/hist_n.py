import sys
import math
import numpy as np
import scipy.stats 
import matplotlib.pyplot as plt

# this script calculates the histogram of n.
# it is compared with Gaussian and Poisson distributions.
# screen output as well as file output are generated.

##############
# parameters #
##############

# datafile contains number densities n of the cells at different times.
datafile = "fort.10"

# average value of n (n_av) and the volume of a cell (dV) are needed
#  to calculate Gaussian and Poisson distributions. 
# for 1d and 2d, dV=dx*dy*cross_section.
n_av = 1. 
dV = 5.    

# the range of histogram [n_min,n_max] is determined by [N_min/dV,N_max/dV].
# number of bins (nbin) is determined by nbin_in:
#  if nbin_in=0, nbin = N_max-N_min+1. 
#  otherwise, nbin = nbin_in.
N_min = -5
N_max = 15
nbin_in = 0

# flags for showing semi-log-scale or linear plot on the screen
plot_semilogy = True
plot_linear = True

# output filenames for (n_hist,Gauss,Stirling), (Gaussian,Stirling), or Poisson
# if "" is given, that file is omitted.  
output_hist = "res.hist" 
output_cont = "res.hist_cont" 
output_poiss = "res.hist_poiss" 

# if additional arguments are given
if (len(sys.argv)!=1):
  if (len(sys.argv)==4):    # for three additional arguments
    output_hist = sys.argv[1]
    output_cont = sys.argv[2]
    output_poiss = sys.argv[3]
  elif (len(sys.argv)==6):  # for five additonal arguments
    output_hist = sys.argv[1]
    output_cont = sys.argv[2]
    output_poiss = sys.argv[3]
    n_av = float(sys.argv[4])
    dV = float(sys.argv[5])
  else:                     # otherwise, generate error
    print "Error: only three or five additonal arguments can be given."
    sys.exit()

#####################
# histogram setting #
#####################

n_min = N_min/dV        # center value of left-most bin (=smallest value)
n_max = N_max/dV        # center value of right-most bin (=largest value)

if (nbin_in==0):        # in this case, consider all possible integers for n*dV
  nbin = N_max-N_min+1
else:                   # for any positive integer value
  nbin = nbin_in

dn = (n_max-n_min)/(nbin-1)        # width of each bin
bins = n_min+dn*np.arange(nbin)    # center values of bins 

bin_edges = list(bins-0.5*dn)      # note: len(bins) = nbin 
bin_edges.append(bins[-1]+0.5*dn)  #       len(bin_edges) = nbin + 1

###################################
# read data | calculate histogram #
###################################

# read data
f_data = []
with open(datafile) as inf:
  for line in inf:
    f_data.append(float(line.split()[1]))  # read the second number
f_data = np.array(f_data)

# calculate histogram (actually, density)
n_hist = np.histogram(f_data,bin_edges,density=True)[0]

###########################
# theoretical predictions #
###########################

# n values for discrete distribution (that is, Poisson)
nbin_disc = N_max-N_min+1
dn_disc = (n_max-n_min)/(nbin_disc-1)
bins_disc = n_min+dn_disc*np.arange(nbin_disc) 

# n values for continuous distributions (that is, Gaussian and Stirling approximation) 
nbin_cont = 10000
dn_cont = (n_max-n_min)/(nbin_cont-1)
bins_cont = n_min+dn_cont*np.arange(nbin_cont) 

# Poisson 

def fnc_Poisson(x):
  if (round(x*dV)>=0):
    return scipy.stats.poisson.pmf(round(x*dV),n_av*dV)
  else:
    return 0

Poisson_disc = np.array([ fnc_Poisson(x) for x in bins_disc ])
Poisson_disc_normalized = Poisson_disc/dn_disc 

# Gaussian

def fnc_Gaussian(x):
  n_std = math.sqrt(n_av*dV)/dV        
  return math.exp(-0.5*(x-n_av)**2/n_std**2)/math.sqrt(2*math.pi)/n_std

Gaussian = np.array([ fnc_Gaussian(x) for x in bins ])
Gaussian_cont = np.array([ fnc_Gaussian(x) for x in bins_cont ])

# Stirling's approximation to the Poisson distribution
# This includes a continuity correction to correct the mean

def fnc_Stirling(x):
  if (x>0):
    return math.exp(-(x+1/(2*dV))*(math.log((x+1/(2*dV))/n_av)-1)*dV);
  else:
    return 0

Stirling = np.array([ fnc_Stirling(x) for x in bins ])
Stirling = Stirling/(sum(Stirling)*dn)
Stirling_cont = np.array([ fnc_Stirling(x) for x in bins_cont ])
Stirling_cont = Stirling_cont/(sum(Stirling_cont)*dn_cont)

############################
# show plots on the screen #
############################

# semi-log scale plot

if (plot_semilogy):
  fig, a = plt.subplots()
  plt.yscale('log')
  a.plot(bins,n_hist,'--ok',label="Numerics") 
  a.plot(bins_disc,Poisson_disc_normalized,'sb',label="Poisson (normalized)",mfc='none')
  a.plot(bins_cont,Stirling_cont,'-g',label="Stirling approximation")
  a.plot(bins_cont,Gaussian_cont,'-r',label="Gaussian")
  a.legend(loc=8,numpoints=1,fontsize='small')
  #plt.title("N=%d particles per cell"%round(n_av*dV))
  plt.xlabel("number density n")
  plt.ylabel("probability density")
  #plt.vlines(0,1e-8,1,colors='k',linestyles='dotted')
  plt.show()

if (plot_linear):
  fig, a = plt.subplots()
  plt.yscale('linear')
  a.plot(bins,n_hist,'--ok',label="Numerics") 
  a.plot(bins_disc,Poisson_disc_normalized,'sb',label="Poisson (normalized)",mfc='none')
  a.plot(bins_cont,Stirling_cont,'-g',label="Stirling approximation")
  a.plot(bins_cont,Gaussian_cont,'-r',label="Gaussian")
  a.legend(loc=1,numpoints=1,fontsize='small')
  #plt.title("N=%d particles per cell"%round(n_av*dV))
  plt.xlabel("number density n")
  plt.ylabel("probability density")
  #plt.vlines(0,1e-8,1,colors='k',linestyles='dotted')
  plt.show()

################
# output files # 
################

# output_hist
if (output_hist!=""):
  out = open(output_hist,'w')
  out.write("# bin_val numerics Stirling Gaussian\n");

  for i in range(nbin):
    out.write("%g\t%g\t%g\t%g\n" % (bins[i], n_hist[i], Stirling[i], Gaussian[i]))
  out.close() 

  print "output_hist file \"%s\" generated." % output_hist

# output_cont
if (output_cont!=""):
  out = open(output_cont,'w')
  out.write("# bin_val Stirling Gaussian\n");

  for i in range(nbin_cont):
    out.write("%g\t%g\t%g\n" % (bins_cont[i], Stirling_cont[i], Gaussian_cont[i]))
  out.close() 

  print "output_cont file \"%s\" generated." % output_cont

# output_poiss
if (output_poiss!=""):
  out = open(output_poiss,'w')
  out.write("# bin_val Poisson(normalized)\n");

  for i in range(nbin_disc):
    out.write("%g\t%g\n" % (bins_disc[i], Poisson_disc_normalized[i]))
  out.close() 

  print "output_poiss file \"%s\" generated." % output_poiss
