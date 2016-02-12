import sys
import math
import numpy as np
import scipy.stats 
import matplotlib.pyplot as plt

# this script calculates the histogram of n as well as the mean and variance of n.
# screen output as well as file output are generated.

# usage:
#   python hist_n.py 
#   python hist_n.py n_av dV show_plot(yes/no) save_plot(yes/no)

##############
# parameters #
##############

# datafile contains number densities n of the cells at different times.
datafile = "fort.10"

# output filenames for (n_hist,Gauss,Stirling), (Gaussian,Stirling), or Poisson
# if "none" is given, that file is omitted.  
output_hist = "res.hist" 
output_cont = "res.hist_cont" 
output_poiss = "res.hist_poiss" 
# output filenames for figures
output_fig_hist1 = "hist_semilogy.png"
output_fig_hist2 = "hist_linear.png"
# output filename for mean and variance
output_n_stat = "res.n_stat"

# average value of n (n_av) and the volume of a cell (dV) are needed
#  to calculate Gaussian and Poisson distributions. 
# for 1d and 2d, dV=dx*dy*cross_section.
n_av = 1. 
dV = 10.

# flags for showing semi-log-scale or linear plot on the screen
show_plot = True
save_plot = True

# if additional arguments are given
if (len(sys.argv)!=1):
  if (len(sys.argv)==5):    # for four additonal arguments
    n_av = float(sys.argv[1])
    dV = float(sys.argv[2])
    if (sys.argv[3]=="no"):
      show_plot = False
    if (sys.argv[4]=="no"):
      save_plot = False
  else:                     # otherwise, generate error
    print "Error: only four additonal arguments can be given."
    sys.exit()

# the range of histogram [n_min,n_max] is determined by [N_min/dV,N_max/dV].
# number of bins (nbin) is determined by nbin_in:
#  if nbin_in=0, nbin = N_max-N_min+1. 
#  otherwise, nbin = nbin_in.
N_min = int(n_av*dV-6.*math.sqrt(n_av*dV))
N_max = int(n_av*dV+7.*math.sqrt(n_av*dV)) 
nbin_in = 0

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
    nval = float(line.split()[1])          # read the second number

    if (nval<n_min or nval>n_max):
      print "** Warning: nval=%g is out of range [n_min=%g:n_max=%g]" % (nval,n_min,n_max)

    f_data.append(nval)                   
f_data = np.array(f_data)

# calculate histogram (actually, density)
n_hist = np.histogram(f_data,bin_edges,density=True)[0]

data_mean = np.mean(f_data)
data_var = np.var(f_data)

print "<n>= %f" % data_mean 
print "Var[n]= %f" % data_var

out = open(output_n_stat,"w")
out.write("<n>= %f\n" % data_mean)
out.write("Var[n] = %f\n" % data_var)
out.close()
print "%s generated." % output_n_stat

###########################
# theoretical predictions #
###########################

# n values for discrete distribution (that is, Poisson)
nbin_disc = N_max-N_min+1
dn_disc = (n_max-n_min)/(nbin_disc-1)
bins_disc = n_min+dn_disc*np.arange(nbin_disc) 

# n values for continuous distributions (that is, Gaussian and Stirling approximation) 
nbin_cont = (N_max-N_min)*10+1
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

overflowerror_Stirling = False

def fnc_Stirling(x):
  global overflowerror_Stirling
  if (overflowerror_Stirling):
    return 0.
  elif (x>0):
    try:
      return math.exp(-(x+1/(2*dV))*(math.log((x+1/(2*dV))/n_av)-1)*dV)
    except OverflowError:
      overflowerror_Stirling = True
      print "** Warning: overflow error in fnc_Stirling"
      return 0.
  else:
    return 0.

Stirling = np.array([ fnc_Stirling(x) for x in bins ])
Stirling_cont = np.array([ fnc_Stirling(x) for x in bins_cont ])

if (not overflowerror_Stirling):
  Stirling = Stirling/(sum(Stirling)*dn)
  Stirling_cont = Stirling_cont/(sum(Stirling_cont)*dn_cont)

############################
# show plots on the screen #
############################

# semi-log scale plot

fig, a = plt.subplots()
plt.yscale('log')
a.plot(bins,n_hist,'--ok',label="Numerics") 
a.plot(bins_disc,Poisson_disc_normalized,'sb',label="Poisson (normalized)",mfc='none')

if (not overflowerror_Stirling):
  a.plot(bins_cont,Stirling_cont,'-g',label="Stirling approximation")

a.plot(bins_cont,Gaussian_cont,'-r',label="Gaussian")
a.legend(loc=8,numpoints=1,fontsize='small')
#plt.title("N=%d particles per cell"%round(n_av*dV))
plt.xlabel("number density n")
plt.ylabel("probability density")
#plt.vlines(0,1e-8,1,colors='k',linestyles='dotted')
if (show_plot):
  plt.show()
if (save_plot):
  fig.savefig(output_fig_hist1)
  print "%s generated." % output_fig_hist1

fig, a = plt.subplots()
plt.yscale('linear')
a.plot(bins,n_hist,'--ok',label="Numerics") 
a.plot(bins_disc,Poisson_disc_normalized,'sb',label="Poisson (normalized)",mfc='none')

if (not overflowerror_Stirling):
  a.plot(bins_cont,Stirling_cont,'-g',label="Stirling approximation")

a.plot(bins_cont,Gaussian_cont,'-r',label="Gaussian")
a.legend(loc=1,numpoints=1,fontsize='small')
#plt.title("N=%d particles per cell"%round(n_av*dV))
plt.xlabel("number density n")
plt.ylabel("probability density")
#plt.vlines(0,1e-8,1,colors='k',linestyles='dotted')
if (show_plot):
  plt.show()
if (save_plot):
  fig.savefig(output_fig_hist2)
  print "%s generated." % output_fig_hist2

################
# output files # 
################

# output_hist
if (output_hist!="none"):

  out = open(output_hist,'w')
  
  if (overflowerror_Stirling):
    out.write("# bin_val numerics Gaussian\n")
  else:
    out.write("# bin_val numerics Gaussian Poisson(Stirling approx)\n")

  for i in range(nbin):
    if (overflowerror_Stirling):
      out.write("%g\t%g\t%g\n" % (bins[i],n_hist[i],Gaussian[i]))
    else:
      out.write("%g\t%g\t%g\t%g\n" % (bins[i],n_hist[i],Gaussian[i], Stirling[i]))

  out.close() 

  print "output_hist file \"%s\" generated." % output_hist

# output_cont
if (output_cont!="none"):
  out = open(output_cont,'w')

  if (overflowerror_Stirling):
    out.write("# bin_val Gaussian\n")
  else:
    out.write("# bin_val Gaussian Poisson(Stirling approx)\n")

  for i in range(nbin_cont):
    if (overflowerror_Stirling):
      out.write("%g\t%g\n" % (bins_cont[i],Gaussian_cont[i]))
    else:
      out.write("%g\t%g\t%g\n" % (bins_cont[i],Gaussian_cont[i],Stirling_cont[i]))
  out.close() 

  print "output_cont file \"%s\" generated." % output_cont

# output_poiss
if (output_poiss!="none"):
  out = open(output_poiss,'w')
  out.write("# bin_val Poisson(normalized)\n");

  for i in range(nbin_disc):
    out.write("%g\t%g\n" % (bins_disc[i], Poisson_disc_normalized[i]))
  out.close() 

  print "output_poiss file \"%s\" generated." % output_poiss
