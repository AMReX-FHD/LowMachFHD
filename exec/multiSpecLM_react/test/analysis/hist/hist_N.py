import sys
import math
import numpy as np
import scipy.stats 
import matplotlib.pyplot as plt

# this script computes the following statistics on the number of molecules
#  of a chosen species in a cell (= N):   
# 1. histogram of N (compared with Gaussian and Poisson statistics)
#   a. semi-log-y plot 
#   b. linear plot
# 2. mean and variance of N 
# 3. probability of negative density
# 4. KL-divergence with respect to the Poisson distribution
# screen outputs as well as file outputs are generated.
#
# usage:
#   python hist_N.py 
#   python hist_N.py input_column norm_factor Nav Nmin Nmax show_plot(yes/no) 
#     save_plot(yes/no)

##############
# parameters #
##############

# datafile contains (number or mass) densities of each cell at different times.
datafile = "fort.10"

# output filenames for histogram (Numeric,Poisson,Gaussian) and continuous distribution (Gaussian)
# if "none" is given, that file is omitted.  
output_hist = "res.hist" 
output_cont = "res.hist_cont" 
# output filenames for figures
output_fig_hist1 = "hist_semilogy.png"
output_fig_hist2 = "hist_linear.png"
# output filename for mean, variance, negative density probability, and the KL-divergence
output_N_stat = "res.N_stat"

# input_column indicates which column of densities is to be analyzed.
# note: c convention is used. e.g. first column = 0, second column = 1, ...
input_column = 1

# norm_factor is required to convert densities into the number of molecules.
# hence, this quantity corresponds to the value of density representing a single molecule.
norm_factor = 3.e-9

# Nav = expected mean of N
Nav = 10.

# Nmin/Nmax: min/max values of histogram bins
Nmin = -2
Nmax = 30

# flag for showing plots on the screen
show_plot = True
# flag for saving plots on a file 
save_plot = True

# if additional arguments are given
if (len(sys.argv)!=1):
  if (len(sys.argv)==8):                # for seven additonal arguments
    input_column = int(sys.argv[1])
    norm_factor = float(sys.argv[2])
    Nav = float(sys.argv[3])
    Nmin = int(sys.argv[4])
    Nmax = int(sys.argv[5])
    if (sys.argv[6]=="no"):
      show_plot = False
    if (sys.argv[7]=="no"):
      save_plot = False
  else:                                 # otherwise, abort 
    print "Error: only seven additonal arguments can be given."
    print "-> input_column norm_factor Nav Nmin Nmax show_plot(yes/no) save_plot(yes/no)"
    sys.exit()

#####################
# histogram setting #
#####################

nbin = Nmax-Nmin+1
bins = Nmin+np.arange(nbin)             # center values of bins 

bin_edges = list(bins-0.5)              # note: len(bins) = nbin 
bin_edges.append(bins[-1]+0.5)          #       len(bin_edges) = nbin + 1

#############
# read data #
#############

print "** reading %s..." % datafile

f_data = []
with open(datafile) as inf:
  for line in inf:
    # read the number in input_column and then convert it into the number of molecules
    val = float(line.split()[input_column])/norm_factor

    if (val<Nmin or val>Nmax):
      print "** Warning: val=%f is out of range [Nmin=%d:Nmax=%d]" % (val,Nmin,Nmax)

    f_data.append(val)                   
f_data = np.array(f_data)

#######################
# calculate histogram #
#######################

# actually, we calculate the probability density function.
Nhist = np.histogram(f_data,bin_edges,density=True)[0]

###########################
# theoretical predictions #
###########################

# Poisson 
def fnc_Poisson(x):
  if (round(x)>=0):
    return scipy.stats.poisson.pmf(round(x),Nav)
  else:
    return 0

Nhist_poiss = np.array([ fnc_Poisson(x) for x in bins ])

# Gaussian
def fnc_Gaussian(x):
  return math.exp(-0.5*(x-Nav)**2/Nav)/math.sqrt(2*math.pi*Nav)

Nhist_gauss = np.array([ fnc_Gaussian(x) for x in bins ])

# for continuous distribution (=Gaussian) 
nbin_cont = (Nmax-Nmin)*10+1
dN_cont = float(Nmax-Nmin)/(nbin_cont-1)
bins_cont = Nmin+dN_cont*np.arange(nbin_cont) 

Nhist_gauss_cont = np.array([ fnc_Gaussian(x) for x in bins_cont ])

############################
# show plots on the screen #
############################

# semi-log scale plot

fig, a = plt.subplots()
plt.yscale('log')

a.plot(bins,Nhist,':ob',label="Numerics") 
a.plot(bins,Nhist_poiss,'sb',label="Poisson",mfc='none')
a.plot(bins_cont,Nhist_gauss_cont,'--r',label="Gaussian")

a.legend(loc=8,numpoints=1,fontsize='small')
plt.title("N=%.2f molecules per cell"%Nav)
plt.xlabel("number of molecules N")
plt.ylabel("probability density")

plt.ylim(ymin=Nhist_poiss[-1]/10)

if (show_plot):
  plt.show()
if (save_plot):
  fig.savefig(output_fig_hist1)
  print "** %s generated." % output_fig_hist1

# linear scale plot

fig, a = plt.subplots()
plt.yscale('linear')

a.plot(bins,Nhist,':ob',label="Numerics") 
a.plot(bins,Nhist_poiss,'sb',label="Poisson",mfc='none')
a.plot(bins_cont,Nhist_gauss_cont,'--r',label="Gaussian")

a.legend(loc=1,numpoints=1,fontsize='small')
plt.title("N=%.2f molecules per cell"%Nav)
plt.xlabel("number of molecules N")
plt.ylabel("probability density")

if (show_plot):
  plt.show()
if (save_plot):
  fig.savefig(output_fig_hist2)
  print "** %s generated." % output_fig_hist2

##############################
# output files for histogram # 
##############################

# output_hist
if (output_hist!="none"):
  out = open(output_hist,'w')
  
  out.write("# bin_val numerics Poisson Gaussian\n")

  for i in range(nbin):
    out.write("%g\t%g\t%g\t%g\n" % (bins[i],Nhist[i],Nhist_poiss[i],Nhist_gauss[i]))

  out.close() 

  print "** %s generated." % output_hist

# output_cont
if (output_cont!="none"):
  out = open(output_cont,'w')

  out.write("# bin_val Gaussian\n")

  for i in range(nbin_cont):
    out.write("%g\t%g\n" % (bins_cont[i],Nhist_gauss_cont[i]))
  out.close() 

  print "** %s generated." % output_cont

#####################################################################
# calculate mean, variance, prob of neg dens, and the KL-divergence # 
#####################################################################

# mean and variance
data_mean = np.mean(f_data)
data_var = np.var(f_data)

print "<N>= %f" % data_mean 
print "Var[N]= %f" % data_var

# probability of negative density
neg_cnt = 0
for N in f_data:
  if (N<0.):
    neg_cnt += 1
data_prob_neg = float(neg_cnt)/len(f_data)

print "prob_of_neg_dens= %e" % data_prob_neg 

# Kullback-Leibler divergence
sum_KLdiv = 0.
for i in range(nbin):
  if (Nhist_poiss[i]>0. and Nhist[i]>0.):
    # D_KL(P||Q) = sum_i P_i*log(P_i/Q_i)
    # sum is over Q_i>0.
    # if P_i=0, its contribution becomes zero.
    # note that Poisson_disc is used for P_i, whereas Poisson_disc_normalized for P_i/Q_i. 
    sum_KLdiv += Nhist_poiss[i]*math.log(Nhist_poiss[i]/Nhist[i])

print "KL-divergence= %e" % sum_KLdiv

###########################################################################
# file output for mean, variance, prob of neg dens, and the KL-divergence #
########################################################################### 

out = open(output_N_stat,"w")
out.write("<N>= %f\n" % data_mean)
out.write("Var[N]= %f\n" % data_var)
out.write("prob_of_neg_dens= %e\n" % data_prob_neg)
out.write("KL-divergence= %e\n" % sum_KLdiv)
out.close()
print "** %s generated." % output_N_stat
