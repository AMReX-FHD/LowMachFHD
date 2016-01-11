import math
import numpy as np
import scipy.stats 
import matplotlib.pyplot as plt

u_av = 1.                                  # average concentration
dx = 5.                                    # grid size (cell volume in higher dimensions)
datafile = "fort.10"

###################
# histogram range #
###################

u_std = math.sqrt(u_av*dx)/dx              # expected standard deviation
u_std_int = round(math.sqrt(u_av*dx))/dx   # ensure that the bins are centered on actual integers

range1 = u_av-4*u_std_int                  # center value of left-most bin in histogram
range2 = u_av+6*u_std_int                  # center value of right-most bin in histogram 

Nbins = int(round((range2-range1)*dx)+1)   # number of bins for numerical results
du = (range2-range1)/(Nbins-1)             # width of each bin
bins = range1+du*np.arange(Nbins)          # center values of bins 

bin_edges = list(bins-0.5*du)              # note: len(bins) = Nbins 
bin_edges.append(bins[-1]+0.5*du)          #       len(bin_edges) = Nbins + 1

Nbins_th = 10000                           # number of bins for theoretical predictions
du_th = (range2-range1)/(Nbins_th-1)       # (compute the theory on a finer grid to get closer to zero)  
bins_th = range1+du_th*np.arange(Nbins_th)

###########################
# theoretical predictions #
###########################

# Gaussian approximation to the Poisson distribution
# Gaussian for bins_th | Gaussian0 for bins
Gaussian = np.exp(-(bins_th-u_av)**2/(2*u_std**2))/math.sqrt(2*math.pi)/u_std;
Gaussian0 = np.exp(-(bins-u_av)**2/(2*u_std**2))/math.sqrt(2*math.pi)/u_std;

# Stirling's approximation to the Poisson distribution
# This includes a continuity correction to correct the mean
# Stirling for bins_th | Stirling for bins
def fnc_Stirling(x):
  if (x>0):
    return math.exp(-(x+1/(2*dx))*(math.log((x+1/(2*dx))/u_av)-1)*dx);
  else:
    return 0

Stirling = np.array([ fnc_Stirling(x) for x in bins_th ])
Stirling = Stirling/(sum(Stirling)*du_th)
Stirling0 = np.array([ fnc_Stirling(x) for x in bins ])
Stirling0 = Stirling0/(sum(Stirling0)*du)

# True Poisson distribution
def fnc_Poisson(x):
  if (round(x*dx)>=0):
    return scipy.stats.poisson.pmf(round(x*dx),u_av*dx)
  else:
    return 0

Poisson = np.array([ fnc_Poisson(x) for x in bins ])
Poisson = Poisson/(sum(Poisson)*du) 

################################################################################
# read data file | calculate histogram | create plot | write histogram on file # 
################################################################################

compare = 0                        # whether to plot many curves at once or just a single run

if (compare):                      # compare difference curves (to be implemented later)

  print "The compare feature has not been implemented yet..."

else:                              # plot results from a single run

  # read data
  f_data = []
  with open('fort.10') as inf:
    for line in inf:
      parts = line.split()
      f_data.append(float(parts[1]))   # remove the time label
  f_data = np.array(f_data)

  # calculate histogram
  u_hist = np.histogram(f_data,bin_edges)[0]
  u_hist_std = np.sqrt(u_hist)

  Z = 1/(sum(u_hist)*du)
  u_hist = Z*u_hist
  u_hist_std = Z*u_hist_std

  u_mean_num = np.mean(f_data)
  u_std_num = np.std(f_data)

  # create plot on screen
  fig, a = plt.subplots()
  plt.yscale('log')
  a.plot(bins,u_hist,'ok',label="Numerics",mfc='none') 
  a.plot(bins,Poisson,'-b',label="Poisson")
  a.plot(bins_th,Stirling,'--g',label="Stirling")
  a.plot(bins_th,Gaussian,'--r',label="Gaussian")
  a.legend(loc=3,numpoints=1)
  plt.title("N=%d particles per cell"%round(u_av*dx))
  plt.show()

  # write on file

  outputfile = datafile+"_hist"
  out = open(outputfile,'w')
  out.write("# bin_val numerics Poisson Stirling Gaussian\n");

  for i in range(Nbins):
    out.write("%g\t%g\t%g\t%g\t%g\n" % (bins[i], u_hist[i], Poisson[i], Stirling0[i], Gaussian0[i]))
  out.close() 

  print "histogram has been saved in %s" % outputfile
