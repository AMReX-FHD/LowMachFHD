import numpy as np
import scipy.optimize as optimization
import matplotlib.pyplot as plt

datafile="fort.21"
factor=32*32*0.5
figfile="n_avg.png"
reportfile="res.fit"

# read data
xdata=[]
ydata=[]
with open(datafile) as inf:
  for line in inf:
    parts = line.split()
    xdata.append(float(parts[0]))
    ydata.append(float(parts[2])/factor)
xdata=np.array(xdata)
ydata=np.array(ydata)

# initial guess
init_a = np.array([2000,250,1000,0.01,3,75,1400])

# fitting function
def f(x,a0,a1,a2,a3,a4,a5,a6):
  return (1-np.tanh((x-a0)/a2))*(a1*np.sin(a3*x+a4)+a5)+a6

# fitting
[popt,pcov] = optimization.curve_fit(f,xdata,ydata,init_a,None)

# resulting function
def f_fitted(x):
  return (1-np.tanh((x-popt[0])/popt[2]))*(popt[1]*np.sin(popt[3]*x+popt[4])+popt[5])+popt[6]

# data from the fit 
xdata2 = np.arange(xdata[0],xdata[-1],0.1)
yfit = map(f_fitted,xdata)
yfit2 = map(f_fitted,xdata2)

# fitting error (root-mean-squared error)
err = np.sqrt(sum((ydata-yfit)**2)/len(ydata))

# plot
plt.plot(xdata,ydata,"r-")
plt.plot(xdata2,yfit2,"b--")
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.title(r'Root-Mean-Squared Error $\sqrt{\frac{1}{N}\sum{(y_i-\hat{y}_i)^2}}$ = %.2f' % err)
#plt.show()
plt.savefig(figfile)
print "%s generated" % figfile

# report
out = open(reportfile,"w")
out.write("%g\t" % popt[0])
out.write("%g\t" % popt[1])
out.write("%g\t" % popt[2])
out.write("%g\t" % popt[3])
out.write("%g\t" % popt[4])
out.write("%g\t" % popt[5])
out.write("%g\n" % popt[6])
out.close()
print "%s generated" % reportfile
