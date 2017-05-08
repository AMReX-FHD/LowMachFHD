import matplotlib.pyplot as plt
import numpy as np

########
# hist #
########

## read data

DATA1 = "res.hist_stat"                 # numerics
DATA2 = "res.hist1"                     # Poisson (third column)
DATA3 = "res.hist_cont"                 # Gaussian
DATA4 = "res.hist_near_zero_stat"

OUTPUT1 = "hist_semilogy.png"
OUTPUT2 = "hist_linear.png"
OUTPUT3 = "hist_near_zero.png"

x1,y1,y1err2 = np.loadtxt(DATA1,unpack=True,usecols=[0,1,2])
y1err2 *= 2

x2,y2 = np.loadtxt(DATA2,unpack=True,usecols=[0,2],comments='#')

x3,y3 = np.loadtxt(DATA3,unpack=True,usecols=[0,1],comments='#')

x4,y4,y4err2 = np.loadtxt(DATA4,unpack=True,usecols=[0,1,2])
y4err2 *= 2

## hist: semi-log plot

fig = plt.figure()
ax = fig.gca()

ax.set_yscale('log')

ax.errorbar(x1,y1,yerr=y1err2,fmt=':ob',label="Numerics")
ax.plot(x2,y2,'sb',mfc='none',label="Poisson")
ax.plot(x3,y3,'--r',label="Gaussian")

ax.legend(loc=8,numpoints=1,fontsize='small')
ax.set_xlabel("Number of molecules N")
ax.set_ylabel("probability density")

ax.set_ylim(ymin=y2[-1]/10)

plt.show()
fig.savefig(OUTPUT1)

## hist: linear plot

fig = plt.figure()
ax = fig.gca()

ax.set_yscale('linear')

ax.errorbar(x1,y1,yerr=y1err2,fmt=":ob",label="Numerics")
ax.plot(x2,y2,'sb',mfc='none',label="Poisson")
ax.plot(x3,y3,'--r',label="Gaussian")

ax.legend(loc=1,numpoints=1,fontsize='small')
ax.set_xlabel("Number of molecules N")
ax.set_ylabel("probability density")

ax.set_ylim(ymin=0.)

plt.show()
fig.savefig(OUTPUT2)

## hist_near_zero

fig = plt.figure()
ax = fig.gca()

ax.errorbar(x4,y4,yerr=y4err2,fmt='x')

ax.set_xlabel("number density n")
ax.set_ylabel("probability density")
ax.set_ylim(bottom=0)

plt.show()
fig.savefig(OUTPUT3)
