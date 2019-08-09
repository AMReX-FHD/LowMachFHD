import matplotlib.pyplot as plt
import numpy as np

########
# hist #
########

## read data

DATA1 = "cell.hist_stat"
DATA2 = "cell.hist_poiss"
DATA3 = "cell.hist_cont"
#DATA4 = "cell.Sk_stat"
DATA5 = "cell.hist_near_zero_stat"

OUTPUT1 = "cell_hist_semilogy.png"
OUTPUT2 = "cell_hist_linear.png"
OUTPUT3 = "cell_Sk.png"
OUTPUT4 = "cell_hist_near_zero.png"

x1,y1,y1err2 = np.loadtxt(DATA1,unpack=True,usecols=[0,1,2])
y1err2 *= 2

x2,y2 = np.loadtxt(DATA2,unpack=True,usecols=[0,1],comments='#')

x3,y3 = np.loadtxt(DATA3,unpack=True,usecols=[0,1],comments='#')

#x4,y4,y4err2 = np.loadtxt(DATA4,unpack=True,usecols=[0,1,2])
#y4err2 *= 2 

x5,y5,y5err2 = np.loadtxt(DATA5,unpack=True,usecols=[0,1,2])
y5err2 *= 2

## hist: semi-log plot

fig = plt.figure()
ax = fig.gca()

ax.set_yscale('log')

if (len(y1)<=100):
  (_, caps, _) = ax.errorbar(x1,y1,yerr=y1err2,label="Numerics",fmt='ro',mfc='none',mec='r')
  for cap in caps:
     cap.set_markersize(10)
     cap.set_markeredgewidth(2)
else:
  ax.errorbar(x1,y1,yerr=y1err2,label="Numerics",fmt='r')

if (len(y2)<=100):
  ax.plot(x2,y2,'bs',label="Poisson (normalized)",mfc='b')
else:
  ax.plot(x2,y2,'b-',label="Poisson (normalized)",mfc='none')

ax.plot(x3,y3,'k:',label="Gaussian")

ax.legend(loc=8,numpoints=1,fontsize='small')
ax.set_ylim([1e-6,10])
ax.set_xlabel("number density n")
ax.set_ylabel("probability density")

plt.show()
fig.savefig(OUTPUT1)

## hist: linear plot

fig = plt.figure()
ax = fig.gca()

ax.set_yscale('linear')

if (len(y1)<=100):
  (_, caps, _) = ax.errorbar(x1,y1,yerr=y1err2,label="Numerics",fmt='ro',mfc='none',mec='r')
  for cap in caps:
     cap.set_markersize(10)
     cap.set_markeredgewidth(2)
else:
  ax.errorbar(x1,y1,yerr=y1err2,label="Numerics",fmt='r')

if (len(y2)<=100):
  ax.plot(x2,y2,'bs',label="Poisson (normalized)",mfc='b')
else:
  ax.plot(x2,y2,'b-',label="Poisson (normalized)",mfc='none')

ax.plot(x3,y3,'k:',label="Gaussian")

ax.legend(loc=1,numpoints=1,fontsize='small')
ax.set_ylim([0,1])
ax.set_xlabel("number density n")
ax.set_ylabel("probability density")

plt.show()
fig.savefig(OUTPUT2)

## S(k)

#fig = plt.figure()
#ax = fig.gca()

#ax.errorbar(x4,y4,yerr=y4err2,fmt='s')

#ax.set_xlabel("k")
#ax.set_ylabel("structure factor S(k)")

#plt.show()
#fig.savefig(OUTPUT3)

## hist_near_zero

fig = plt.figure()
ax = fig.gca()

ax.errorbar(x5,y5,yerr=y5err2,fmt='x')

ax.set_xlabel("number density n")
ax.set_ylabel("probability density")
ax.set_ylim([0,0.6])
#ax.set_ylim(bottom=0)

plt.show()
fig.savefig(OUTPUT4)
