import matplotlib.pyplot as plt
import numpy as np

########
# hist #
########

## read data

DATA1 = "res.hist_stat"
DATA2 = "res.hist_poiss"
DATA3 = "res.hist_cont"
DATA4 = "res.Sk_stat"
DATA5 = "res.hist_near_zero_stat"

OUTPUT1 = "hist_semilogy.png"
OUTPUT2 = "hist_linear.png"
OUTPUT3 = "Sk.png"
OUTPUT4 = "hist_near_zero.png"

x1,y1,y1err2 = np.loadtxt(DATA1,unpack=True,usecols=[0,1,2])
y1err2 *= 2

x2,y2 = np.loadtxt(DATA2,unpack=True,usecols=[0,1],comments='#')

x3,y3 = np.loadtxt(DATA3,unpack=True,usecols=[0,1],comments='#')

x4,y4,y4err2 = np.loadtxt(DATA4,unpack=True,usecols=[0,1,2])
y4err2 *= 2 

x5,y5,y5err2 = np.loadtxt(DATA5,unpack=True,usecols=[0,1,2])
y5err2 *= 2

## hist: semi-log plot

fig = plt.figure()
ax = fig.gca()

ax.set_yscale('log')

if (len(y1)<=100):
  ax.errorbar(x1,y1,yerr=y1err2,label="Numerics",fmt='--o')
else:
  ax.errorbar(x1,y1,yerr=y1err2,label="Numerics",fmt='k')

if (len(y2)<=100):
  ax.plot(x2,y2,'bs',label="Poisson (normalized)",mfc='none')
else:
  ax.plot(x2,y2,'b-',label="Poisson (normalized)",mfc='none')

ax.plot(x3,y3,'r-',label="Gaussian")

ax.legend(loc=8,numpoints=1,fontsize='small')
ax.set_xlabel("number density n")
ax.set_ylabel("probability density")

plt.show()
fig.savefig(OUTPUT1)

## hist: linear plot

fig = plt.figure()
ax = fig.gca()

ax.set_yscale('linear')

if (len(y1)<=100):
  ax.errorbar(x1,y1,yerr=y1err2,label="Numerics",fmt='--o')
else:
  ax.errorbar(x1,y1,yerr=y1err2,label="Numerics",fmt='k')

if (len(y2)<=100):
  ax.plot(x2,y2,'bs',label="Poisson (normalized)",mfc='none')
else:
  ax.plot(x2,y2,'b-',label="Poisson (normalized)",mfc='none')

ax.plot(x3,y3,'r-',label="Gaussian")

ax.legend(loc=1,numpoints=1,fontsize='small')
ax.set_xlabel("number density n")
ax.set_ylabel("probability density")

plt.show()
fig.savefig(OUTPUT2)

## S(k)

fig = plt.figure()
ax = fig.gca()

ax.errorbar(x4,y4,yerr=y4err2,fmt='s')

ax.set_xlabel("k")
ax.set_ylabel("structure factor S(k)")

plt.show()
fig.savefig(OUTPUT3)

## hist_near_zero

fig = plt.figure()
ax = fig.gca()

ax.errorbar(x5,y5,yerr=y5err2,fmt='x')

ax.set_xlabel("number density n")
ax.set_ylabel("probability density")
ax.set_ylim(bottom=0)

plt.show()
fig.savefig(OUTPUT4)
