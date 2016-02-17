import matplotlib.pyplot as plt
import numpy as np
import sys

# if additional argument is given, normalize the sample std by the number
if (len(sys.argv)==2):
  NRUN = int(sys.argv[1])
else:
  NRUN = 1

########
# hist #
########

DATA1 = "res.hist_stat"
DATA2 = "res.hist_poiss"
DATA3 = "res.hist_cont"

x1 = []
y1 = []
y1err2 = []
with open(DATA1) as inf:
  for line in inf:
    parts = line.split()
    x1.append(float(parts[0]))
    y1.append(float(parts[1]))
    y1err2.append(float(parts[2]))
y1err2=np.array(y1err2)
y1err2=2*np.sqrt(y1err2/NRUN)  # normalize by NRUN (errorbars with 2 std)

x2 = []
y2 = []
lcnt = 0
with open(DATA2) as inf:
  for line in inf:
    lcnt += 1
    if (lcnt<=1):
      continue
    parts = line.split()
    x2.append(float(parts[0]))
    y2.append(float(parts[1]))

x3 = []
y31 = []
y32 = []
lcnt = 0
with open(DATA3) as inf:
  for line in inf:
    lcnt += 1
    if (lcnt<=1):
      continue
    parts = line.split()
    x3.append(float(parts[0]))
    if (len(parts)==3):
      y31.append(float(parts[1]))
      y32.append(float(parts[2]))
    else:
      y31.append(float(parts[1]))

# semi-log plot
fig, a = plt.subplots()
plt.yscale('log')

if (len(y1)<=100):
  a.errorbar(x1,y1,yerr=y1err2,label="Numerics",fmt='--o')
else:
  a.errorbar(x1,y1,yerr=y1err2,label="Numerics",fmt='k')

if (len(y2)<=100):
  a.plot(x2,y2,'sb',label="Poisson (normalized)",mfc='none')
else:
  a.plot(x2,y2,'-b',label="Poisson (normalized)",mfc='none')

if (len(y32)>0):
  a.plot(x3,y32,'-g',label="Stirling approximation")

a.plot(x3,y31,'-r',label="Gaussian")
a.legend(loc=8,numpoints=1,fontsize='small')
#plt.title("N=%d particles per cell"%round(n_av*dV))
plt.xlabel("number density n")
plt.ylabel("probability density")
#plt.vlines(0,1e-8,1,colors='k',linestyles='dotted')
plt.show()
fig.savefig("hist_semilogy.png")

# linear plot
fig, a = plt.subplots()
plt.yscale('linear')

if (len(y1)<=100):
  a.errorbar(x1,y1,yerr=y1err2,label="Numerics",fmt='--o')
else:
  a.errorbar(x1,y1,yerr=y1err2,label="Numerics",fmt='k')

if (len(y2)<=100):
  a.plot(x2,y2,'sb',label="Poisson (normalized)",mfc='none')
else:
  a.plot(x2,y2,'-b',label="Poisson (normalized)",mfc='none')

if (len(y32)>0):
  a.plot(x3,y32,'-g',label="Stirling approximation")

a.plot(x3,y31,'-r',label="Gaussian")
a.legend(loc=1,numpoints=1,fontsize='small')
#plt.title("N=%d particles per cell"%round(n_av*dV))
plt.xlabel("number density n")
plt.ylabel("probability density")
#plt.vlines(0,1e-8,1,colors='k',linestyles='dotted')
plt.show()
fig.savefig("hist_linear.png")

######
# Sk #
######

DATA4 = "res.Sk_stat"

x4 = []
y4 = []
y4err2 = []
with open(DATA4) as inf:
  for line in inf:
    parts = line.split()
    x4.append(float(parts[0]))
    y4.append(float(parts[1]))
    y4err2.append(float(parts[2]))
y4err2=np.array(y4err2)
y4err2=2*np.sqrt(y4err2/NRUN)  # normalize by NRUN (errorbars with 2 std)

# plot
fig, a = plt.subplots()
a.errorbar(x4,y4,yerr=y4err2,label="Numerics",fmt='s')
#a.legend(loc=1,numpoints=1,fontsize='small')
#plt.title("N=%d particles per cell"%round(n_av*dV))
plt.xlabel("k")
plt.ylabel("structure factor S(k)")
plt.show()
fig.savefig("Sk.png")
