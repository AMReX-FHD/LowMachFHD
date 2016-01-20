import matplotlib.pyplot as plt
import numpy as np
import sys

DATA1 = "res.hist_stat"
DATA2 = "res.hist_poiss"
DATA3 = "res.hist_cont"

# if additional argument is given, normalize the sample std by the number
if (len(sys.argv)==2):
  NRUN = int(sys.argv[1])
else:
  NRUN = 1

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
    y31.append(float(parts[1]))
    y32.append(float(parts[2]))

# semi-log plot
fig, a = plt.subplots()
plt.yscale('log')
a.errorbar(x1,y1,yerr=y1err2,label="Numerics",fmt='--o')
a.plot(x2,y2,'sb',label="Poisson (normalized)",mfc='none')
a.plot(x3,y32,'-g',label="Stirling approximation")
a.plot(x3,y31,'-r',label="Gaussian")
a.legend(loc=8,numpoints=1,fontsize='small')
#plt.title("N=%d particles per cell"%round(n_av*dV))
plt.xlabel("number density n")
plt.ylabel("probability density")
#plt.vlines(0,1e-8,1,colors='k',linestyles='dotted')
plt.show()

# linear plot
fig, a = plt.subplots()
plt.yscale('linear')
a.errorbar(x1,y1,yerr=y1err2,label="Numerics",fmt='--o')
a.plot(x2,y2,'sb',label="Poisson (normalized)",mfc='none')
a.plot(x3,y32,'-g',label="Stirling approximation")
a.plot(x3,y31,'-r',label="Gaussian")
a.legend(loc=1,numpoints=1,fontsize='small')
#plt.title("N=%d particles per cell"%round(n_av*dV))
plt.xlabel("number density n")
plt.ylabel("probability density")
#plt.vlines(0,1e-8,1,colors='k',linestyles='dotted')
plt.show()
