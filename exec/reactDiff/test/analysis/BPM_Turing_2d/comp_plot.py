import sys
import matplotlib.pyplot as plt
import numpy as np

# check command-line input

if len(sys.argv)<3:
  print "ERROR: python this_script fit_data1 fit_data2 ..."
  exit(0)

nfitdata = len(sys.argv)-1

fitdata = []
for i in range(nfitdata):
  fitdata.append(sys.argv[i+1])

# read data

nparam = 8

aaa = []
for i in range(nfitdata):
  filename = fitdata[i]+".fit"
  print "reading %s" % filename 

  aa = []
  with open(filename) as inf:
    for line in inf:
      linesplit=line.split()

      if len(linesplit)!=nparam:
        print "ERROR: number of parameters should be 8 (including fitting error)"

      a = []
      for n in range(len(linesplit)):
        val = float(linesplit[n])
        a.append(val)
  
      aa.append(a)
 
  aaa.append(aa)

# plot

fig = plt.figure(figsize=(12,12),dpi=96)

for n in range(nparam):
  ax = fig.add_subplot(4,2,n+1)

  if n!=7:
    ax.set_title("a%i"%n)
  else:
    ax.set_title("RMS Fit Error")

  ax.set_xlim([0.5,nfitdata+0.5])
  ax.set_xticks(range(1,nfitdata+1))

  for i in range(nfitdata):
    nsample = len(aaa[i])
    xx = (i+1)*np.ones(nsample)
 
    yy = []
    for j in range(nsample):
      yy.append(aaa[i][j][n]) 

    ax.plot(xx,yy,'o')

legend=""
for i in range(nfitdata):
  tmp = "%i: %s      " % (i+1,fitdata[i])
  legend+=tmp
  if ((i+1)%4==0):
    legend+="\n"
fig.suptitle(legend)

plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.show()
