import sys
import matplotlib.pyplot as plt
import numpy as np

colperline = 4
figfile = "test.png"

# check command-line input

if len(sys.argv)<4:
  print "ERROR: python this_script colperline figfile fit_data1 fit_data2 ..."
  exit(0)

colperline = int(sys.argv[1])
figfile = sys.argv[2]

print "colperline=%d"%colperline
print "figfile=%s"%figfile

nfitdata = len(sys.argv)-3

fitdata = []
for i in range(nfitdata):
  fitdata.append(sys.argv[i+3])

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

plt.rc('font',family="serif")
plt.rc('font',size=10)

fig = plt.figure()

for n in range(nparam):
  ax = fig.add_subplot(4,2,n+1)

  if n+1 != 8:
    tmp="$a_%i$"%(n+1)
    ax.set_title(r"%s"%tmp,fontsize=10)
  else:
    ax.set_title("RMS Fit Error",fontsize=10)

  if n+1 == 2:
    prevn = 2
  elif n+1 == 3:
    prevn = 1
  else:
    prevn = n

  ax.set_xlim([0.5,nfitdata+0.5])
  ax.set_xticks(range(1,nfitdata+1))

  for i in range(nfitdata):
    nsample = len(aaa[i])
    xx = (i+1)*np.ones(nsample)
 
    yy = []
    for j in range(nsample):
      yy.append(aaa[i][j][prevn]) 

    ax.plot(xx,yy,'o')

legend="$\\left(1-\\tanh\\frac{t-a_1}{a_2}\\right)\\left(a_3\\sin(a_4 t+a_5)+a_6\\right)+a_7$\n"
for i in range(nfitdata):
  tmp = "$%i:$ %s  " % (i+1,fitdata[i])
  legend+=tmp
  if ((i+1)%colperline==0):
    legend+="\n"
fig.suptitle(r"%s"%legend,fontsize=10)

fig.set_size_inches(10,10)
fig.tight_layout()
fig.subplots_adjust(top=0.8)
fig.savefig(figfile,dpi=200)
