import numpy as np

data="D_cs5_dx1_Nc512_EulerTau_avg12_dt0.001.hist_near_zero_stat"
output=data+"_intg"

[nn,hh] = np.loadtxt(data,unpack=True,usecols=[0,1])

dx = nn[1]-nn[0]

ss = []
for i in range(len(nn)):
  if (i==0):
    ss.append(dx*hh[i])
  else:
    ss.append(ss[i-1]+dx*hh[i])
ss = np.array(ss)

np.savetxt(output,np.transpose([nn,hh,ss]))
