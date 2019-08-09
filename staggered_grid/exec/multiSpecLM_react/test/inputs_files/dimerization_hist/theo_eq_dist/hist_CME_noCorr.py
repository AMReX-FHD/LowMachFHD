import numpy as np

N0 = 40
#K = 0.7755     # from detailed balance
K = 0.7625      # without correction, this gives N2eq=10

PP = np.ones(N0/2+1)
NN2 = np.arange(N0/2+1)

def x1sq(n0,n2):
  tmp = (n0-2.*n2)/(n0-n2)
  tmp *= tmp                            # no correction
#  tmp *= (n0-2.*n2-1)/(n0-n2-1)        # correction -> detailed balance
  return tmp

for i in range(N0/2):
  N2 = NN2[i] 
  PP[i+1] = PP[i]*K*x1sq(N0,N2)*(N0-N2-1)/(N2+1)

PP /= sum(PP)

print "K=%f\t<N2>= %e" % (K,sum(PP*NN2))

np.savetxt("res.hist_CME_noCorr_K%.4f" % K,np.transpose([NN2,PP]))
