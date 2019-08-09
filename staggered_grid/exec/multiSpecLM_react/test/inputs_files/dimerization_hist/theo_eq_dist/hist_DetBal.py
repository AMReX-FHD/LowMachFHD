import numpy as np

N0 = 40
K = 0.7755

PP = np.ones(N0/2+1)
NN2 = np.arange(N0/2+1)

for i in range(N0/2):
  N2 = NN2[i] 
  PP[i+1] = PP[i]*K*(N0-2*N2)*(N0-2*N2-1)/(N0-N2)/(N2+1)

PP /= sum(PP)

print "K=%f\t<N2>= %e" % (K,sum(PP*NN2))

np.savetxt("res.hist_DetBal",np.transpose([NN2,PP]))
