import numpy as np
import matplotlib.pyplot as plt

num_pnts = 1000
T = 4.0e-4 
phi_0 = 34.
xvals = np.linspace(0., 8.*T, num_pnts, endpoint=True)
yvals = np.zeros(num_pnts)
yvals2 = np.zeros(num_pnts)
yvals3 = np.zeros(num_pnts)
def beta(x):
	return x - np.floor(x/T)*T

def ac(x):
	omega = 0.05*T
	#return -1.*phi_0*np.tanh((x-T/2)/omega)
	return -1.*phi_0*np.tanh((x-T/2)/omega)

def ac2(x):
	if int(np.floor(x/T)) % 2 == 0:
		return ac(x - np.floor(x/T)*T)
	else: 
		return -1.*ac(x - np.floor(x/T)*T)

 
for i in range(num_pnts): 
	yvals[i] = ac2(xvals[i])
	#yvals2[i] = -1.*ac(xvals[i]-T)
	#yvals3[i] = ac(xvals[i]-2.*T)


charge_time = 1.e-4

plt.figure(1) 
plt.plot(xvals, yvals, 'b-o')
plt.plot(np.ones(100)*charge_time, np.linspace(0.0, phi_0, 100), 'g-o')

plt.show()
