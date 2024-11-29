import os
import urllib2
import numpy as np
import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import pylab


##==========================================
## Get time averaged profiles
##

#debye_len_case = '1_50'
Ly = 1.28e-4
debye_len_case = '1_10'

if (debye_len_case == '1_50'): 
	y_vals, u, v, w, c1, c2, c3, epot = np.loadtxt('hstat00064000_debye_len_1_50', unpack = True, usecols = [0,1,2,3,5,6,7,9])
elif (debye_len_case == '1_10'): 
	y_vals, u, v, w, c1, c2, c3, epot = np.loadtxt('hstat00056500_debye_len_1_10', unpack = True, usecols = [0,1,2,3,5,6,7,9])

u_det = np.loadtxt('u_deterministic_' + debye_len_case + '.txt', unpack = True, usecols = [0])
c1_det = np.loadtxt('c1_deterministic_' + debye_len_case + '.txt', unpack = True, usecols = [0])
c2_det = np.loadtxt('c2_deterministic_' + debye_len_case + '.txt', unpack = True, usecols = [0])
c3_det = np.loadtxt('c3_deterministic_' + debye_len_case + '.txt', unpack = True, usecols = [0])
epot_det = np.loadtxt('epot_deterministic_' + debye_len_case + '.txt', unpack = True, usecols = [0])


##==========================================
## Plot results 
##

plt.figure(1) 
plt.plot(y_vals, u, 'b-o', label = '3d averaged')
plt.plot(y_vals, u_det, 'r-o', label = '2d')
plt.title('$u$ v. wall-normal $y$')
plt.legend(loc='best')

plt.figure(7) 
plt.plot(y_vals, v, 'b-o', label = '3d averaged')
plt.title('$v$ v. wall-normal $y$')
plt.legend(loc='best')

plt.figure(9) 
plt.plot(y_vals, w, 'b-o', label = '3d averaged')
plt.title('$w$ v. wall-normal $y$')
plt.legend(loc='best')

plt.figure(2) 
plt.plot(y_vals, c1, 'b-o', label = '3d averaged')
plt.plot(y_vals, c1_det, 'r-o', label = '2d')
plt.title('$c_1$ v. wall-normal $y$')
plt.legend(loc='best')

plt.figure(3) 
plt.plot(y_vals, c2, 'b-o', label = '3d averaged')
plt.plot(y_vals, c2_det, 'r-o', label = '2d')
plt.title('$c_2$ v. wall-normal $y$')
plt.legend(loc='best')

plt.figure(4) 
plt.plot(y_vals, c3, 'b-o', label = '3d averaged')
plt.plot(y_vals, c3_det, 'r-o', label = '2d')
plt.title('$c_3$ v. wall-normal $y$')
plt.legend(loc='best')

plt.figure(5) 
plt.plot(y_vals, epot, 'b-o', label = '3d averaged')
plt.plot(y_vals, epot_det, 'r-o', label = '2d')
plt.title('$\phi$ v. wall-normal $y$')
plt.legend(loc='best')

### show ###
plt.show()


