import numpy as np
import matplotlib.pyplot as plt

debye_len_case = '1_50'
Ly = 1.28e-4
#debye_len_case = '1_10'

if (debye_len_case == '1_50'):
	E_ext = 3.e9
	dphi_dn_wall = 794201762623.94604
elif (debye_len_case == '1_10'): 
	E_ext = 3.e10
	dphi_dn_wall = 31767657979.612209

# import deterministic potential and velocity

u_det = np.loadtxt('u_deterministic_' + debye_len_case + '.txt', unpack = True, usecols = [0])
epot = np.loadtxt('epot_deterministic_' + debye_len_case + '.txt', unpack = True, usecols = [0])
c1 = np.loadtxt('c1_deterministic_' + debye_len_case + '.txt', unpack = True, usecols = [0])
c2 = np.loadtxt('c2_deterministic_' + debye_len_case + '.txt', unpack = True, usecols = [0])
c3 = np.loadtxt('c3_deterministic_' + debye_len_case + '.txt', unpack = True, usecols = [0])


#####################
##  need to calculuate lambda_d based on centerline concentration values: 
#####################
mu = 1.05e-2                           # momentum diffusion coefficient
dielectric_const = 6.91e-19            # dielectric constant
k_B = 1.3806488e-16
T = 300.

q1 = 4.2e3
q2 = -2.72e3 

m1 = 3.82e-23
m2 = 5.89e-23
m3 = 3.35e-23
rho = 1.
c1_c = c1[len(c1)/2]
c2_c = c2[len(c2)/2]
c3_c = c3[len(c3)/2]
epot_c = epot[len(epot)/2]
u_c = u_det[len(u_det)/2]

######################
## check concentration profile against Boltzmann distribution
######################
c1_boltz = c1_c*np.exp(-m1*q1/k_B/T*(epot-epot_c))
c1_approx = c1_c*(1. - m1*q1/k_B/T*(epot-epot_c))

c2_boltz = c2_c*np.exp(-m2*q2/k_B/T*(epot-epot_c))
c2_approx = c2_c*(1. - m2*q2/k_B/T*(epot-epot_c))


#####################
##  check relationship: 
##
##  u = \epsilon E_x/ \mu (\phi - \phi_{wall})
## 
#####################
shift = epot_c - u_c*mu/dielectric_const/E_ext
epot_wall = epot[0]
u_theory = dielectric_const*E_ext/mu * (epot - shift)

######################
## compute debye length based on concentrations at channel center
######################
lambda_d = np.sqrt(dielectric_const*k_B*T/(rho*(c1_c*m1*q1**2 + c2_c*m2*q2**2)))
print 'debye len is: ', lambda_d
print 'ratio Ly/lambda_d: ', Ly/lambda_d




######################
## compare centerline velocity with D-H theory
######################
print 'Debye Huckel slip velocity: ', lambda_d*dielectric_const*dphi_dn_wall*E_ext/mu
print 'Deterministic calculation vel: ', u_c



plt.figure(1)
plt.plot(u_det, 'b-o', label = 'simulation $u$')
plt.plot(u_theory, 'r-o', label = '$\epsilon E_x (\phi-\phi_w)/\mu$')
plt.title('Velocity')
plt.legend(loc = 'best')

plt.figure(5)
plt.plot(u_det-u_theory, 'g-o')
plt.title('Difference between velocities')

plt.figure(2)
plt.plot(c1, 'g-o', label = 'simulation')
plt.plot(c1_boltz, 'r-o', label = 'theory') 
plt.plot(c1_approx, 'b-o', label = 'approximation') 
plt.title('Concentration 1')
plt.legend(loc='best')

plt.figure(3)
plt.plot(c2, 'g-o', label = 'simulation')
plt.plot(c2_boltz, 'r-o', label = 'theory') 
plt.plot(c2_approx, 'b-o', label = 'approximation') 
plt.title('Concentration 2')
plt.legend(loc='best')

plt.show()
