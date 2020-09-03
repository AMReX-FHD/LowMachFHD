!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!! Given a bunch of parameters from an input file,
!! this script computes the Debye length of the simulation, based
!! on the initial concentrations. 
!!
!! The script also computes the value of Epot_wall necessary for the 
!! compatibility condition to be satisfied (if we prescribe a Neumann 
!! boundary condition for the electric potential at the solid walls). 
!!
!! IMPORTANT TO NOTE:
!!       The parameter 'code_tot_q' is crucial to get correct, in order to 
!!       ensure Epot_wall actually satisfies the compatibility condition. 
!!       The best way to get the correct value is to run the simulation once
!!       with the prescribed concentrations and charge_per_mass values. Although
!!       the gmres solver will crash of the bat, the code prints off 'total charge' 
!!       and this is the value one should use for code_tot_q. 
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program compute_debye_len_and_e_pot


double precision,  parameter:: rho = 1.000000d0       !Density 
double precision, parameter::c1 = 4.2629d-6            !concentration 1
double precision, parameter::c2 = 3.4299d-6           !concentration 2
double precision, parameter::q1 = 4.2d3                  !charge per unit mass 1
double precision, parameter::q2 = -2.72d3                 !charge per unit mass 2
real(kind=8), parameter::m1 = 3.82d-23                 !molar mass 1
real(kind=8), parameter::m2 = 5.89d-23                 !molar mass 2
real(kind=8), parameter::m3 = 3.35d-23                 !molar mass 3

real(kind=8), parameter::T = 300.d0                      !Temperature
real(kind=8), parameter::k_B = 1.3806488d-16           !Boltzmann's constant
real(kind=8), parameter::eps  = 6.91d-19               !dielectric constant

real(kind=8):: c3, lambda, e_stat_stab                 !concentration 3 and debye len and electrostatic stability parameter
real(kind=8), parameter::Lx = 5.12d-4                  !domain length in x 
real(kind=8), parameter::Ly = 1.28d-4                  !domain length in y 
real(kind=8), parameter::Lz = 1.28d-4                  !domain length in z
real(kind=8), parameter::code_tot_q  = 7.1931450880846915d-014 ! need this from an initial, preliminary run
                                                               ! to compute E_pot to 16 decimal places... :-(
real(kind=8), dimension(1:3) ::  Dbar
real(kind=8), parameter::Dbar1 = 1.17d-5
real(kind=8), parameter::Dbar2 = 1.33d-5
real(kind=8), parameter::Dbar3 = 2.03d-5

c3 = 1.d0 - c1 - c2
Dbar(1) = Dbar1
Dbar(2) = Dbar2
Dbar(3) = Dbar3

print*, 'total charge should roughly be      : ', rho*(q1*c1 + q2*c2)*Lx*Ly*Lz
print*, 'code total charge is                : ', code_tot_q

lambda = sqrt(eps*k_B*T/(rho*(c1*m1*q1*q1 + c2*m2*q2*q2)))  ! Debye length
print*, 'Debye length is                     : ', lambda
print*, 'Ratio of Ly to Debye length is      : ', Ly/lambda

e_stat_stab = lambda**2.d0/maxval(Dbar)

print*, 'Electrostatic stability requires dt < ', e_stat_stab

print*, 'E_pot needs to be                   : ', code_tot_q/eps/2.d0/Lx/Lz

end program 
