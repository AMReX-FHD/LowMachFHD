########################################################
##
## This is a README for the electro-osmosis 3d input files. 
##
########################################################

It should be the case that dx=dy=dz for the two input files. 

The main differences between the two electro-osmosis 3d files are the values for (i) c_init, (ii) E_ext, and (iii) Epot_wall

(i) c_init values are changed so as to change the debye length for the different simulations,
    since we fix charge_per_mass and molmass

(ii) E_ext varies so as to maintain a flow velocity that is on the order of 0.1 cm/s, since
     in the deterministic Debye-Huckel theory, the flow velocity is directly proportional to 
     both deybe_len and E_ext. 

(iii) Epot_wall is the value of d(phi)/dn at the wall and must satisfy the compatibility 
      condition for the Poisson equation: -\Delta \phi = q_e/\epsilon 
                                           \implies 
                                           -\int_{\partial \Omega} \partial \phi/ \partial n ds = \int_{\Omega} q/\epsilon dx  (*)

      where \Omega is the domain, and \partial \Omega is the domain boundary. 
      Since \partial \phi/\partial n is constant along the channel walls, (*)  
      let's us determine the value of Epot_wall that makes the compatibility
      condition hold, once we're given the total charge in the simulation.

#####################
## SCRIPT FOR DETERMINING Epot_wall:
#####################

      The script '../electroosmosis_2d_inputs/compute_debye_len_and_epot.f90' prints out the correct value of 
      Epot_wall (and debye length based on initial concentrations) once all 
      other simulation parameters are set. 


With these differences in mind, the inputs files correspond to the following three cases: 

NAME                                  (approximate) DEBYE_LENGTH                 (approximate) ratio: Ly/lambda_d
inputs_electroosmosis1_3d_stochastic                       2.56d-6                                           50
inputs_electroosmosis2_3d_stochastic                       1.28d-5                                           10

where the last column assumes Ly = 1.28d-4. 

#####################
##
## NOTE:
##
#####################
The following parameters are NOT meant to be fixed and may need to be adjusted, depending on the goal of the simulation. 

num_cells
prob_hi
fixed_dt
plot_int, chk_int, print_int
variance_coef_mom, variance_coef_mass
stats_int, n_steps_skip
plot_umac_tavg, plot_Epot_tavg, plot_rho_tavg, plot_avg_gradPhiApprox, plot_mass_fluxes, plot_mass_fluxes_tavg, plot_charge_fluxes, plot_charge_fluxes_tavg








