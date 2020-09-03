########################################################
##
## This is a README for the induced-charge electro-osmosis 2d input file. 
##
########################################################

There are three important differences between the inputs_iceo file and the 
inputs_electroosmosis scripts:

1) everything is symmetric in the ICEO file--the molar masses and initial concentrations
   are the exact same for the charged species, and they have equal and opposite values for
   the charge per unit mass.


2) the BC are different--for EO, we specify constant neumann conditions for the electric
   potential in the Poisson solve. For ICEO, we want to impose constant Neumann condition 
   on the upper wall, and mixed Dirichlet/Neumann conditions on the lower wall. 
   The constant, homogeneous Neumann condition on the lower wall is enforced by artificially setting
   eps=0 (where eps is the dielectric constant). 

   The Dirichlet BC portion of the lower wall requires a bit of care--since we impose an external 
   electric field that is constant in the x-direction (which implies the potential associated with 
   the external electric field is linear), we have to prescribe a linear Dirichlet BC for the computed
   potential that is equal and oppposite in magnitude to the one due to the external E-field. 

   This is all handled in src_lowMach/inhomogeneous_bc_val.f90, but, one has to make sure that:

   induced_charge_eo = T
  
   Lastly, one can control how large the portion of the lower wall that enforces a Dirichlet condition
   is with the parameters zero_eps_on_wall_left_end and zero_eps_on_wall_right_start. For instance, 
   if zero_eps_on_wall_left_end is set to 0.25, then eps will be set to 0 on the wall from 0*Lx --> 0.25*Lx
   and if zero_eps_on_wall_right_start is set to 0.75, then eps will be set to 0 on the wall from 0.75*Lx --> 1.*Lx. 

3) In ICEO, one can impose an external electric field that varies in time, simply set 

   ac_iceo = T

   The time-dependent E-field looks like a mollified square wave, with amplitude set by E_ext_value. Currently 
   both the frequency and the period of the wave is hard-coded in the function alternating_current_efield,
   located in src_lowMach/inhomogeneous_bc_val.f90. 
