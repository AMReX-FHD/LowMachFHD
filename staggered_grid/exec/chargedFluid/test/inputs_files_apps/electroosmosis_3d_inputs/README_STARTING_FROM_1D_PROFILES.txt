########################################################
##
## This is a README for running a 3d simulation starting from 1d profiles. 
##
########################################################

-Firstly, note: (i) only nspecies = 3 is supported.  
               (ii) only dimension=3 is supported. 
               (iii) the 1d profiles must represent functions of y, not of x or z. 
               (iv) the velocity profile must be for the first (streamwise) component, u, of the
                    the velocity vector. 

-See the comments associated with prob_type=16 in: 
      src_lowMach/init_lowmach.f90. 


-To actually start a simulation from 1d profiles for the species concentration fields
and streamwise velocity, you must set 

      prob_type = 16

in the inputs file.

-The profiles should be stored in .txt files that must be named: 
      c1_1d_vals.txt
      c2_1d_vals.txt
      c3_1d_vals.txt
      u_1d_vals.txt
and should not contain any comments at the top of the file. They should ONLY contain
a single list of numbers corresponding to the values of the field at y-locations, starting
from the first grid cell to the last. 


#################################
##
## Profiles that are stored:
##
#################################

Currently stored in this directory are profiles for the two 3d cases described in the 
other README file in this directory. They are in 
     ./1d_profiles_debye_len_1_10
     ./1d_profiles_debye_len_1_50




