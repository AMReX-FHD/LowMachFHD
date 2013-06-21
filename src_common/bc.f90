module bc_module

  implicit none

  integer, save :: vel_bc_comp, pres_bc_comp, scal_bc_comp, num_scal_bc
  
  ! These are physical boundary condition types
  ! We set these in the inputs file and they get translated
  ! by define_bc_tower.f90 into the definitions below
  integer, parameter, public ::            INLET  = 11 ! not currently used but needed by BoxLib
  integer, parameter, public ::            OUTLET = 12 ! not currently used but needed by BoxLib
  integer, parameter, public ::          SYMMETRY = 13 ! not currently used but needed by BoxLib
  integer, parameter, public ::         SLIP_WALL = 16
  integer, parameter, public ::      NO_SLIP_WALL = 17
  integer, parameter, public ::    SLIP_RESERVOIR = 18
  integer, parameter, public :: NO_SLIP_RESERVOIR = 19

  ! These specify boundary conditions on phi for the cc mg solver  
  integer, parameter, public :: BC_PER       = -1
  integer, parameter, public :: BC_INT       = 0
  integer, parameter, public :: BC_DIR       = 1 ! not currently used but needed by BoxLib
  integer, parameter, public :: BC_NEU       = 2

  ! These specify boundary conditions for the velocity, pressure, and scalars
  integer, parameter, public :: PERIODIC     = -1
  integer, parameter, public :: INTERIOR     =  0
  integer, parameter, public :: EXT_DIR      =  23
  integer, parameter, public :: FOEXTRAP     =  24
  integer, parameter, public :: HOEXTRAP     =  25
  integer, parameter, public :: DIR_VEL      =  26
  integer, parameter, public :: DIR_TRACT    =  27

end module bc_module
