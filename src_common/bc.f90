module bc_module

  implicit none

  integer, parameter, public :: BC_PER       = -1
  integer, parameter, public :: BC_INT       = 0
  integer, parameter, public :: BC_DIR       = 1
  integer, parameter, public :: BC_NEU       = 2

  integer, parameter, public :: PERIODIC     = -1
  integer, parameter, public :: INTERIOR     =  0

  integer, parameter, public :: INLET        = 11
  integer, parameter, public :: OUTLET       = 12
  integer, parameter, public :: SYMMETRY     = 13

  integer, parameter, public ::         SLIP_WALL = 16
  integer, parameter, public ::      NO_SLIP_WALL = 17
  integer, parameter, public ::    SLIP_RESERVOIR = 18
  integer, parameter, public :: NO_SLIP_RESERVOIR = 19

  integer, parameter, public :: EXT_DIR      =  23
  integer, parameter, public :: FOEXTRAP     =  24
  integer, parameter, public :: HOEXTRAP     =  25
  integer, parameter, public :: DIR_VEL      =  26
  integer, parameter, public :: DIR_TRACT    =  27

end module bc_module
