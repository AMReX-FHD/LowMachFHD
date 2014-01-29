! This module stores the runtime parameters.  The probin_init() routine is
! used to initialize the runtime parameters

module probin_module

  use bl_types

  implicit none

  ! For comments and instructions on how to set the input parameters see namelist section below
  !------------------------------------------------------------- 
  integer, save :: barodiffusion_type
  logical, save :: use_bds,use_overdamped

  !------------------------------------------------------------- 
  ! Input parameters controlled via namelist input, with comments
  !------------------------------------------------------------- 

  !----------------------
  namelist /probin/ barodiffusion_type ! 0 = no barodiffusion
                                       ! 1 = fixed gradp from initialization
                                       ! 2 = update gradp each time step
  namelist /probin/ use_bds            ! enable bds advection in space and time
  namelist /probin/ use_overdamped     ! enable overdamped (infinite Sc) algorithm
  !------------------------------------------------------------- 

contains

  subroutine probin_init()

    use f2kcli
    use parallel
    use bc_module
    use bl_IO_module
    use bl_prof_module
    use bl_error_module
    use bl_constants_module
    use cluster_module
    
    integer    :: narg, farg, un

    character(len=128) :: fname

    logical :: lexist

    narg = command_argument_count()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize the runtime parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Defaults
    barodiffusion_type = 0
    use_bds = .false.
    use_overdamped = .false.

    farg = 1
    if (narg >= 1) then
       call get_command_argument(farg, value = fname)
       inquire(file = fname, exist = lexist )
       if ( lexist ) then
          farg = farg + 1
          un = unit_new()
          open(unit=un, file = fname, status = 'old', action = 'read')
          read(unit=un, nml = probin)
          close(unit=un)
       end if
    end if

    ! stuff that can be read in from the command line by appending, e.g., "--prob_hi_x 64.0"
    do while ( farg <= narg )
       call get_command_argument(farg, value = fname)
       select case (fname)

       case ('--barodiffusion_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) barodiffusion_type

       case ('--use_bds')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) use_bds

       case ('--use_overdamped')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) use_overdamped

       case ('--')
          farg = farg + 1
          exit

       case default

       end select

       farg = farg + 1
    end do
    
  end subroutine probin_init

end module probin_module
