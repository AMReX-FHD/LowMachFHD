! This module stores the runtime parameters.  The probin_init() routine is
! used to initialize the runtime parameters

module probin_module

  use bl_types

  implicit none

  ! For comments and instructions on how to set the input parameters see namelist section below
  !------------------------------------------------------------- 
  logical, save :: use_barodiffusion
  logical, save :: use_bds

  !------------------------------------------------------------- 
  ! Input parameters controlled via namelist input, with comments
  !------------------------------------------------------------- 

  !----------------------
  namelist /probin/ use_barodiffusion  ! enable barodiffusion
  namelist /probin/ use_bds            ! enable bds advection in space and time
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
    use_barodiffusion = .false.
    use_bds = .false.

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

       case ('--use_barodiffusion')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) use_barodiffusion

       case ('--use_bds')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) use_bds

       case ('--')
          farg = farg + 1
          exit

       case default
          if ( .not. parallel_q() ) then
             write(*,*) 'UNKNOWN option = ', fname
             call bl_error("MAIN")
          end if
       end select

       farg = farg + 1
    end do
    
  end subroutine probin_init

end module probin_module
