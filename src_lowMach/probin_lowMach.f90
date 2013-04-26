! This module stores the runtime parameters.  The probin_lowmach_init() routine is
! used to initialize the runtime parameters

module probin_lowmach_module

  use bl_types

  implicit none

  ! For comments and instructions on how to set the input parameters see 
  ! namelist section below
  !------------------------------------------------------------- 
  integer   , save :: max_step,nscal
  real(dp_t), save :: rhobar(2)

  !------------------------------------------------------------- 
  ! Input parameters controlled via namelist input, with comments
  !------------------------------------------------------------- 

  !----------------------
  namelist /probin_lowmach/ max_step      ! maximum number of time steps
  namelist /probin_lowmach/ nscal         ! scalars; nscal=2 means we carry rho and rho*c

contains

  subroutine probin_lowmach_init()

    use f2kcli
    use parallel
    use bc_module
    use bl_IO_module
    use bl_prof_module
    use bl_error_module
    use bl_constants_module
    use cluster_module
    
    integer    :: narg, farg

    character(len=128) :: fname

    logical :: lexist

    integer :: un

    narg = command_argument_count()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize the runtime parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Defaults

    max_step = 1
    nscal = 2

    rhobar(1) = 1.d0
    rhobar(2) = 1.d0

    farg = 1
    if (narg >= 1) then
       call get_command_argument(farg, value = fname)
       inquire(file = fname, exist = lexist )
       if ( lexist ) then
          farg = farg + 1
          un = unit_new()
          open(unit=un, file = fname, status = 'old', action = 'read')
          read(unit=un, nml = probin_lowmach)
          close(unit=un)
       end if
    end if

    ! stuff that can be read in from the command line by appending, e.g., "--prob_hi_x 64.0"
    do while ( farg <= narg )
       call get_command_argument(farg, value = fname)
       select case (fname)

       case ('--max_step')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) max_step

       case ('--nscal')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) nscal

       case ('--rhobar_1')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) rhobar(1)
       case ('--rhobar_2')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) rhobar(2)

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
    
  end subroutine probin_lowmach_init

end module probin_lowmach_module
