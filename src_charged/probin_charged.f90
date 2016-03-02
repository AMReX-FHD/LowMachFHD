module probin_charged_module

  use bl_types
  use bl_space
 
  implicit none

  integer, parameter :: max_species=10

  logical            :: use_charged_fluid
  real(kind=dp_t)    :: dielectric_const
  real(kind=dp_t)    :: charge_per_mass(max_species)
  real(kind=dp_t)    :: Epot_wall(1:2,1:3)
  real(kind=dp_t)    :: theta_pot
  integer            :: num_pot_iters
  
  ! for charged fluid
  namelist /probin_charged/ use_charged_fluid
  namelist /probin_charged/ dielectric_const
  namelist /probin_charged/ charge_per_mass
  namelist /probin_charged/ Epot_wall
  namelist /probin_charged/ theta_pot          ! for implicit algorithm_type=1, controls
                                               ! temporal discretization for potential term
  namelist /probin_charged/ num_pot_iters

contains

  subroutine probin_charged_init()

    use f2kcli
    use parallel
    use bl_IO_module
    use bl_prof_module
    use bl_error_module
    use bl_constants_module
    use cluster_module
    
    integer            :: narg, farg
    character(len=128) :: fname
    integer            :: un
    logical            :: lexist,need_inputs
    
    narg = command_argument_count()

    ! defaults - will be overwritten by inputs file or command line
    use_charged_fluid  = .false.
    dielectric_const   = 1.d0
    charge_per_mass(:) = 0.d0
    Epot_wall(:,:)     = 0.d0
    theta_pot          = 1.d0
    num_pot_iters      = 2
 
    ! read from input file 
    need_inputs = .true.
    farg = 1
    if ( need_inputs .AND. narg >= 1 ) then
       call get_command_argument(farg, value = fname)
       inquire(file = fname, exist = lexist )
       if ( lexist ) then
          farg = farg + 1
          un = unit_new()
          open(unit=un, file = fname, status = 'old', action = 'read')
          read(unit=un, nml = probin_charged)
          close(unit=un)
          need_inputs = .false.
       end if
    end if

    ! also can be read in from the command line by appending 
    do while ( farg <= narg )
       call get_command_argument(farg, value = fname)
       select case (fname)

       case ('--use_charged_fluid')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) use_charged_fluid

       case ('--dielectric_const')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) dielectric_const

       case ('--theta_pot')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) theta_pot

       case ('--num_pot_iters')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) num_pot_iters

       case ('--')
          farg = farg + 1
          exit
       case default

       end select
       farg = farg + 1
    end do
    
  end subroutine probin_charged_init

end module probin_charged_module
