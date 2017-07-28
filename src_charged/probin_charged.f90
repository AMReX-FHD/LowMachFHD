module probin_charged_module

  use bl_types
  use probin_common_module, only: MAX_SPECIES
 
  implicit none

  logical            :: use_charged_fluid
  real(kind=dp_t)    :: dielectric_const
  integer            :: dielectric_type
  real(kind=dp_t)    :: charge_per_mass(MAX_SPECIES)
  real(kind=dp_t)    :: Epot_wall_bc_type
  real(kind=dp_t)    :: Epot_wall(1:2,1:3)
  real(kind=dp_t)    :: theta_pot
  integer            :: num_pot_iters
  real(kind=dp_t)    :: dpdt_factor
  integer            :: E_ext_type
  real(kind=dp_t)    :: E_ext_value(1:3)
  logical            :: electroneutral

  ! for charged fluid
  namelist /probin_charged/ use_charged_fluid
  namelist /probin_charged/ dielectric_const
  namelist /probin_charged/ dielectric_type
  namelist /probin_charged/ charge_per_mass
  namelist /probin_charged/ Epot_wall_bc_type  ! 1 = Dirichlet (fixed potential)
                                               ! 2 = Neumann (fixed charge density)
  namelist /probin_charged/ Epot_wall          ! Dirichlet or Neumann condition
  namelist /probin_charged/ theta_pot          ! for implicit algorithm_type=3, controls
                                               ! temporal discretization for potential term
  namelist /probin_charged/ num_pot_iters
  namelist /probin_charged/ dpdt_factor
  namelist /probin_charged/ E_ext_type         ! external electric field
  namelist /probin_charged/ E_ext_value

  namelist /probin_charged/ electroneutral     ! use electroneutral diffusion fluxes

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
    dielectric_type    = 0      ! 0 = assumes constant epsilon
                                ! 1 = (1+c1)*dielectric_const
                                ! see fluid_charge.f90:compute_permittivity()
    charge_per_mass(:) = 0.d0
    Epot_wall_bc_type  = 1
    Epot_wall(:,:)     = 0.d0
    theta_pot          = 0.5d0
    num_pot_iters      = 2
    dpdt_factor        = 0.d0
    E_ext_type         = 0
    E_ext_value(:)     = 0.d0

    electroneutral = .false.
 
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

       case ('--dielectric_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) dielectric_type

       case ('--Epot_wall_bc_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) Epot_wall_bc_type

       case ('--theta_pot')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) theta_pot

       case ('--num_pot_iters')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) num_pot_iters

       case ('--dpdt_factor')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) dpdt_factor

       case ('--E_ext_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) E_ext_type

       case ('--electroneutral')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) electroneutral

       case ('--')
          farg = farg + 1
          exit
       case default

       end select
       farg = farg + 1
    end do
    
  end subroutine probin_charged_init

end module probin_charged_module
