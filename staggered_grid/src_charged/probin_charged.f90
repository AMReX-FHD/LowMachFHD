module probin_charged_module

  use bl_types
  use bl_space
  use probin_common_module, only: MAX_SPECIES
 
  implicit none

  logical            :: use_charged_fluid, print_debye_len
  real(kind=dp_t)    :: dielectric_const
  integer            :: dielectric_type
  real(kind=dp_t)    :: charge_per_mass(MAX_SPECIES)
  real(kind=dp_t)    :: Epot_wall_bc_type(1:2,MAX_SPACEDIM)
  real(kind=dp_t)    :: Epot_wall(1:2,MAX_SPACEDIM)
  real(kind=dp_t)    :: theta_pot
  integer            :: num_pot_iters
  real(kind=dp_t)    :: dpdt_factor
  real(kind=dp_t)    :: relxn_param_charge = 1.0_dp_t
  integer            :: E_ext_type
  real(kind=dp_t)    :: E_ext_value(1:3)
  logical            :: electroneutral
  logical            :: induced_charge_eo, ac_iceo
  integer            :: zero_eps_on_wall_type
  integer            :: zero_charge_on_wall_type
  real(kind=dp_t)    :: zero_eps_on_wall_left_end, zero_eps_on_wall_right_start
  integer            :: epot_mg_verbose
  real(kind=dp_t)    :: epot_mg_abs_tol, epot_mg_rel_tol
  integer            :: bc_function_type
  real(kind=dp_t)    :: L_pos, L_trans, L_zero

  ! for charged fluid
  namelist /probin_charged/ use_charged_fluid
  namelist /probin_charged/ print_debye_len
  namelist /probin_charged/ dielectric_const
  namelist /probin_charged/ dielectric_type
  namelist /probin_charged/ charge_per_mass
  namelist /probin_charged/ Epot_wall_bc_type  ! 1 = Dirichlet (fixed potential)
                                               ! 2 = Neumann (fixed charge density)

  namelist /probin_charged/ bc_function_type   ! 0 = constant 
                                               ! 1 = cubic, see description below

  namelist /probin_charged/ L_pos              ! length of part of boundary where there is positive charge flux, if cubic function is imposed
  namelist /probin_charged/ L_trans            ! length of transition part of boundary, where the value varies like a cubic
  namelist /probin_charged/ L_zero             ! length of part of boundary where there is zero charge flux, if cubic function is imposed

  namelist /probin_charged/ Epot_wall          ! Dirichlet or Neumann condition
  namelist /probin_charged/ theta_pot          ! for implicit algorithm_type=3, controls
                                               ! temporal discretization for potential term
  namelist /probin_charged/ num_pot_iters
  namelist /probin_charged/ dpdt_factor
  namelist /probin_charged/ E_ext_type         ! if 1, sets an external E field to E_ext_value
  namelist /probin_charged/ E_ext_value        ! spacedim-vector specifying external E field

  namelist /probin_charged/ electroneutral     ! use electroneutral diffusion fluxes
  namelist /probin_charged/ induced_charge_eo     ! Are we simulating ICEO?  
  namelist /probin_charged/ ac_iceo               ! If so, do we want to drive the flow with direct or alternating current? 
  namelist /probin_charged/ relxn_param_charge ! Used to prevent slow buildup of charge for electroneutral, keep at 1.0
  namelist /probin_charged/ zero_eps_on_wall_type ! set eps=0 on certain Dirichlet walls
                                                  ! if we want homogeneous Neumann bc's on 
                                                  ! phi for part of a Dirichlet wall
  namelist /probin_charged/ zero_charge_on_wall_type ! set sigma=0 on certain Neumann walls
                                                     ! if we want homogeneous Neumann bc's on 
                                                     ! phi for part of a Dirichlet wall
  namelist /probin_charged/ zero_eps_on_wall_left_end   ! eg if set to 0.25, then eps will be set to 0 on the wall from 0*Lx --> 0.25*Lx
  namelist /probin_charged/ zero_eps_on_wall_right_start     ! eg if set to 0.75, then eps will be set to 0 on the wall from 0.75*Lx --> 1.*Lx
  namelist /probin_charged/ epot_mg_verbose    ! verbosity for poisson solve
  namelist /probin_charged/ epot_mg_abs_tol    ! absolute tolerance for poisson solve
  namelist /probin_charged/ epot_mg_rel_tol    ! relative tolerance for poisson solve



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
    print_debye_len    = .false.
    dielectric_const   = 1.d0
    dielectric_type    = 0      ! 0 = assumes constant epsilon
                                ! 1 = (1+c1)*dielectric_const
                                ! see fluid_charge.f90:compute_permittivity()
    charge_per_mass(:)     = 0.d0
    Epot_wall_bc_type(:,:) = 1
    bc_function_type       = 0 
    L_pos                  = 0.d0
    L_trans                = 0.d0
    L_zero                 = 0.d0
    Epot_wall(:,:)         = 0.d0
    theta_pot              = 0.5d0
    num_pot_iters          = 2
    dpdt_factor            = 0.d0
    E_ext_type             = 0
    E_ext_value(:)         = 0.d0
    epot_mg_verbose        = 0
    epot_mg_abs_tol        = 0.d0
    epot_mg_rel_tol        = 1.d-10

    electroneutral = .false.
    induced_charge_eo = .false.
    ac_iceo = .false. 
    relxn_param_charge = 1.0d0
    zero_eps_on_wall_type = 0
    zero_charge_on_wall_type = 0
    zero_eps_on_wall_left_end = 0.25d0
    zero_eps_on_wall_right_start = 0.75d0
 
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

       case ('--print_debye_len')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) print_debye_len 

       case ('--dielectric_const')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) dielectric_const

       case ('--dielectric_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) dielectric_type

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

       case ('--epot_mg_verbose')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) epot_mg_verbose

       case ('--epot_mg_abs_tol')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) epot_mg_abs_tol

       case ('--epot_mg_rel_tol')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) epot_mg_rel_tol

       case ('--electroneutral')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) electroneutral

       case ('--induced_charge_eo')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) induced_charge_eo

       case ('--ac_iceo')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) ac_iceo 

       case ('--zero_eps_on_wall_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) zero_eps_on_wall_type

       case ('--zero_charge_on_wall_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) zero_charge_on_wall_type

       case ('--')
          farg = farg + 1
          exit
       case default

       end select
       farg = farg + 1
    end do
    
  end subroutine probin_charged_init

end module probin_charged_module
