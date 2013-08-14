! This module stores the runtime parameters.  The probin_lowmach_init() routine is
! used to initialize the runtime parameters

module probin_lowmach_module

  use bl_types

  implicit none

  ! For comments and instructions on how to set the input parameters see 
  ! namelist section below
  !------------------------------------------------------------- 
  integer   , save :: max_step,nscal,print_int
  real(dp_t), save :: rhobar(2),visc_coef,diff_coef,smoothing_width
  real(dp_t), save :: variance_coef,conc_scal,c_init(2),grav(3)
  real(dp_t), save :: mol_mass(2),kT
  integer   , save :: stoch_stress_form,filtering_width
  integer   , save :: project_eos_int
  real(dp_t), save :: material_properties(2,3),c_bc(3,2)

  !------------------------------------------------------------- 
  ! Input parameters controlled via namelist input, with comments
  !------------------------------------------------------------- 

  ! problem setup
  namelist /probin_lowmach/ smoothing_width   ! scale factor for smoothing initial profile
  namelist /probin_lowmach/ c_init            ! controls initial concentration range
  namelist /probin_lowmach/ c_bc              ! c boundary conditions (dir,face)
  namelist /probin_lowmach/ grav              ! gravity vector (negative is downwards)

  ! simulation parameters
  namelist /probin_lowmach/ max_step          ! maximum number of time steps
  namelist /probin_lowmach/ print_int         ! how often to output EOS drift and sum of conserved quantities
  namelist /probin_lowmach/ project_eos_int   ! how often to call project_onto_eos

  ! fluid properties
  namelist /probin_lowmach/ nscal             ! scalars; nscal=2 means we carry rho and rho*c
  namelist /probin_lowmach/ rhobar            ! rho1bar and rho2bar
  namelist /probin_lowmach/ visc_coef         ! momentum diffusion coefficient 'eta'   
  namelist /probin_lowmach/ diff_coef         ! concentration diffusion coefficient 'chi'
  namelist /probin_lowmach/ mol_mass          ! molar mass of species
  namelist /probin_lowmach/ kT                ! temperature
  namelist /probin_lowmach/ material_properties ! a/b for chi/eta/kappa

  ! stochastic properties
  namelist /probin_lowmach/ variance_coef     ! global scaling epsilon for stochastic forcing
  namelist /probin_lowmach/ conc_scal         ! Scaling for concentration stochastic forcing is variance_coeff*conc_scal
  namelist /probin_lowmach/ filtering_width   ! If positive the random numbers will be filtered to smooth out the fields
  namelist /probin_lowmach/ stoch_stress_form ! 0=nonsymmetric (div(v)=0), 1=symmetric (no bulk)

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

    smoothing_width = 1.d0
    c_init(1:2) = 1.d0
    c_bc(1:3,1:2) = 1.d0
    grav(1:3) = 0.d0

    max_step = 1
    print_int = 0
    project_eos_int = 1

    nscal = 2
    rhobar(1) = 1.1d0 
    rhobar(2) = 0.9d0
    visc_coef = 1.d0
    diff_coef = 1.d0
    mol_mass(1:2) = 1.d0
    kT = 1.d0
    material_properties(1:2,1:3) = 0.d0

    variance_coef = 1.d0
    conc_scal = 1.d0
    filtering_width = 0
    stoch_stress_form = 1    

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

       case ('--smoothing_width')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) smoothing_width

       case ('--c_init_1')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_init(1)
       case ('--c_init_2')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_init(2)

       case ('--c_bc_x_lo')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_bc(1,1)
       case ('--c_bc_x_hi')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_bc(1,2)
       case ('--c_bc_y_lo')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_bc(2,1)
       case ('--c_bc_y_hi')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_bc(2,2)
       case ('--c_bc_z_lo')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_bc(3,1)
       case ('--c_bc_z_hi')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_bc(3,2)

       case ('--grav_x')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) grav(1)
       case ('--grav_y')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) grav(2)
       case ('--grav_z')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) grav(3)

       case ('--max_step')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) max_step

       case ('--print_int')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) print_int

       case ('--project_eos_int')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) project_eos_int

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

       case ('--visc_coef')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) visc_coef

       case ('--diff_coef')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) diff_coef

       case ('--mol_mass_1')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mol_mass(1)
       case ('--mol_mass_2')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mol_mass(2)

       case ('--material_properties_1_1')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) material_properties(1,1)
       case ('--material_properties_2_1')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) material_properties(2,1)
       case ('--material_properties_1_2')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) material_properties(1,2)
       case ('--material_properties_2_2')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) material_properties(2,2)
       case ('--material_properties_1_3')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) material_properties(1,3)
       case ('--material_properties_2_3')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) material_properties(2,3)

       case ('--variance_coef')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) variance_coef

       case ('--conc_scal')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) conc_scal

       case ('--filtering_width')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) filtering_width

       case ('--stoch_stress_form')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stoch_stress_form

       case ('--')
          farg = farg + 1
          exit

       case default

       end select

       farg = farg + 1
    end do
    
  end subroutine probin_lowmach_init

end module probin_lowmach_module
