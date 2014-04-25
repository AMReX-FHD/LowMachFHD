! This module stores the runtime parameters.  The probin_binarylm_init() routine is
! used to initialize the runtime parameters

module probin_binarylm_module

  use bl_types

  implicit none

  ! For comments and instructions on how to set the input parameters see 
  ! namelist section below
  !------------------------------------------------------------- 
  integer   , save :: restart,max_step,print_int
  real(dp_t), save :: rhobar(2),diff_coef,smoothing_width
  real(dp_t), save :: initial_variance,conc_scal,c_init(2),u_init(2),grav(3)
  real(dp_t), save :: mol_mass(2),temperature
  integer   , save :: project_eos_int
  real(dp_t), save :: material_properties(3,3),c_bc(3,2)
  integer   , save :: algorithm_type,barodiffusion_type
  logical, save :: analyze_binary

  !------------------------------------------------------------- 
  ! Input parameters controlled via namelist input, with comments
  !------------------------------------------------------------- 

  ! problem setup
  namelist /probin_binarylm/ smoothing_width   ! scale factor for smoothing initial profile
  namelist /probin_binarylm/ c_init            ! controls initial concentration range
  namelist /probin_binarylm/ u_init            ! controls initial velocity
  namelist /probin_binarylm/ c_bc              ! c boundary conditions (dir,face).
                                               ! Dirichlet for RESERVOIR; Neumann for WALL
  namelist /probin_binarylm/ grav              ! gravity vector (negative is downwards)

  ! simulation parameters
  namelist /probin_binarylm/ restart           ! checkpoint restart number
  namelist /probin_binarylm/ max_step          ! maximum number of time steps
  namelist /probin_binarylm/ print_int         ! how often to output EOS drift and sum of conserved quantities
  namelist /probin_binarylm/ project_eos_int   ! how often to call project_onto_eos

  ! fluid properties
  namelist /probin_binarylm/ rhobar            ! rho1bar and rho2bar
  namelist /probin_binarylm/ diff_coef         ! concentration diffusion coefficient 'chi'
  namelist /probin_binarylm/ mol_mass          ! molar mass of species
  namelist /probin_binarylm/ temperature       ! temperature
  namelist /probin_binarylm/ material_properties ! Coefficients A/B/C for chi/eta/kappa
     ! Formula is: indx=1 for chi, indx=2 for eta, indx=3 for kappa (NOT implemented yet)
     ! coeff=coeff0*(material_properties(1,indx) + material_properties(2,indx)*c) / (1.d0 + material_properties(3,indx)*c)

  ! stochastic properties
  namelist /probin_binarylm/ initial_variance  ! multiplicative factor for initial fluctuations
                                               ! (if negative, total momentum is set to zero)
  namelist /probin_binarylm/ conc_scal         ! Scaling for concentration stochastic forcing is variance_coeff*conc_scal

  namelist /probin_binarylm/ barodiffusion_type ! 0 = no barodiffusion
                                                ! 1 = fixed gradp from initialization
                                                ! 2 = update gradp each time step
  namelist /probin_binarylm/ algorithm_type     ! 0 = Inertial algorithm
                                                ! 1 = Overdamped with 1 RNG
                                                ! 2 = Overdamped with 2 RNGs
  namelist /probin_binarylm/ analyze_binary
                             ! Call the older analyze_spectra_binary or the new analyze_spectra?  

contains

  subroutine probin_binarylm_init()

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
    u_init(1:2) = 0.d0
    c_bc(1:3,1:2) = 0.d0
    grav(1:3) = 0.d0

    restart = -1
    max_step = 1
    print_int = 0
    project_eos_int = 1

    rhobar(1) = 1.1d0 
    rhobar(2) = 0.9d0
    diff_coef = 1.d0
    mol_mass(1:2) = 1.d0
    temperature = 1.d0
    material_properties(1:3,1:3) = 0.d0

    initial_variance = 0.d0
    conc_scal = 1.d0

    barodiffusion_type = 0
    algorithm_type = 0
    
    analyze_binary=.true.

    farg = 1
    if (narg >= 1) then
       call get_command_argument(farg, value = fname)
       inquire(file = fname, exist = lexist )
       if ( lexist ) then
          farg = farg + 1
          un = unit_new()
          open(unit=un, file = fname, status = 'old', action = 'read')
          read(unit=un, nml = probin_binarylm)
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

       case ('--u_init_1')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) u_init(1)
       case ('--u_init_2')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) u_init(2)

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

       case ('--restart')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) restart

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

       case ('--rhobar_1')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) rhobar(1)
       case ('--rhobar_2')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) rhobar(2)

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
       case ('--material_properties_3_1')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) material_properties(3,1)

       case ('--material_properties_1_2')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) material_properties(1,2)
       case ('--material_properties_2_2')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) material_properties(2,2)
       case ('--material_properties_3_2')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) material_properties(3,2)

       case ('--material_properties_1_3')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) material_properties(1,3)
       case ('--material_properties_2_3')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) material_properties(2,3)
       case ('--material_properties_3_3')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) material_properties(3,3)

       case ('--initial_variance')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) initial_variance

       case ('--conc_scal')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) conc_scal

       case ('--barodiffusion_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) barodiffusion_type

       case ('--algorithm_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) algorithm_type

       case ('--analyze_binary')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) analyze_binary

       case ('--')
          farg = farg + 1
          exit

       case default

       end select

       farg = farg + 1
    end do
    
  end subroutine probin_binarylm_init

end module probin_binarylm_module
