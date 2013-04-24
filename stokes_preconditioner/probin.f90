! This module stores the runtime parameters.  The probin_init() routine is
! used to initialize the runtime parameters

module probin_module

  use bl_types
  use bl_space

  implicit none

  ! For comments and instructions on how to set the input parameters see namelist section below
  !------------------------------------------------------------- 
  integer,save    :: dim_in,plot_int
  real(dp_t),save :: fixed_dt,theta_fac,p_norm_weight
  integer,save    :: visc_type,bc_lo(MAX_SPACEDIM),bc_hi(MAX_SPACEDIM)
  integer,save    :: seed,mode_coefs(2)
  integer,save    :: n_cells(MAX_SPACEDIM),max_grid_size(MAX_SPACEDIM)
  real(dp_t),save :: prob_lo(MAX_SPACEDIM),prob_hi(MAX_SPACEDIM)
  integer,save    :: prob_coeff, prob_sol, prob_dir, test_type
  real(dp_t),save :: smoothing_width,var_coeff_mag(3),coeff_mag(3), coeff_ratio(3), ABC_coefs(3)

  !------------------------------------------------------------- 
  ! Input parameters controlled via namelist input, with comments
  !------------------------------------------------------------- 

  ! Problem specification
  !----------------------
  namelist /probin/ dim_in        ! 2D or 3D  

  ! Initial conditions and forcing functions (for exact_solutions.f90)
  !----------------------
  namelist /probin/ prob_coeff    ! Initial setup/forcing etc selection 
  namelist /probin/ prob_sol      ! Initial setup/forcing etc selection 
  namelist /probin/ test_type     ! Check gmres, preconditioner, or discretization accuracy
                                  ! (see main_driver.f90 for key)

  namelist /probin/ prob_lo       ! physical lo coordinate
  namelist /probin/ prob_hi       ! physical hi coordinate
  namelist /probin/ n_cells       ! number of cells in domain
  namelist /probin/ max_grid_size ! max number of cells in a box
  namelist /probin/ fixed_dt      ! time step
  namelist /probin/ plot_int      ! Interval for writing a plotfile (for visit/amrvis)

  ! For interfaces, we smooth the concentration with tanh profile to avoid a Gibbs phenomenon:
  namelist /probin/ smoothing_width ! The smoothing width is in units of cell lengths, not physical length!

  ! Wave-indices for incompressible modes inside a box
  ! First component corresponds to prob_dir direction
  ! Second component corresponds to mod(prob_dir+1,DIM)
  namelist /probin/ mode_coefs
  namelist /probin/ prob_dir   ! Which direction is vorticity along
  namelist /probin/ ABC_coefs  ! ABC flow in 3D (a/b/c) (like Taylor vortices in 2D):

  ! Physical parameters
  !----------------------
  ! Every coefficient is of the form
  ! coeff_mag*(deterministic + var_coeff_mag*rand())
  ! The coefficients are indexed from 1-3, index 1 is density, 2 is shear viscosity, and 3 is bulk viscosity
  namelist /probin/ var_coeff_mag ! The magnitude of density/viscosity variation, for random coefficient
  namelist /probin/ coeff_mag     ! The magnitude of density/viscosity
  
  ! For two-phase system, the coefficients for phase 2 are scaled by coeff_ratio
  ! So phase 1 is: coeff_mag*(deterministic + var_coeff_mag*rand())
  ! So phase 2 is: coeff_ratio*coeff_mag*(deterministic + var_coeff_mag*rand())
  namelist /probin/ coeff_ratio ! The ratio for alpha, beta and gamma

  ! Algorithm control / selection
  !----------------------

  ! random number seed
  ! 0        = unpredictable seed based on clock
  ! positive = fixed seed
  namelist /probin/ seed

  ! 0.0 = time-independent
  ! 1.0 = time-dependent
  namelist /probin/ theta_fac

  ! L phi operator
  ! if abs(visc_type) = 1, L = div beta grad
  ! if abs(visc_type) = 2, L = div [ beta (grad + grad^T) ]
  ! if abs(visc_type) = 3, L = div [ beta (grad + grad^T) + I (gamma - (2/3)*beta) div ]
  namelist /probin/ visc_type

  ! Boundary conditions
  !----------------------
  ! BC specifications:
  ! -1 = periodic
  ! 16 = slip    (Dirichlet velocity condition for normal; Dirichlet traction condition for trans)
  ! 17 = no-slip (Dirichlet velocity condition for normal; Dirichlet velocity condition for trans)
  namelist /probin/ bc_lo
  namelist /probin/ bc_hi
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
    
    integer    :: narg, farg

    character(len=128) :: fname
    character(len=128) :: probin_env

    logical :: lexist
    logical :: need_inputs

    integer :: un, ierr

    narg = command_argument_count()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize the runtime parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Defaults

    dim_in = 2

    prob_coeff = 2
    prob_sol = 2
    test_type = 1

    prob_lo(1:MAX_SPACEDIM) = 0.d0
    prob_hi(1:MAX_SPACEDIM) = 1.d0
    n_cells(1:MAX_SPACEDIM) = 64
    max_grid_size(1:MAX_SPACEDIM) = 64
    fixed_dt = 1.d0
    plot_int = 0

    smoothing_width = 1.d0
    mode_coefs(1:2) = 1
    prob_dir = 1
    ABC_coefs(1:3) = 1.d0

    var_coeff_mag(1:3) = 0.d0
    coeff_mag(1:3) = 1.d0
    coeff_ratio(1:3) = 1.d0

    seed = 0
    theta_fac = 1.d0
    visc_type = 1
    p_norm_weight = 1.d0

    bc_lo(1:MAX_SPACEDIM) = PERIODIC
    bc_hi(1:MAX_SPACEDIM) = PERIODIC

    need_inputs = .true.

    farg = 1
    if ( need_inputs .AND. narg >= 1 ) then
       call get_command_argument(farg, value = fname)
       inquire(file = fname, exist = lexist )
       if ( lexist ) then
          farg = farg + 1
          un = unit_new()
          open(unit=un, file = fname, status = 'old', action = 'read')
          read(unit=un, nml = probin)
          close(unit=un)
          need_inputs = .false.
       end if
    end if

    ! stuff that can be read in from the command line by appending, e.g., "--prob_hi_x 64.0"
    do while ( farg <= narg )
       call get_command_argument(farg, value = fname)
       select case (fname)

       case ('--dim_in')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) dim_in

       case ('--prob_coeff')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) prob_coeff

       case ('--prob_sol')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) prob_sol

       case ('--prob_dir')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) prob_dir

       case ('--test_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) test_type

       case ('--smoothing_width')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) smoothing_width

       case ('--var_coeff_mag_alpha')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) var_coeff_mag(1)

       case ('--var_coeff_mag_beta')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) var_coeff_mag(2)

       case ('--var_coeff_mag_gamma')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) var_coeff_mag(3)

       case ('--coeff_mag_alpha')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) coeff_mag(1)

       case ('--coeff_mag_beta')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) coeff_mag(2)

       case ('--coeff_mag_gamma')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) coeff_mag(3)

       case ('--coeff_ratio_alpha')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) coeff_ratio(1)

       case ('--coeff_ratio_beta')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) coeff_ratio(2)

       case ('--coeff_ratio_gamma')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) coeff_ratio(3)

       case ('--mode_coefs_x')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mode_coefs(1)

       case ('--mode_coefs_y')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mode_coefs(2)

       case ('--ABC_coefs_x')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) ABC_coefs(1)

       case ('--ABC_coefs_y')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) ABC_coefs(2)

       case ('--ABC_coefs_z')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) ABC_coefs(3)

       case ('--prob_lo_x')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) prob_lo(1)
       case ('--prob_lo_y')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) prob_lo(2)
       case ('--prob_lo_z')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) prob_lo(3)

       case ('--prob_hi_x')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) prob_hi(1)
       case ('--prob_hi_y')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) prob_hi(2)
       case ('--prob_hi_z')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) prob_hi(3)

       case ('--n_cells_x')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) n_cells(1)
       case ('--n_cells_y')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) n_cells(2)
       case ('--n_cells_z')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) n_cells(3)

       case ('--max_grid_size_x')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) max_grid_size(1)
       case ('--max_grid_size_y')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) max_grid_size(2)
       case ('--max_grid_size_z')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) max_grid_size(3)

       case ('--fixed_dt')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) fixed_dt

       case ('--plot_int')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) plot_int

       case ('--seed')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) seed

       case ('--theta_fac')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) theta_fac

       case ('--visc_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) visc_type

       case ('--p_norm_weight')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) p_norm_weight

       case ('--bc_lo_x')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) bc_lo(1)
       case ('--bc_lo_y')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) bc_lo(2)
       case ('--bc_lo_z')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) bc_lo(3)
       case ('--bc_hi_x')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) bc_hi(1)
       case ('--bc_hi_y')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) bc_hi(2)
       case ('--bc_hi_z')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) bc_hi(3)

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
