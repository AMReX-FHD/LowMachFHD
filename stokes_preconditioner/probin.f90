! This module stores the runtime parameters.  The probin_init() routine is
! used to initialize the runtime parameters

module probin_module

  use bl_types
  use bl_space

  implicit none

  ! For comments and instructions on how to set the input parameters see namelist section below
  !------------------------------------------------------------- 
  integer,save    :: dim_in,plot_int,mg_verbose,cg_verbose,mg_max_vcycles,precon_type
  integer,save    :: mg_bottom_solver,mg_nsmooths_down,mg_nsmooths_up,mg_nsmooths_bottom,mg_minwidth
  real(dp_t),save :: mg_rel_tol,gmres_rel_tol, gmres_abs_tol
  integer,save    :: mg_max_bottom_nlevels,stag_mg_verbosity,stag_mg_max_vcycles
  integer,save    :: stag_mg_maxlevs,stag_mg_minwidth,stag_mg_nsmooths_down
  integer,save    :: stag_mg_nsmooths_bottom
  integer,save    :: stag_mg_nsmooths_up,gmres_verbose,gmres_max_outer,gmres_max_inner,gmres_max_iter,gmres_min_iter
  real(dp_t),save :: stag_mg_rel_tol
  real(dp_t),save :: stag_mg_omega,fixed_dt,theta_fac,p_norm_weight
  integer,save    :: stag_mg_smoother,visc_type,bc_lo(MAX_SPACEDIM),bc_hi(MAX_SPACEDIM)
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

  ! preconditioner type
  ! 1 = projection preconditioner
  !-1 = projection preconditioner with expensive pressure update
  ! 2 = lower triangular preconditioner
  !-2 = lower triangular preconditioner with negative sign
  ! 3 = upper triangular preconditioner
  !-3 = upper triangular preconditioner with negative sign
  ! 4 = Block diagonal preconditioner
  !-4 = Block diagonal preconditioner with negative sign
  namelist /probin/ precon_type

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

  ! weighting of pressure when computing norms and inner products
  namelist /probin/ p_norm_weight

  ! Boundary conditions
  !----------------------
  ! BC specifications:
  ! -1 = periodic
  ! 16 = slip    (Dirichlet velocity condition for normal; Dirichlet traction condition for trans)
  ! 17 = no-slip (Dirichlet velocity condition for normal; Dirichlet velocity condition for trans)
  namelist /probin/ bc_lo
  namelist /probin/ bc_hi

  ! MAC projection solver parameters:
  namelist /probin/ mg_verbose            ! multigrid verbosity
  namelist /probin/ cg_verbose            ! BiCGStab (mg_bottom_solver=1) verbosity
  namelist /probin/ mg_max_vcycles        ! maximum number of V-cycles
  namelist /probin/ mg_minwidth           ! length of box at coarsest multigrid level
  namelist /probin/ mg_bottom_solver      ! bottom solver type
                                          ! 0 = smooths only, controlled by mg_nsmooths_bottom
                                          ! 1 = BiCGStab
                                          ! 4 = Fancy bottom solve that coarsens as far as possible 
                                          !     and then applies BiCGStab
  namelist /probin/ mg_nsmooths_down      ! number of smooths at each level on the way down
  namelist /probin/ mg_nsmooths_up        ! number of smooths at each level on the way up
  namelist /probin/ mg_nsmooths_bottom    ! number of smooths at the bottom (only if mg_bottom_solver=0)
  namelist /probin/ mg_max_bottom_nlevels ! for mg_bottom_solver=4, number of additional levels of multigrid
  namelist /probin/ mg_rel_tol            ! relative tolerance stopping criteria

  ! Staggered multigrid solver parameters
  namelist /probin/ stag_mg_verbosity       ! verbosity
  namelist /probin/ stag_mg_max_vcycles     ! max number of v-cycles
  namelist /probin/ stag_mg_maxlevs         ! max number of multigrid levels
  namelist /probin/ stag_mg_minwidth        ! length of box at coarsest multigrid level
  namelist /probin/ stag_mg_nsmooths_down   ! number of smooths at each level on the way down
  namelist /probin/ stag_mg_nsmooths_up     ! number of smooths at each level on the way up
  namelist /probin/ stag_mg_nsmooths_bottom ! number of smooths at the bottom
  namelist /probin/ stag_mg_omega           ! weighted-jacobi omega coefficient
  namelist /probin/ stag_mg_smoother        ! 0 = jacobi; 1 = 2*dm-color Gauss-Seidel
  namelist /probin/ stag_mg_rel_tol         ! relative tolerance stopping criteria

  ! GMRES solver parameters
  namelist /probin/ gmres_rel_tol         ! relative tolerance stopping criteria
  namelist /probin/ gmres_abs_tol         ! absolute tolerance stopping criteria
  namelist /probin/ gmres_verbose         ! gmres verbosity; if greater than 1, more residuals will be printed out
  namelist /probin/ gmres_max_outer       ! max number of outer iterations
  namelist /probin/ gmres_max_inner       ! max number of inner iterations, or restart number
  namelist /probin/ gmres_max_iter        ! max number of gmres iterations
  namelist /probin/ gmres_min_iter        ! min number of gmres iterations
  !------------------------------------------------------------- 

contains

  subroutine probin_init(namelist_file)

    use f2kcli
    use parallel
    use bc_module
    use bl_IO_module
    use bl_prof_module
    use bl_error_module
    use bl_constants_module
    use cluster_module
    integer, intent(out), optional :: namelist_file
    
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

    precon_type = 1
    seed = 0
    theta_fac = 1.d0
    visc_type = 1
    p_norm_weight = 1.d0

    bc_lo(1:MAX_SPACEDIM) = PERIODIC
    bc_hi(1:MAX_SPACEDIM) = PERIODIC

    mg_verbose = 0
    cg_verbose = 0
    mg_max_vcycles = 100
    mg_minwidth = 2
    mg_bottom_solver = 0
    mg_nsmooths_down = 2
    mg_nsmooths_up = 2
    mg_nsmooths_bottom = 8
    mg_max_bottom_nlevels = 1000
    mg_rel_tol = 1.d-9

    stag_mg_verbosity = 0
    stag_mg_max_vcycles = 100
    stag_mg_maxlevs = 100
    stag_mg_minwidth = 2
    stag_mg_nsmooths_down = 2
    stag_mg_nsmooths_up = 2
    stag_mg_nsmooths_bottom = 8
    stag_mg_omega = 1.d0
    stag_mg_smoother = 1
    stag_mg_rel_tol = 1.d-9

    gmres_rel_tol = 1.d-9
    gmres_abs_tol = 0.d0
    gmres_verbose = 1
    gmres_max_outer = 20
    gmres_max_inner = 5
    gmres_max_iter = 100
    gmres_min_iter = 1

    need_inputs = .true.
    
    call get_environment_variable('PROBIN', probin_env, status = ierr)
    if ( need_inputs .AND. ierr == 0 ) then
       un = unit_new()
       open(unit=un, file = probin_env, status = 'old', action = 'read')
       call read_nml()
       need_inputs = .false.
    end if

    farg = 1
    if ( need_inputs .AND. narg >= 1 ) then
       call get_command_argument(farg, value = fname)
       inquire(file = fname, exist = lexist )
       if ( lexist ) then
          farg = farg + 1
          un = unit_new()
          open(unit=un, file = fname, status = 'old', action = 'read')
          call read_nml()
          need_inputs = .false.
       end if
    end if

    inquire(file = 'inputs_varden', exist = lexist)
    if ( need_inputs .AND. lexist ) then
       un = unit_new()
       open(unit=un, file = 'inputs_varden', status = 'old', action = 'read')
       call read_nml()
       need_inputs = .false.
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

       case ('--precon_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) precon_type

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

       case ('--mg_verbose')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mg_verbose

       case ('--cg_verbose')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) cg_verbose

       case ('--mg_max_vcycles')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mg_max_vcycles

       case ('--mg_minwidth')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mg_minwidth

       case ('--mg_bottom_solver')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mg_bottom_solver

       case ('--mg_nsmooths_down')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mg_nsmooths_down

       case ('--mg_nsmooths_up')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mg_nsmooths_up

       case ('--mg_nsmooths_bottom')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mg_nsmooths_bottom

       case ('--mg_max_bottom_nlevels')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mg_max_bottom_nlevels

       case ('--mg_rel_tol')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mg_rel_tol

       case ('--stag_mg_verbosity')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stag_mg_verbosity

       case ('--stag_mg_max_vcycles')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stag_mg_max_vcycles

       case ('--stag_mg_maxlevs')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stag_mg_maxlevs

       case ('--stag_mg_minwidth')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stag_mg_minwidth

       case ('--stag_mg_nsmooths_down')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stag_mg_nsmooths_down

       case ('--stag_mg_nsmooths_up')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stag_mg_nsmooths_up

       case ('--stag_mg_nsmooths_bottom')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stag_mg_nsmooths_bottom

       case ('--stag_mg_omega')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stag_mg_omega

       case ('--stag_mg_smoother')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stag_mg_smoother

       case ('--gmres_rel_tol')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) gmres_rel_tol

       case ('--gmres_abs_tol')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) gmres_rel_tol

       case ('--gmres_verbose')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) gmres_verbose

       case ('--gmres_max_outer')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) gmres_max_outer

       case ('--gmres_max_inner')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) gmres_max_inner
      
       case ('--gmres_max_iter')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) gmres_max_iter

      case ('--gmres_min_iter')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) gmres_min_iter

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

  contains
  
    subroutine read_nml()
       read(unit=un, nml = probin)       
       if(present(namelist_file)) then
         namelist_file=un ! We need to read more stuff, it will be closed by someone else
       else
         close(unit=un)
       end if  
    end subroutine
    
  end subroutine probin_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine probin_close()

  end subroutine probin_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module probin_module
