! This module stores the runtime parameters.  The probin_init() routine is
! used to initialize the runtime parameters

module probin_module

  use bl_types
  use bl_space

  implicit none

  integer, parameter :: N_MAX_SPECIES=2 ! For now the code can only deal with binary mixtures
  
  ! For comments and instructions on how to set the input parameters see namelist section below
  !------------------------------------------------------------- 
  integer,save    :: narg, farg
  integer,save    :: dim_in, nlevs, max_levs
  integer,save    :: max_step, n_steps_skip
  integer,save    :: plot_int, chk_int, regrid_int, proj_int, hydro_grid_int, stats_int
  integer,save    :: print_int
  integer,save    :: amr_buf_width
  integer,save    :: verbose, mg_verbose, cg_verbose
  integer,save    :: mg_max_vcycles,mg_bottom_solver
  real(dp_t),save :: poisson_rel_tol=1.d-9 
  integer,save    :: hg_bottom_solver
  integer,save    :: max_mg_bottom_nlevels
  integer,save    :: stag_mg_verbosity,stag_mg_max_vcycles,stag_mg_maxlevs,stag_mg_minwidth
  integer,save    :: stag_mg_nsmooths_down,stag_mg_nsmooths_up
  real(dp_t),save :: stag_mg_omega
  integer,save    :: restart,prob_type,prob_dir
  real(dp_t),save :: cflfac,init_shrink,fixed_dt
  real(dp_t),save :: visc_coef,diff_coef,bulk_visc
  integer,save    :: visc_type
  real(dp_t),save :: variance_coeff=1,conc_scal=1,initial_variance=0
  real(dp_t),save :: rhobar(N_MAX_SPECIES)
  real(dp_t),save :: vel_init(MAX_SPACEDIM),triangle_coeff
  real(dp_t),save :: stop_time
  real(dp_t),save :: grav
  real(dp_t),save :: smoothing_width
  real(dp_t),save :: wallspeed_lo(MAX_SPACEDIM-1,MAX_SPACEDIM), wallspeed_hi(MAX_SPACEDIM-1,MAX_SPACEDIM)
  integer,save    :: bc_lo(MAX_SPACEDIM), bc_hi(MAX_SPACEDIM)
  integer,save    :: ng_mom, ng_scal
  integer,save    :: n_cells(MAX_SPACEDIM)
  integer,save    :: ref_ratio
  integer,save    :: max_grid_size(MAX_SPACEDIM), max_grid_projection(MAX_SPACEDIM-1)
  integer,save    :: nscal
  integer,save    :: stencil_order
  integer,save    :: temporal_scheme
  logical,save    :: include_concentration=.true.
  logical,save    :: use_stochastic_forcing
  integer,save    :: stoch_forcing_stop
  integer,save    :: stoch_stress_form = 1
  logical, save   :: fake_projection=.false.
  logical,save    :: print_error_norms,print_eos_error,print_conserved,print_cminmax
  logical,save    :: boyce_bc
  logical,save    :: plot_stag
  logical,save    :: enforce_eos
  real(dp_t),save :: prob_lo(MAX_SPACEDIM),prob_hi(MAX_SPACEDIM)
  real(dp_t),save :: c_wall(MAX_SPACEDIM,2)
  real(dp_t),save :: max_dt_growth
  integer,save    :: min_width
  real(dp_t),save :: min_eff
  real(dp_t),save :: c_init(2) 
  real(dp_t),save :: ABC_coefs(3)=1.0d0
  integer, save :: mode_coefs(2)=1
  real(dp_t), save :: kT=1.0d0, mol_mass(N_MAX_SPECIES)=1.0d0
  real(dp_t), save :: material_properties(2,3)=0.0d0 ! a/b for chi/eta/kappa
  integer, save :: filtering_width=0
  integer :: seed=0 
  integer :: n_steps_save_stats=-1
  integer :: project_dir = 0
  logical :: center_snapshots=.false.
  logical :: analyze_conserved=.false.

  ! This will be allocated and defined below
  logical   , allocatable, save :: nodal(:)
  logical   , allocatable, save :: pmask(:)

  integer, parameter :: MAX_ALLOWED_LEVS = 10

  integer :: i

  character(len=128), save :: fixed_grids

  !------------------------------------------------------------- 
  ! Input parameters controlled via namelist input, with comments
  !------------------------------------------------------------- 

  ! Problem type (equations to be solved)
  !----------------------
  namelist /probin/ dim_in ! 2D or 3D  
  namelist /probin/ include_concentration ! This HAS TO BE TRUE for now
  namelist /probin/ use_stochastic_forcing ! Include thermal fluctuations?
  namelist /probin/ stoch_forcing_stop ! if positive, stop the stochastic forcing after this step

  ! 0=nonsymmetric (div(v)=0), 1=symmetric (no bulk)
  ! FIXME: Add option 3=symm. traceless (with bulk)
  namelist /probin/ stoch_stress_form

  ! Delete advection entirely and only do diffusion for both momentum and concentration
  ! This will keep density constant and not enforce the eos -- it is for testing purposes
  namelist /probin/ fake_projection

  ! Physical parameters
  !----------------------
  namelist /probin/ visc_coef ! dynamic viscosity
  namelist /probin/ diff_coef ! mass diffusion
  namelist /probin/ bulk_visc ! bulk viscosity
  namelist /probin/ rhobar ! Pure component densities at global pressure and temperature

  ! In some cases we make chi/eta/kappa quadratic functions of concentration, with coefficients
  ! property = property0*(1+a*c+b*c^2)
  namelist /probin/ material_properties ! a/b for chi/eta/kappa

  namelist /probin/ grav ! Gravitational constant g
  
  ! For stochastic forcing
  namelist /probin/ kT ! k_B*T (temperature) 
  namelist /probin/ mol_mass ! Molecular masses     
  namelist /probin/ variance_coeff ! Global scaling epsilon for stochastic forcing
  namelist /probin/ conc_scal ! Scaling for concentration stochastic forcing is variance_coeff*conc_scal
  namelist /probin/ initial_variance ! Scaling for fluctuations in initial state (if negative, total momentum is set to zero)
  namelist /probin/ filtering_width ! If positive the random numbers will be filtered to smooth out the fields a bit
  
  ! Spatio-temporal discretization
  !----------------------
  namelist /probin/ n_cells ! Grid size

  ! Domain boundaries:
  namelist /probin/ prob_lo
  namelist /probin/ prob_hi
  
  namelist /probin/ fixed_dt ! If positive, fix dt at absolute value

  ! Parallelization parameters:
  namelist /probin/ max_grid_size
  namelist /probin/ max_grid_projection
  
  ! Adaptive time step control (not really tested)
  namelist /probin/ cflfac
  namelist /probin/ init_shrink
  namelist /probin/ max_step
  namelist /probin/ max_dt_growth

  ! AMR is NOT IMPLEMENTED
  namelist /probin/ max_levs
  namelist /probin/ ref_ratio
  namelist /probin/ amr_buf_width
  namelist /probin/ fixed_grids
  namelist /probin/ min_eff
  namelist /probin/ min_width
    
  ! Run parameters
  !----------------------
  namelist /probin/ restart ! New run or restart of previous run?
  namelist /probin/ stop_time ! Do not exceed this time

  namelist /probin/ n_steps_skip ! How many steps to skip

  namelist /probin/ seed ! Positive seed for the RNG (0=use clock)

  ! How often to perform output-related tasks (if negative, never)
  namelist /probin/ print_int ! Print some statistics
  namelist /probin/ plot_int ! Write a plotfile (for visit/amrvis)
  namelist /probin/ chk_int ! Write a checkpoint file
  namelist /probin/ proj_int ! Project umac onto div(umac)=0
  namelist /probin/ stats_int ! Project grid for analysis
  namelist /probin/ regrid_int ! Regridding is NOT IMPLEMENTED

  namelist /probin/ plot_stag ! Write plot file with staggered velocities

  ! Output control:
  namelist /probin/ print_error_norms
  namelist /probin/ print_eos_error
  namelist /probin/ print_conserved
  namelist /probin/ print_cminmax

  ! Verbosity level  
  namelist /probin/ verbose
  namelist /probin/ mg_verbose
  namelist /probin/ cg_verbose
  
  ! Boundary conditions
  !----------------------
  ! BC specifications: -1=periodic,
  ! 11=inlet (reservoir) with Dirichlet for concentration and no-slip,
  ! Or use Neumann conditions for concentration and for velocity
  ! 14=slip wall, 15=noslip wall
  namelist /probin/ bc_lo
  namelist /probin/ bc_hi

  ! For inlet/reservoir boundaries Dirichlet conditions for concentration: 
  namelist /probin/ c_wall ! Shape is (1-DIM=direction, 2=lohi)
  
  ! Each no-slip wall may be moving with a specified tangential velocity along the tangential directions:
  namelist /probin/ wallspeed_lo
  namelist /probin/ wallspeed_hi
  
  namelist /probin/ boyce_bc ! Regularize lid-driven cavity BC for testing (see Boyce Griffith's paper)
  
  ! Analysis of data by projecting along certain directions (see analyze_spectra.f90)
  !----------------------
  namelist /probin/ analyze_conserved ! Should we use conserved variables for the analysis (does not work well)
  namelist /probin/ center_snapshots ! Should we use cell-centered momenta for the analysis (will smooth fluctuations)

  namelist /probin/ hydro_grid_int ! How often to call updateHydroGrid
      ! 0 if never, negative for projectHydroGrid custom analysis, positive for updateHydroGrid
  namelist /probin/ n_steps_save_stats ! How often to dump HydroGrid output files
  
  namelist /probin/ project_dir ! Projection direction (1=x, 2=y, 3=z)
  ! Meaning: 0=analyze 3D data only (no projection needed for HydroGrid, but still need projection if stats_int>0)
  ! +dim=project along dim then analyze 2D only,
  ! -dim=analyze 3D and then project along dim so we also analyze 2D data
  ! It is better to use the conserved variables but it does not quite work for staggered
  
  ! Initial conditions and forcing functions (for exact_solutions.f90)
  !----------------------
  namelist /probin/ prob_type ! Initial setup/forcing etc selection (see exact_solutions.f90 for key)

  namelist /probin/ c_init ! Used for initializing concentration
  namelist /probin/ vel_init ! Used for initializing velocity

  ! For interfaces, we smooth the concentration with tanh profile to avoid a Gibbs phenomenon:
  namelist /probin/ smoothing_width ! The smoothing width is in units of cell lengths, not physical length!

  ! ABC flow in 3D (a/b/c) (like Taylor vortices in 2D):
  namelist /probin/ prob_dir ! Which direction is vorticity along
  namelist /probin/ ABC_coefs
  
  ! Zig-zag (triangle) wave forced concentration with periodic BCs
  namelist /probin/ triangle_coeff
 
  ! Wave-indices for incompressible modes inside a box
  ! First component corresponds to prob_dir direction
  ! Second component corresponds to mod(prob_dir+1,DIM)
  namelist /probin/ mode_coefs
  
  ! Algorithm control / selection
  !----------------------

  ! 0=Euler, 1=Trapezoidal explicit predictor-corrector, 2=Midpoint explicit predictor-corrector,
  ! 3=RK3 with OLD weights, 4=RK3 with NEW weights (preferred!)
  ! If negative, the same random increments are used in all stages
  ! Otherwise the random increments are recalculated each stage
  namelist /probin/ temporal_scheme

  ! L phi operator
  ! if abs(visc_type) = 1, L = div beta grad
  ! if abs(visc_type) = 2, L = div [ beta (grad + grad^T) ]
  ! if abs(visc_type) = 3, L = div [ beta (grad + grad^T) + I (gamma - (2/3)*beta) div ]
  namelist /probin/ visc_type

  namelist /probin/ enforce_eos ! Do an L2 projection onto the EOS after every time step
  
  ! MAC projection solver parameters:
  namelist /probin/ mg_max_vcycles
  namelist /probin/ mg_bottom_solver
  namelist /probin/ poisson_rel_tol
  namelist /probin/ hg_bottom_solver
  namelist /probin/ max_mg_bottom_nlevels

  ! Staggered multigrid solver parameters
  namelist /probin/ stag_mg_verbosity     ! verbosity
  namelist /probin/ stag_mg_max_vcycles   ! max number of v-cycles
  namelist /probin/ stag_mg_maxlevs       ! max number of multigrid levels
  namelist /probin/ stag_mg_minwidth      ! length of box at coarsest multigrid level
  namelist /probin/ stag_mg_nsmooths_down ! number of smooths at each level on the way down
  namelist /probin/ stag_mg_nsmooths_up   ! number of smooths at each level on the way up
  namelist /probin/ stag_mg_omega         ! weighted-jacobi omega coefficient

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
    n_cells=1
    
    grav = 0.d0

    smoothing_width = 1.d0

    rhobar(1) = 1.d0
    rhobar(2) = 1.d0
   
    triangle_coeff = 0.d0

    c_init(1) = 0.25d0
    c_init(2) = 0.75d0

    vel_init = 0.d0

    max_step  = 1
    n_steps_skip = 0
    stop_time = -1.d0

    ref_ratio = 2
    ng_mom  = 1
    ng_scal = 2

    max_levs = 1
    nlevs = -1

    max_grid_size = 256
    max_grid_projection = 128

    stencil_order = 2

    plot_int   = 0
    chk_int    = 0
    regrid_int = -1
    print_int = 1 ! How often to print stats/progress
    proj_int   = -1 ! If positive, how often to project umac onto div(umac)=0
    stats_int = -1 ! If positive, how often to compute mean and standard deviation over reduced dimensions

    hydro_grid_int = 0 ! If non-zero, absolute value telsl us how often to call HydroGrid

    amr_buf_width = -1

    min_eff   = 0.7
    min_width = 1

    prob_lo = 0.d0
    prob_hi = 0.d0

    c_wall = 0.d0

    verbose = 0
    mg_verbose = 0
    cg_verbose = 0
    
    poisson_rel_tol=1.d-10  
    mg_max_vcycles = 50
    mg_bottom_solver = -1
    hg_bottom_solver = -1
    max_mg_bottom_nlevels = 1000

    stag_mg_verbosity = 0
    stag_mg_max_vcycles = 100
    stag_mg_maxlevs = 100
    stag_mg_minwidth = 2
    stag_mg_nsmooths_down = 3
    stag_mg_nsmooths_up = 3
    stag_mg_omega = 1.d0

    init_shrink =  1.0
    fixed_dt    = -1.0

    visc_coef = 0.d0
    diff_coef = 0.d0
    bulk_visc = 0.d0
    material_properties = 0.0d0

    need_inputs = .true.
    fixed_grids = ''
    restart  = -1
  
    prob_type = 2
    prob_dir = 1

    wallspeed_lo = 0.d0
    wallspeed_hi = 0.d0

    bc_lo = PERIODIC
    bc_hi = PERIODIC
  
    max_dt_growth = 1.1d0

    temporal_scheme = 1
    visc_type = 2
    fake_projection = .false.
    use_stochastic_forcing = .false.
    stoch_forcing_stop = -1
    
    print_error_norms = .false.
    print_eos_error = .false.
    print_conserved = .false.
    print_cminmax = .false.
    plot_stag = .false.
    enforce_eos = .false.

    boyce_bc = .false.

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

    do while ( farg <= narg )
       call get_command_argument(farg, value = fname)
       select case (fname)

       case ('--dim_in')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) dim_in

       case ('--stop_time')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stop_time

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

       case ('--max_step')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) max_step

       case ('--n_steps_skip')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) n_steps_skip

       case ('--plot_int')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) plot_int

       case ('--chk_int')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) chk_int

       case ('--regrid_int')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) regrid_int

       case ('--proj_int')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) proj_int

       case ('--print_int')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) print_int

       case ('--amr_buf_width')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) amr_buf_width

       case ('--cflfac')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) cflfac

       case ('--init_shrink')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) init_shrink

       case ('--fixed_dt')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) fixed_dt

       case ('--max_dt_growth')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) max_dt_growth

       case ('--visc_coef')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) visc_coef

       case ('--diff_coef')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) diff_coef

       case ('--bulk_visc')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) bulk_visc

       case ('--visc_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) visc_type

       case ('--variance_coeff')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) variance_coeff

       case ('--conc_scal')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) conc_scal

       case ('--rhobar_1')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) rhobar(1)
       case ('--rhobar_2')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) rhobar(2)

       case ('--triangle_coeff')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) triangle_coeff

       case ('--c_bc_1')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_init(1)
       case ('--c_bc_2')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_init(2)

       case ('--vel_init_x')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) vel_init(1)
       case ('--vel_init_y')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) vel_init(2)
       case ('--vel_init_z')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) vel_init(3)

       case ('--fixed_grids')
          farg = farg + 1
          call get_command_argument(farg, value = fixed_grids)

       case ('--restart')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) restart

       case ('--prob_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) prob_type

       case ('--prob_dir')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) prob_dir

       case ('--wallspeed_lo_yvel_xface')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) wallspeed_lo(1,1)
       case ('--wallspeed_lo_zvel_xface')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) wallspeed_lo(2,1)
       case ('--wallspeed_lo_xvel_yface')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) wallspeed_lo(1,2)
       case ('--wallspeed_lo_zvel_yface')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) wallspeed_lo(2,2)
       case ('--wallspeed_lo_xvel_zface')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) wallspeed_lo(1,3)
       case ('--wallspeed_lo_yvel_zface')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) wallspeed_lo(2,3)

       case ('--wallspeed_hi_yvel_xface')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) wallspeed_hi(1,1)
       case ('--wallspeed_hi_zvel_xface')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) wallspeed_hi(2,1)
       case ('--wallspeed_hi_xvel_yface')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) wallspeed_hi(1,2)
       case ('--wallspeed_hi_zvel_yface')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) wallspeed_hi(2,2)
       case ('--wallspeed_hi_xvel_zface')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) wallspeed_hi(1,3)
       case ('--wallspeed_hi_yvel_zface')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) wallspeed_hi(2,3)

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

       case ('--verbose')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) verbose

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

       case ('--mg_bottom_solver')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mg_bottom_solver

       case ('--poisson_rel_tol')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) poisson_rel_tol

       case ('--hg_bottom_solver')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) hg_bottom_solver

       case ('--max_mg_bottom_nlevels')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) max_mg_bottom_nlevels

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

       case ('--stag_mg_omega')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stag_mg_omega

       case ('--grav')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) grav

       case ('--smoothing_width')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) smoothing_width

       case ('--temporal_scheme')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) temporal_scheme

       case ('--include_concentration')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) include_concentration

       case ('--fake_projection')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) fake_projection

       case ('--use_stochastic_forcing')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) use_stochastic_forcing

       case ('--stoch_forcing_stop')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stoch_forcing_stop

       case ('--stoch_stress_form')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stoch_stress_form

       case ('--hydro_grid_int')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) hydro_grid_int

       case ('--stats_int')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stats_int

       case ('--print_error_norms')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) print_error_norms

       case ('--print_eos_error')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) print_eos_error

       case ('--print_conserved')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) print_conserved

       case ('--plot_stag')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) plot_stag

       case ('--enforce_eos')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) enforce_eos

       case ('--boyce_bc')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) boyce_bc

       case ('--max_levs')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) max_levs

       case ('--n_cellx')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) n_cells(1)
       case ('--n_celly')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) n_cells(2)
       case ('--n_cellz')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) n_cells(3)

       case ('--ref_ratio')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) ref_ratio

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

       case ('--max_grid_projection_1')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) max_grid_projection(1)
       case ('--max_grid_projection_2')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) max_grid_projection(2)

       case ('--min_eff')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) min_eff

       case ('--ABC_coefs')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) ABC_coefs

       case ('--mode_coefs')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mode_coefs

       case ('--min_width')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) min_width

       case ('--seed')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) seed

       case ('--n_steps_save_stats')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) n_steps_save_stats

       case ('--center_snapshots')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) center_snapshots

       case ('--project_dir')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) project_dir

       case ('--analyze_conserved')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) analyze_conserved

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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Error checking
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! At present we only support a single concentration field in the EOS and diffusion implementation
    if(include_concentration) then
      nscal = 2
    else
      call bl_error('at present concentration must be included')
      nscal = 1
    end if

    if (use_stochastic_forcing .and. fixed_dt .eq. -1.0) then
       call bl_error('use_stochastic_forcing=T requires fixed_dt')
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Need to specify regrid_int if max_levs > 1 and not 'fixed grids'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (regrid_int .ne. -1) then
       call bl_error('regrid_int not implemented yet')
    end if

    if (max_levs > 1) then
       if (fixed_grids == '' .and. regrid_int < 1) then
          call bl_error('regrid_int must be specified if max_levs > 1')
       else if (fixed_grids /= '' .and. regrid_int > 0) then
          call bl_warn('Note: regrid_int will be ignored')
       end if
    end if
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Make sure that the buffer width for tagging is at least as big as
    ! the regrid interval
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (regrid_int > 0 .and. amr_buf_width < regrid_int) then
       if (parallel_IOProcessor()) then
          print *,"************************************************************************"
          print *,"WARNING: regrid_int > 0 but amr_buf_width < regrid_int"
          print *,"         setting amr_buf_width = regrid_int"
          print *,"************************************************************************"
       endif
       amr_buf_width = regrid_int
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialize nodal, prob_lo, prob_hi, and pmask
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(nodal(dim_in))
    nodal = .true.

    allocate(pmask(dim_in))
    pmask = .false.
    do i=1,dim_in
       if (bc_lo(i) .eq. PERIODIC .and. bc_hi(i) .eq. PERIODIC) then
          pmask(i) = .true.
       end if
    end do
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialize min_eff
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call cluster_set_min_eff(min_eff)
    call cluster_set_minwidth(min_width)

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

     deallocate(nodal)
     deallocate(pmask)

  end subroutine probin_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module probin_module
