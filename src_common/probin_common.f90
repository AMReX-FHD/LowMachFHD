! This module stores the runtime parameters.  The probin_init() routine is
! used to initialize the runtime parameters

module probin_common_module

  use bl_types
  use bl_space

  implicit none

  ! For comments and instructions on how to set the input parameters see namelist section below
  !------------------------------------------------------------- 
  integer,save    :: dim_in,plot_int,chk_int,prob_type,advection_type
  real(dp_t),save :: fixed_dt,cfl,grav(3)
  integer,save    :: visc_type,diff_type,bc_lo(MAX_SPACEDIM),bc_hi(MAX_SPACEDIM),seed
  integer,save    :: n_cells(MAX_SPACEDIM),max_grid_size(MAX_SPACEDIM)
  real(dp_t),save :: prob_lo(MAX_SPACEDIM),prob_hi(MAX_SPACEDIM)
  real(dp_t),save :: wallspeed_lo(MAX_SPACEDIM-1,MAX_SPACEDIM)
  real(dp_t),save :: wallspeed_hi(MAX_SPACEDIM-1,MAX_SPACEDIM)
  integer,save    :: hydro_grid_int,project_dir,max_grid_projection(2)
  integer,save    :: stats_int,n_steps_save_stats,n_steps_skip
  logical,save    :: analyze_conserved,center_snapshots
  real(dp_t),save :: variance_coef,k_B,visc_coef
  integer   , save :: stoch_stress_form,filtering_width
  real(dp_t), save :: boussinesq_beta

  !------------------------------------------------------------- 
  ! Input parameters controlled via namelist input, with comments
  !------------------------------------------------------------- 

  ! Problem specification
  !----------------------
  namelist /probin_common/ dim_in        ! 2D or 3D  
  namelist /probin_common/ prob_lo       ! physical lo coordinate
  namelist /probin_common/ prob_hi       ! physical hi coordinate
  namelist /probin_common/ n_cells       ! number of cells in domain
  namelist /probin_common/ max_grid_size ! max number of cells in a box
  namelist /probin_common/ fixed_dt      ! time step
  namelist /probin_common/ cfl           ! cfl number
  namelist /probin_common/ plot_int      ! Interval for writing a plotfile (for visit/amrvis)
  namelist /probin_common/ chk_int       ! Interval for writing a checkpoint
  namelist /probin_common/ prob_type     ! sets scalars, m, coefficients (see init.f90)

  ! Algorithm control / selection
  !----------------------

  ! random number seed
  ! 0        = unpredictable seed based on clock
  ! positive = fixed seed
  namelist /probin_common/ seed

  namelist /probin_common/ grav               ! gravity vector (negative is downwards)
  namelist /probin_common/ boussinesq_beta    ! beta for boussinesq gravity

  ! L phi operator
  ! if abs(visc_type) = 1, L = div beta grad
  ! if abs(visc_type) = 2, L = div [ beta (grad + grad^T) ]
  ! if abs(visc_type) = 3, L = div [ beta (grad + grad^T) + I (gamma - (2/3)*beta) div ]
  ! positive = assume constant coefficients
  ! negative = assume spatially-varying coefficients
  namelist /probin_common/ visc_type
  namelist /probin_common/ visc_coef         ! momentum diffusion coefficient 'eta'   

  ! 1 = constant coefficients
  ! -1 = spatially-varing coefficients
  namelist /probin_common/ diff_type

  namelist /probin_common/ advection_type ! 0 = centered explicit
                                          ! 1 = unlimited bilinear bds in space and time
                                          ! 2 = limited bliniear bds in space and time
                                          ! 3 = unlimited quadratic bds in space and time

  ! Boundary conditions
  !----------------------
  ! BC specifications:
  ! -1 = periodic
  ! 100 = no-slip wall      (Dir condition for normal vel; Dir velocity condition for trans vel)
  ! 101 = no-slip reservoir (Dir condition for normal vel; Dir velocity condition for trans vel)
  ! 200 = slip wall         (Dir condition for normal vel; Dir traction condition for trans vel)
  ! 201 = slip reservoir    (Dir condition for normal vel; Dir traction condition for trans vel)
  ! For a complete list see bc.f90
  namelist /probin_common/ bc_lo
  namelist /probin_common/ bc_hi

  ! Each no-slip wall may be moving with a specified tangential 
  ! velocity along the tangential directions
  ! In 2D:
  ! wallspeed_lo/hi(1,1) - yvel on x-face
  ! wallspeed_lo/hi(1,2) - xvel on y-face
  ! In 3D:
  ! wallspeed_lo/hi(1,1) - yvel on x-face
  ! wallspeed_lo/hi(2,1) - zvel on x-face
  ! wallspeed_lo/hi(1,2) - xvel on y-face
  ! wallspeed_lo/hi(2,2) - zvel on y-face
  ! wallspeed_lo/hi(1,3) - xvel on z-face
  ! wallspeed_lo/hi(2,3) - yvel on z-face
  namelist /probin_common/ wallspeed_lo
  namelist /probin_common/ wallspeed_hi

  ! Control for analyze_spectra.90
  !----------------------

  namelist /probin_common/ hydro_grid_int     ! How often to call updateHydroGrid
                                               ! 0 if never
                                               ! negative for projectHydroGrid custom analysis
                                               ! positive for updateHydroGrid

  namelist /probin_common/ project_dir         ! Projection direction (1=x, 2=y, 3=z)
  ! Meaning: 0=analyze 3D data only (no projection needed for HydroGrid, 
  !          but still need projection if stats_int>0)
  ! +dim=project along dim then analyze 2D only,
  ! -dim=analyze 3D and then project along dim so we also analyze 2D data
  ! It is better to use the conserved variables but it does not quite work for staggered

  namelist /probin_common/ max_grid_projection ! parallelization parameters
  namelist /probin_common/ stats_int           ! Project grid for analysis
                                                ! If positive, how often to compute mean and 
                                                ! standard deviation over reduced dimensions
  namelist /probin_common/ n_steps_save_stats  ! How often to dump HydroGrid output files
  namelist /probin_common/ n_steps_skip        ! How many steps to skip
  namelist /probin_common/ analyze_conserved   ! Should we use conserved variables for the analysis
                                                ! (does not work well)
  namelist /probin_common/ center_snapshots    ! Should we use cell-centered momenta for the analysis
                                                ! (will smooth fluctuations)

  ! stochastic properties
  namelist /probin_common/ variance_coef     ! global scaling epsilon for stochastic forcing
  namelist /probin_common/ k_B               ! Boltzmann's constant
  namelist /probin_common/ filtering_width   ! If positive the random numbers will be filtered to smooth out the fields
  namelist /probin_common/ stoch_stress_form ! 0=nonsymmetric (div(v)=0), 1=symmetric (no bulk)

  !------------------------------------------------------------- 

contains

  subroutine probin_common_init()

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
    logical :: need_inputs

    narg = command_argument_count()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize the runtime parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Defaults

    dim_in = 2

    prob_lo(1:MAX_SPACEDIM) = 0.d0
    prob_hi(1:MAX_SPACEDIM) = 1.d0
    n_cells(1:MAX_SPACEDIM) = 64
    max_grid_size(1:MAX_SPACEDIM) = 64
    fixed_dt = 1.d0
    cfl = 0.5d0
    plot_int = 0
    chk_int = 0
    prob_type = 1

    seed = 1

    grav(1:3) = 0.d0
    boussinesq_beta = 0.d0

    visc_type = 1
    visc_coef = 1.d0
    diff_type = 1

    advection_type = 0

    bc_lo(1:MAX_SPACEDIM) = PERIODIC
    bc_hi(1:MAX_SPACEDIM) = PERIODIC

    wallspeed_lo(1:MAX_SPACEDIM-1,1:MAX_SPACEDIM) = 0.d0
    wallspeed_hi(1:MAX_SPACEDIM-1,1:MAX_SPACEDIM) = 0.d0   

    hydro_grid_int = 0
    project_dir = 0
    max_grid_projection = 128
    stats_int = -1
    n_steps_save_stats = -1
    n_steps_skip = 0
    analyze_conserved = .false.
    center_snapshots = .false.

    variance_coef = 1.d0
    k_B = 1.d0
    filtering_width = 0
    stoch_stress_form = 1

    need_inputs = .true.

    farg = 1
    if ( need_inputs .AND. narg >= 1 ) then
       call get_command_argument(farg, value = fname)
       inquire(file = fname, exist = lexist )
       if ( lexist ) then
          farg = farg + 1
          un = unit_new()
          open(unit=un, file = fname, status = 'old', action = 'read')
          read(unit=un, nml = probin_common)
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

       case ('--cfl')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) cfl

       case ('--plot_int')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) plot_int

       case ('--chk_int')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) chk_int

       case ('--prob_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) prob_type

       case ('--seed')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) seed

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

       case ('--boussinesq_beta')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) boussinesq_beta

       case ('--visc_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) visc_type

       case ('--visc_coef')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) visc_coef

       case ('--diff_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) diff_type

       case ('--advection_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) advection_type

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

       case ('--hydro_grid_int')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) hydro_grid_int

       case ('--project_dir')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) project_dir

       case ('--max_grid_projection_1')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) max_grid_projection(1)
       case ('--max_grid_projection_2')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) max_grid_projection(2)

       case ('--stats_int')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stats_int

       case ('--n_steps_save_stats')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) n_steps_save_stats

       case ('--n_steps_skip')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) n_steps_skip

       case ('--analyze_conserved')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) analyze_conserved

       case ('--center_snapshots')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) center_snapshots

       case ('--variance_coef')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) variance_coef

       case ('--k_B')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) k_B

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
    
  end subroutine probin_common_init

end module probin_common_module
