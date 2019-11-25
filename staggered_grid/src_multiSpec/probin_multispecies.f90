module probin_multispecies_module

  use bl_types
  use probin_common_module, only: MAX_SPECIES
 
  implicit none

  integer, parameter :: max_element=MAX_SPECIES*(MAX_SPECIES-1)/2  

  integer, save      :: inverse_type,temp_type,chi_iterations
  real(kind=dp_t)    :: start_time
  real(kind=dp_t)    :: T_init(2) 
  real(kind=dp_t)    :: Dbar(max_element)
  real(kind=dp_t)    :: Dtherm(MAX_SPECIES)
  real(kind=dp_t)    :: H_offdiag(max_element)
  real(kind=dp_t)    :: H_diag(MAX_SPECIES)
  real(kind=dp_t)    :: fraction_tolerance
  logical            :: correct_flux,print_error_norms,plot_stag
  logical            :: is_nonisothermal,is_ideal_mixture,use_lapack
  logical            :: use_multiphase
  real(kind=dp_t)    :: kc_tension, alpha_gex
  real(kind=dp_t)    :: c_init(2,MAX_SPECIES)
  real(kind=dp_t)    :: c_bc(3,2,MAX_SPECIES)
  real(kind=dp_t)    :: alpha1,beta,delta,sigma     ! manufactured solution parameters populated in init
  integer            :: midpoint_stoch_mass_flux_type
  integer            :: avg_type
  integer            :: mixture_type

  real(kind=dp_t)    :: rhotot_bc(3,2)              ! used for boussinesq algorithm - not part of namelist
                                                    ! code sets this in src_lowMach/define_bc_tower.f90

  ! Physical properties:
  !----------------------
  namelist /probin_multispecies/ fraction_tolerance ! For roundoff errors in mass and mole fractions
  namelist /probin_multispecies/ start_time
  namelist /probin_multispecies/ inverse_type       ! Only for LAPACK:  1=inverse, 2=pseudo inverse
  namelist /probin_multispecies/ correct_flux       ! Manually ensure mass is conserved to roundoff 
  namelist /probin_multispecies/ print_error_norms   
  namelist /probin_multispecies/ is_ideal_mixture   ! If T assume Gamma=I (H=0) and simplify
  namelist /probin_multispecies/ is_nonisothermal   ! If T Soret effect will be included
  namelist /probin_multispecies/ use_lapack         ! Use LAPACK or iterative method for diffusion matrix (recommend False)
  namelist /probin_multispecies/ chi_iterations     ! number of iterations used in Dbar2chi_iterative
  namelist /probin_multispecies/ use_multiphase     ! whether to enable multiphase terms
  namelist /probin_multispecies/ kc_tension         ! interfacial tension parameter
  namelist /probin_multispecies/ alpha_gex          ! parameter for excess gibbs free energy

  ! Initial and boundary conditions 
  !----------------------
  namelist /probin_multispecies/ T_init     ! initial values for temperature (bottom/top, inside/outside circle, etc.)
  namelist /probin_multispecies/ temp_type  ! for initializing temperature
  namelist /probin_multispecies/ c_init     ! initial values for c
  namelist /probin_multispecies/ c_bc       ! c_i boundary conditions (dir,lohi,species)
  
  ! Thermodynamic and transport properties:
  !----------------------

  ! These are lower-triangules of symmetric matrices represented as vectors
  ! Number of elements is (nspecies*(nspecies-1)/2)
  ! The values are read row by row starting from top going down (this allows easy addition/deletion of new species/rows)
  ! So D_12; D_13, D_23; D_14, D_24, D_34; ...
  namelist /probin_multispecies/ Dbar       ! Maxwell-Stefan diffusion constant  
  namelist /probin_multispecies/ Dtherm     ! thermo-diffusion coefficients, only differences among elements matter
  namelist /probin_multispecies/ H_offdiag
  namelist /probin_multispecies/ H_diag     ! Diagonal of H=d^2F/dx^2, these are vectors of length nspecies
  namelist /probin_multispecies/ plot_stag  ! plot staggered velocities in separate plotfile

  ! Algorithm control
  !----------------------
  namelist /probin_multispecies/ midpoint_stoch_mass_flux_type  ! 1 = Strato
                                                                ! 2 = Ito

  namelist /probin_multispecies/ avg_type   ! how to compute stochastc_mass_fluxdiv
                                            ! 1=arithmetic (with C0-Heaviside), 2=geometric, 3=harmonic
                                            ! 10=arithmetic average with discontinuous Heaviside function
                                            ! 11=arithmetic average with C1-smoothed Heaviside function
                                            ! 12=arithmetic average with C2-smoothed Heaviside function

  namelist /probin_multispecies/ mixture_type ! Model for how transport and thermodynamic coefficients depend on composition
                                              ! See compute_mixture_properties.f90 for values supported at present
                                              ! The default mixture_type=0 means no dependence on composition

contains

  subroutine probin_multispecies_init()

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

    ! here we set some random values to be replaced from the input file
    fraction_tolerance = 1.d-14  ! must be larger than machine eps of else the W=(1,0) case fails)
    start_time         = 0.0d0 
    inverse_type       = 1
    correct_flux       = .true.
    print_error_norms  = .true.
    is_ideal_mixture   = .true.
    is_nonisothermal   = .false.
    use_lapack         = .false.
    use_multiphase     = .true.
    kc_tension         = 0.d0
    alpha_gex          = 0.d0
    chi_iterations     = 10
    T_init             = 1.0d0
    temp_type          = 0
    c_init             = 1.0d0
    c_bc               = 0.d0
    Dbar               = 1.0d0
    Dtherm             = 0.0d0
    H_offdiag          = 0.0d0
    H_diag             = 0.0d0
    plot_stag          = .false.
    midpoint_stoch_mass_flux_type = 1
    avg_type           = 1
    mixture_type       = 0
 
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
          read(unit=un, nml = probin_multispecies)
          close(unit=un)
          need_inputs = .false.
       end if
    end if

    ! also can be read in from the command line by appending 
    do while ( farg <= narg )
       call get_command_argument(farg, value = fname)
       select case (fname)

       case ('--fraction_tolerance')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) fraction_tolerance

       case ('--start_time')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) start_time

       case ('--inverse_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) inverse_type

       case ('--correct_flux')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) correct_flux

       case ('--print_error_norms')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) print_error_norms

       case ('--is_ideal_mixture')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) is_ideal_mixture

       case ('--is_nonisothermal')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) is_nonisothermal

       case ('--use_lapack')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) use_lapack

       case ('--T_init_1')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) T_init(1)

       case ('--T_init_2')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) T_init(2)

       case ('--chi_iterations')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) chi_iterations

       case ('--temp_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) temp_type

       case ('--plot_stag')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) plot_stag

       case ('--c_init1_1')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_init(1,1)

       case ('--c_init1_2')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_init(1,2)

       case ('--c_init1_3')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_init(1,3)

       case ('--c_init1_4')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_init(1,4)

       case ('--c_init2_1')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_init(2,1)

       case ('--c_init2_2')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_init(2,2)

       case ('--c_init2_3')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_init(2,3)

       case ('--c_init2_4')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_init(2,4)

       case ('--midpoint_stoch_mass_flux_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) midpoint_stoch_mass_flux_type

       case ('--avg_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) avg_type

       case ('--mixture_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mixture_type

       case ('--')
          farg = farg + 1
          exit
       case default

       end select
       farg = farg + 1
    end do
    
  end subroutine probin_multispecies_init

end module probin_multispecies_module
