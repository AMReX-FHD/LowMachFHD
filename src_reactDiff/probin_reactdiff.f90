module probin_reactdiff_module

  use bl_types
  use bl_space
  use probin_common_module, only: dim_in
 
  implicit none

  integer, parameter :: max_species=10
  integer, parameter :: max_reactions=20

  ! Problem description
  !----------------------
  integer, save         :: nspecies = 2             ! number of species
  integer, save         :: nreactions = 1           ! number of reactions
  
  ! Control of algorithm
  !----------------------
  integer, save   :: diffusion_type = 0             ! 0=explicit trapezoidal predictor/corrector
                                                    ! 1=Crank-Nicolson semi-implicit
                                                    ! 2=explicit midpoint
  integer, save   :: reaction_type = 0              ! 0=first-order tau leaping or CLE
                                                    ! 1=second-order tau leaping or CLE
                                                    ! 2=SSA
  integer, save   :: use_Poisson_rng = 1            ! Only used for reaction_type = 0, 1
                                                    ! If -1 do deterministic chemistry
                                                    ! If 1 do tau leaping (Poisson increments),
                                                    ! If 0 do Chemical Langevin Equation (CLE) (Gaussian increments)
  integer, save   :: splitting_type = 0             ! 0=D + R (first-order splitting)
                                                    ! 1=(1/2)R + D + (1/2)R (Strang option 1)
                                                    ! 2=(1/2)D + R + (1/2)D (Strang option 2)
  logical, save   :: inhomogeneous_bc_fix = .false. ! use the Einkemmer boundary condition fix
  integer, save   :: avg_type = 3                   ! how to compute n on faces for stochastic weighting
                                                    ! 1=arithmetic, 2=geometric, 3=harmonic
  
  ! Initial and boundary conditions
  !----------------------
  real(kind=dp_t), save :: n_init_in(2,max_species) = 1.d0 ! Initial values to be used in init_n.f90
  real(kind=dp_t), save :: n_bc(3,2,max_species) = 0.d0    ! n_i boundary conditions (dir,lohi,species)

  ! Diffusion     
  !----------------------                          
  real(kind=dp_t), save :: D_Fick(max_species) = 1.d0          ! Fickian diffusion coeffs
  integer, save         :: diffusion_stencil_order = 1         ! diffusion boundary stencil order
  integer, save         :: mg_verbose = 0                      ! implicit diffusion solve verbosity
  integer, save         :: cg_verbose = 0                      ! implicit diffusion solve bottom solver verbosity
  real(kind=dp_t), save :: implicit_diffusion_rel_eps = 1.d-10 ! relative eps for implicit diffusion solve
  real(kind=dp_t), save :: implicit_diffusion_abs_eps = -1.d0  ! absolute eps for implicit diffusion solve
  
  ! Chemical reactions
  !----------------------
  real(kind=dp_t), save :: cross_section = 1.d0 ! thickness (in 2D) of cell

  ! Whether to compute chemical rates using classical LMA or integer-based one
  logical, save         :: include_discrete_LMA_correction = .true. 

  ! LMA chemical reaction rate for each reaction (assuming Law of Mass holds)
  real(kind=dp_t), save :: chemical_rates(max_reactions) = 0.0d0, rate_multiplier=1.0d0

  ! stoichiometric factors for each reaction (species,LHS(1)/RHS(2),reaction)
  ! Example: For N1 + 2*N2 -> N3 use
  ! stoichiometric_factors(1:3,1,1) = 1 2 0
  ! stoichiometric_factors(1:3,2,1) = 0 0 1
  integer, save         :: stoichiometric_factors(max_species,2,max_reactions) = 0 
  
  namelist /probin_reactdiff/ nspecies, nreactions
  namelist /probin_reactdiff/ diffusion_type, reaction_type, use_Poisson_rng, splitting_type, avg_type
  namelist /probin_reactdiff/ inhomogeneous_bc_fix, n_init_in, n_bc
  namelist /probin_reactdiff/ D_Fick, diffusion_stencil_order, mg_verbose, cg_verbose
  namelist /probin_reactdiff/ implicit_diffusion_rel_eps, implicit_diffusion_abs_eps
  namelist /probin_reactdiff/ cross_section, include_discrete_LMA_correction
  namelist /probin_reactdiff/ chemical_rates, rate_multiplier, stoichiometric_factors

contains

  subroutine probin_reactdiff_init()

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

    ! You can put default values here if you want, but we have specified them above 
    ! in the variable declaration
 
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
          read(unit=un, nml = probin_reactdiff)
          close(unit=un)
          need_inputs = .false.
       end if
    end if
    
    ! also can be read in from the command line by appending 
    do while ( farg <= narg )
       call get_command_argument(farg, value = fname)
       select case (fname)

       case ('--nspecies')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) nspecies

       case ('--nreactions')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) nreactions

       case ('--diffusion_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) diffusion_type

       case ('--reaction_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) reaction_type

       case ('--use_Poisson_rng')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) use_Poisson_rng

       case ('--splitting_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) splitting_type

       case ('--inhomogeneous_bc_fix')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) inhomogeneous_bc_fix

       case ('--avg_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) avg_type

       case ('--diffusion_stencil_order')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) diffusion_stencil_order

       case ('--mg_verbose')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mg_verbose

       case ('--cg_verbose')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) cg_verbose

       case ('--implicit_diffusion_rel_eps')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) implicit_diffusion_rel_eps

       case ('--implicit_diffusion_abs_eps')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) implicit_diffusion_abs_eps

       case ('--cross_section')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) cross_section

       case ('--include_discrete_LMA_correction')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) include_discrete_LMA_correction

       case ('--')
          farg = farg + 1
          exit
       case default

       end select
       farg = farg + 1
    end do
    
    ! check that nspecies<=max_species, otherwise abort with error message
    if(nspecies.gt.max_species) then 
       call bl_error(" nspecies greater than max_species - Aborting")
       stop
    end if
    
    ! check that nreactions<=max_reactions, otherwise abort with error message
    if(nreactions.gt.max_reactions) then 
       call bl_error(" nreactions greater than max_reactions - Aborting")
       stop
    end if
    
  end subroutine probin_reactdiff_init

end module probin_reactdiff_module
