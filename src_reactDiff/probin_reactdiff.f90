module probin_reactdiff_module

  use bl_types
  use bl_space
  use probin_common_module, only: dim_in
 
  implicit none

  integer, parameter :: max_species=10
  integer, parameter :: max_reactions=20

  integer, save         :: nspecies = 2
  integer, save         :: nreactions = 0
  
  ! Control of algorithm:
  integer, save         :: diffusion_type = 0 ! 0=explicit trapezoidal predictor/corrector
                                              ! 1=Crank-Nicolson semi-implicit
                                              ! 2=explicit midpoint
  integer, save         :: reaction_type = 0  ! 0=first-order tau leaping or CLE
                                              ! 1=second-order tau leaping or CLE
                                              ! 2=SSA
  logical, save         :: use_Poisson_rng = .true. ! Only used for reaction_type = 0, 1
                                                    ! If true do tau leaping (Poisson increments),
                                                    ! otherwise do Chemical Langevin Equation (CLE) (Gaussian increments)
  integer, save         :: splitting_type = 0 ! 0=D + R
                                              ! 1=(1/2)R + D + (1/2)R
                                              ! 2=(1/2)D + R + (1/2)D
  integer, save         :: avg_type = 3 ! compute n on faces for stochastic weighting
                                        ! 1=arithmetic, 2=geometric, 3=harmonic
  
  ! Initial and boundary conditions:
  real(kind=dp_t), save :: n_init_in(2,max_species) = 1.d0 ! Initial values to be used in init_n.f90
  real(kind=dp_t), save :: n_bc(3,2,max_species) = 0.d0 ! n_i boundary conditions (dir,lohi,species)
  
  ! Diffusion                                      
  real(kind=dp_t), save :: D_Fick(max_species) = 1.d0 ! Fickian diffusion coeffs
  integer, save         :: mg_verbose = 0 ! implicit diffusion solve verbosity
  integer, save         :: cg_verbose = 0 ! implicit diffusion solve bottom solver verbosity
  real(kind=dp_t), save :: implicit_diffusion_rel_eps = 1.d-10 ! relative epsilon for implicit diffusion solve
  real(kind=dp_t), save :: implicit_diffusion_abs_eps = -1.d0  ! absolute epsilon for implicit diffusion solve
  
  ! Chemical reactions:
  ! Donev: Added rates here as input and a new flag, see compute_reaction_rates
  logical, save         :: include_discrete_LMA_correction = .true. ! Whether to compute chemical rates using classical LMA or integer-based one
  real(kind=dp_t), save :: chemical_rates(max_reactions)=0.0d0 ! LMA chemical reaction rate for each reaction (assuming Law of Mass holds)
  integer, save         :: stoichiometric_factors(max_species,2,max_reactions) = 0 ! stoichiometric factors for each reaction
                                                                                   ! (species,LHS(1)/RHS(2),reaction)
  
  namelist /probin_reactdiff/ nspecies, nreactions
  namelist /probin_reactdiff/ diffusion_type, reaction_type, splitting_type, avg_type
  namelist /probin_reactdiff/ D_Fick, n_init_in, n_bc
  namelist /probin_reactdiff/ mg_verbose, cg_verbose
  namelist /probin_reactdiff/ implicit_diffusion_rel_eps, implicit_diffusion_abs_eps
  namelist /probin_reactdiff/ stoichiometric_factors, use_Poisson_rng, chemical_rates, include_discrete_LMA_correction

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

       case ('--splitting_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) splitting_type

       case ('--avg_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) avg_type

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
