module probin_chemistry_module

  use bl_types

  use probin_common_module, only: MAX_SPECIES

  implicit none

  integer, parameter    :: max_reactions = 20
  integer, save         :: nreactions = 1           ! number of reactions

  ! stoichiometric factors for each reaction (species,LHS(1)/RHS(2),reaction)
  ! Example: For N1 + 2*N2 -> N3 use
  ! stoichiometric_factors(1:3,1,1) = 1 2 0
  ! stoichiometric_factors(1:3,2,1) = 0 0 1
  integer, save         :: stoichiometric_factors(MAX_SPECIES,2,max_reactions) = 0 

  ! reaction rate constant for each reaction (assuming Law of Mass Action holds)
  ! using rate_multiplier, reaction rates can be changed by the same factor
  ! if include_discrete_LMA_correction, n^2 and n^3 in rate expressions become
  ! n*(n-1/dv) and n*(n-1/dv)*(n-2/dv). 
  real(kind=dp_t), save :: rate_const(max_reactions) = 0.d0
  real(kind=dp_t), save :: rate_multiplier = 1.d0
  logical, save         :: include_discrete_LMA_correction = .true. 

  ! if n is positive, exclude species n (=solvent) when computing reaction rates
  ! in this case, the concentration of the solvent is assumed to be constant,
  ! which should be reflected on rate constants.
  ! if 0, no species is excluded
  ! e.g. U + S -> 2U, if exclude_solvent_comput_rates=0, rate=k*n_U*n_S
  !                   if exclude_solvent_comput_rates=2, rate=k_new*n_U where k_new=k*n_S
  integer, save         :: exclude_solvent_comput_rates = 0
  
  ! how to sample chemical production rates 
  integer, save         :: use_Poisson_rng = 1      ! 2=SSA
                                                    ! 1=do tau leaping (Poisson increments)
                                                    ! 0=do CLE (Gaussian increments)
                                                    ! -1=do deterministic chemistry

  namelist /probin_chemistry/ nreactions, stoichiometric_factors, rate_const, rate_multiplier
  namelist /probin_chemistry/ include_discrete_LMA_correction, exclude_solvent_comput_rates
  namelist /probin_chemistry/ use_Poisson_rng

contains

  subroutine probin_chemistry_init()

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
    logical            :: lexist, need_inputs
    
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
          open(unit = un, file = fname, status = 'old', action = 'read')
          read(unit = un, nml = probin_chemistry)
          close(unit = un)
          need_inputs = .false.
       end if
    end if
    
    ! also can be read in from the command line by appending 
    do while ( farg <= narg )
       call get_command_argument(farg, value = fname)
       select case (fname)

       case ('--nreactions')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) nreactions

       case ('--rate_multiplier')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) rate_multiplier

       case ('--include_discrete_LMA_correction')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) include_discrete_LMA_correction

       case ('--exclude_solvent_comput_rates')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) exclude_solvent_comput_rates 

       case ('--use_Poisson_rng')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) use_Poisson_rng

       case default
          if (parallel_IOProcessor() ) then
             print*,'probin_chemistry: command-line input ',trim(fname),' not read'
          end if

       end select

       farg = farg + 1
    end do
    
    ! check that nreactions<=max_reactions, otherwise abort with error message
    if (nreactions .gt. max_reactions) then 
       call bl_error(" nreactions greater than max_reactions - Aborting")
    end if
    
  end subroutine probin_chemistry_init

end module probin_chemistry_module
