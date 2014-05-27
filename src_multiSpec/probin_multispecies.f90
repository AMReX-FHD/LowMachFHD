module probin_multispecies_module

  use bl_types
  use bl_space
 
  implicit none

  integer, parameter :: max_species=10
  integer, parameter :: max_element=max_species*(max_species-1)/2
  integer, save      :: nspecies,init_type,inverse_type,timeinteg_type
  real(kind=dp_t)    :: diff_coef,Press,start_time ! diff_coef is maximum eigenvalue of diffusion matrix
  real(kind=dp_t)    :: rho_init(2,max_species) ! initial values for concentration, 2 for inside & outside circle
  real(kind=dp_t)    :: T_init(2) ! initial values for temperature (bottom/top, inside/outside circle, etc.)
  real(kind=dp_t)    :: molmass(max_species)     ! molar masses for nspecies
  real(kind=dp_t)    :: rhobar(max_species)      ! rhobar for nspecies
  real(kind=dp_t)    :: Dbar(max_element)        ! SM diffusion constant  
  real(kind=dp_t)    :: Dtherm(max_element)          ! thermo-diffusion coefficients  
  real(kind=dp_t)    :: alpha1,beta,delta,sigma     ! manufactured solution parameters populated in init
  real(kind=dp_t)    :: variance_coef_mass,fraction_tolerance  
  integer            :: rho_part_bc_comp, mol_frac_bc_comp, temp_bc_comp ! not input: populated at main 
  logical            :: correct_flux,use_stoch,print_error_norms
  logical            :: is_nonisothermal,is_ideal_mixture,use_lapack
  real(kind=dp_t)    :: rho_bc(3,2,max_species)
  
  namelist /probin_multispecies/ nspecies
  namelist /probin_multispecies/ diff_coef
  namelist /probin_multispecies/ variance_coef_mass
  namelist /probin_multispecies/ Press
  namelist /probin_multispecies/ fraction_tolerance
  namelist /probin_multispecies/ start_time
  namelist /probin_multispecies/ init_type
  namelist /probin_multispecies/ inverse_type   
  namelist /probin_multispecies/ timeinteg_type   
  namelist /probin_multispecies/ correct_flux   
  namelist /probin_multispecies/ use_stoch   
  namelist /probin_multispecies/ print_error_norms   
  namelist /probin_multispecies/ is_ideal_mixture   
  namelist /probin_multispecies/ is_nonisothermal   
  namelist /probin_multispecies/ use_lapack   
  namelist /probin_multispecies/ rho_init
  namelist /probin_multispecies/ molmass
  namelist /probin_multispecies/ rhobar
  namelist /probin_multispecies/ Dbar
  namelist /probin_multispecies/ Dtherm
  namelist /probin_multispecies/ T_init
  namelist /probin_multispecies/ rho_bc ! rho_i boundary conditions (dir,lohi,species)

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
    nspecies           = 2 
    diff_coef                = 1.0d0
    variance_coef_mass = 1E-3
    Press              = 1.0d0
    fraction_tolerance = 1e-13 
    start_time         = 0.0d0 
    init_type          = 0
    inverse_type       = 1
    timeinteg_type     = 1
    correct_flux       = .true.
    use_stoch          = .true.
    print_error_norms  = .true.
    is_ideal_mixture   = .true.
    is_nonisothermal   = .true.
    use_lapack         = .true.
    rho_init           = 1.0d0
    molmass            = 1.0d0
    rhobar             = 1.d0
    Dbar               = 1.0d0
    Dtherm             = 1.0d0
    T_init             = 1.0d0
    rho_bc             = 0.d0
 
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

       case ('--nspecies')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) nspecies
       case ('--diff_coef')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) diff_coef
       case ('--init_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) init_type
       case ('--inverse_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) inverse_type
       case ('--timeinteg_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) timeinteg_type
       case ('--rho_init')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) rho_init
       case ('--molmass')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) molmass
       case ('--rhobar')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) rhobar
       case ('--Dbar')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) Dbar
       case ('--Dtherm')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) Dtherm
       case ('--Press')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) Press
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
    
  end subroutine probin_multispecies_init

end module probin_multispecies_module
