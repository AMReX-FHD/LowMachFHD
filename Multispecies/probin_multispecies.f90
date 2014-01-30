module probin_multispecies_module

  use bl_types
  use bl_space
 
  implicit none

  integer, parameter :: max_species=10
  integer, parameter :: max_element=max_species*(max_species-1)/2
  integer, save      :: nspecies,max_step,init_type,inverse_type,timeinteg_type
  real(kind=dp_t)    :: cfl,chi,Temp,Press,start_time   ! chi is maximum eigenvalue of diffusion matrix
  real(kind=dp_t)    :: rho_in(2,max_species)     ! initial values for concentration, 2 for inside & outside circle
  real(kind=dp_t)    :: molmass_in(max_species) ! molar masses for nspecies
  real(kind=dp_t)    :: Dbar_in(max_element)    ! SM diffusion constant  
  real(kind=dp_t)    :: alpha,beta,fraction_tolerance  ! manufactured solution parameters
  integer            :: rho_part_bc_comp, mol_frac_bc_comp, diff_coeff_bc_comp 
  logical            :: print_error_norms, is_ideal_mixture, use_lapack
  
  namelist /probin_multispecies/ nspecies
  namelist /probin_multispecies/ max_step
  namelist /probin_multispecies/ cfl
  namelist /probin_multispecies/ chi
  namelist /probin_multispecies/ Temp
  namelist /probin_multispecies/ Press
  namelist /probin_multispecies/ fraction_tolerance
  namelist /probin_multispecies/ start_time
  namelist /probin_multispecies/ init_type
  namelist /probin_multispecies/ inverse_type   
  namelist /probin_multispecies/ timeinteg_type   
  namelist /probin_multispecies/ print_error_norms   
  namelist /probin_multispecies/ is_ideal_mixture   
  namelist /probin_multispecies/ use_lapack   
  namelist /probin_multispecies/ rho_in
  namelist /probin_multispecies/ molmass_in 
  namelist /probin_multispecies/ Dbar_in

contains

  subroutine probin_multispecies_init()

    use f2kcli
    use parallel
    use define_bc_module
    use bc_module
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
    nspecies          = 3
    max_step          = 10000
    cfl               = 1.0d0
    chi               = 1.0d0
    init_type         = 1
    inverse_type      = 1
    timeinteg_type    = 1
    print_error_norms = .true.
    is_ideal_mixture  = .true.
    use_lapack        = .true.
    rho_in            = 1.0d0
    molmass_in        = 1.0d0
    Dbar_in           = 1.0d0
    Temp              = 1.0d0
    Press             = 1.0d0
    fraction_tolerance= 1e-13 
    start_time        = 1.0d0
 
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
       case ('--max_step')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) max_step
       case ('--cfl')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) cfl
       case ('--chi')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) max_step
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
       case ('--rho_in')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) rho_in
       case ('--molmass_in')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) molmass_in
       case ('--Dbar_in')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) Dbar_in
       case ('--Temp')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) Temp
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
    endif
    
  end subroutine probin_multispecies_init

end module probin_multispecies_module
