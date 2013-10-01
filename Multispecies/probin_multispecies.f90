module probin_multispecies_module

  use bl_types
  use bl_space
 
  implicit none

  integer, parameter :: max_species=10
  integer, parameter :: max_element=max_species*(max_species-1)/2
  integer, save      :: nspecies,max_step,init_type,inverse_type,timeinteg_type
  real(kind=dp_t)    :: chi                 !maximum eigenvalue of diffusion matrix 
  real(kind=dp_t)    :: c_bc(2,max_species) !initial values for concentration, 2 for inside & outside circle
  real(kind=dp_t)    :: d_bc(max_species)   !initial values for diffusivities, presently scalar numbers 
  real(kind=dp_t)    :: m_bc(max_species)   ! masses of nspecies
  real(kind=dp_t)    :: Dbar_bc(max_element)! SM diffusion constant  
  integer            :: rho_part_bc_comp, mol_frac_bc_comp, diff_coeff_bc_comp 
  
  namelist /probin_multispecies/ nspecies
  namelist /probin_multispecies/ max_step
  namelist /probin_multispecies/ chi
  namelist /probin_multispecies/ init_type
  namelist /probin_multispecies/ inverse_type   
  namelist /probin_multispecies/ timeinteg_type   
  namelist /probin_multispecies/ c_bc
  namelist /probin_multispecies/ d_bc
  namelist /probin_multispecies/ m_bc
  namelist /probin_multispecies/ Dbar_bc

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
    logical            :: lexist
    logical            :: need_inputs
    integer            :: un

    narg = command_argument_count()
 
    ! here we set some random values to be replaced from the input file
    nspecies  = 3
    max_step  = 10000
    chi       = 1.0d0
    init_type = 1
    inverse_type = 1
    timeinteg_type = 1
    c_bc      = 1.0d0
    m_bc      = 1.0d0
    Dbar_bc   = 1.0d0
 
    ! bc_tower strcuture in memory 
    ! 1-3 = velocity, 4 = Pressure, rho_tot = scal_bc_comp, rho_i = rhot_tot+1,
    ! mol_frac = rho_tot+2, diff_coeff=tran_bc_comp
    rho_part_bc_comp = scal_bc_comp + 1
    mol_frac_bc_comp = scal_bc_comp + 2
    diff_coeff_bc_comp = tran_bc_comp
 
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

    ! stuff that can be read in from the command line by appending 
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
      case ('--init_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) init_type
      case ('--inverse_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) inverse_type
      case ('--integrator_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) timeinteg_type
      case ('--')
          farg = farg + 1
          exit
      case default

      end select

      farg = farg + 1
    end do
    
  end subroutine probin_multispecies_init

end module probin_multispecies_module
