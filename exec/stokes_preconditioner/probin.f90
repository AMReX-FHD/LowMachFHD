! This module stores the runtime parameters.  The probin_init() routine is
! used to initialize the runtime parameters

module probin_module

  use bl_types
  use bl_space

  implicit none

  ! For comments and instructions on how to set the input parameters see namelist section below
  !------------------------------------------------------------- 
  integer,save    :: mode_coefs(2),prob_coeff,prob_sol,prob_dir,test_type
  real(dp_t),save :: smoothing_width,var_coeff_mag(3),coeff_mag(3),coeff_ratio(3),ABC_coefs(3)

  !------------------------------------------------------------- 
  ! Input parameters controlled via namelist input, with comments
  !------------------------------------------------------------- 

  !----------------------
  namelist /probin/ prob_coeff    ! Initial setup/forcing etc selection 
  namelist /probin/ prob_sol      ! Initial setup/forcing etc selection 
  namelist /probin/ test_type     ! Check gmres, preconditioner, or discretization accuracy
                                  ! (see main_driver.f90 for key)

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
  !------------------------------------------------------------- 

contains

  subroutine probin_init()

    use f2kcli
    use parallel
    use bc_module
    use bl_IO_module
    use bl_prof_module
    use bl_error_module
    use bl_constants_module
    use cluster_module
    
    integer    :: narg, farg

    character(len=128) :: fname
    character(len=128) :: probin_env

    logical :: lexist

    integer :: un, ierr

    narg = command_argument_count()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize the runtime parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Defaults
    prob_coeff = 2
    prob_sol = 2
    test_type = 1

    smoothing_width = 1.d0
    mode_coefs(1:2) = 1
    prob_dir = 1
    ABC_coefs(1:3) = 1.d0

    var_coeff_mag(1:3) = 0.d0
    coeff_mag(1:3) = 1.d0
    coeff_ratio(1:3) = 1.d0

    farg = 1
    if (narg >= 1) then
       call get_command_argument(farg, value = fname)
       inquire(file = fname, exist = lexist )
       if ( lexist ) then
          farg = farg + 1
          un = unit_new()
          open(unit=un, file = fname, status = 'old', action = 'read')
          read(unit=un, nml = probin)
          close(unit=un)
       end if
    end if

    ! stuff that can be read in from the command line by appending, e.g., "--prob_hi_x 64.0"
    do while ( farg <= narg )
       call get_command_argument(farg, value = fname)
       select case (fname)

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
    
  end subroutine probin_init

end module probin_module
