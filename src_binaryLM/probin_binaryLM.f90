! This module stores the runtime parameters.  The probin_binarylm_init() routine is
! used to initialize the runtime parameters

module probin_binarylm_module

  use bl_types

  implicit none

  ! For comments and instructions on how to set the input parameters see 
  ! namelist section below
  !------------------------------------------------------------- 
  real(dp_t), save :: c_init(2)
  real(dp_t), save :: temperature
  real(dp_t), save :: material_properties(3,3),c_bc(3,2)
  integer   , save :: barodiffusion_type
  logical   , save :: analyze_binary,plot_stag
  real(dp_t), save :: boussinesq_beta

  !------------------------------------------------------------- 
  ! Input parameters controlled via namelist input, with comments
  !------------------------------------------------------------- 

  ! problem setup
  namelist /probin_binarylm/ c_init            ! controls initial concentration range
  namelist /probin_binarylm/ c_bc              ! c boundary conditions (dir,face).
                                               ! Dirichlet for RESERVOIR; Neumann for WALL
  ! fluid properties
  namelist /probin_binarylm/ temperature       ! temperature
  namelist /probin_binarylm/ material_properties ! Coefficients A/B/C for chi/eta/kappa
     ! Formula is: indx=1 for chi, indx=2 for eta, indx=3 for kappa (NOT implemented yet)
     ! coeff=coeff0*(material_properties(1,indx) + material_properties(2,indx)*c) / (1.d0 + material_properties(3,indx)*c)

  namelist /probin_binarylm/ boussinesq_beta    ! beta for boussinesq gravity

  namelist /probin_binarylm/ barodiffusion_type ! 0 = no barodiffusion
                                                ! 1 = fixed gradp from initialization
                                                ! 2 = update gradp each time step
  namelist /probin_binarylm/ analyze_binary     ! Call the older analyze_spectra_binary or the new analyze_spectra?

  namelist /probin_binarylm/ plot_stag          ! include staggered plotfiles

                             

contains

  subroutine probin_binarylm_init()

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

    logical :: lexist

    integer :: un

    narg = command_argument_count()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize the runtime parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Defaults

    c_init(1:2) = 1.d0
    c_bc(1:3,1:2) = 0.d0

    temperature = 1.d0
    material_properties(1:3,1:3) = 0.d0

    boussinesq_beta = 0.d0

    barodiffusion_type = 0
    
    analyze_binary=.true.
    plot_stag = .false.

    farg = 1
    if (narg >= 1) then
       call get_command_argument(farg, value = fname)
       inquire(file = fname, exist = lexist )
       if ( lexist ) then
          farg = farg + 1
          un = unit_new()
          open(unit=un, file = fname, status = 'old', action = 'read')
          read(unit=un, nml = probin_binarylm)
          close(unit=un)
       end if
    end if

    ! stuff that can be read in from the command line by appending, e.g., "--prob_hi_x 64.0"
    do while ( farg <= narg )
       call get_command_argument(farg, value = fname)
       select case (fname)

       case ('--c_init_1')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_init(1)
       case ('--c_init_2')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_init(2)

       case ('--c_bc_x_lo')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_bc(1,1)
       case ('--c_bc_x_hi')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_bc(1,2)
       case ('--c_bc_y_lo')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_bc(2,1)
       case ('--c_bc_y_hi')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_bc(2,2)
       case ('--c_bc_z_lo')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_bc(3,1)
       case ('--c_bc_z_hi')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) c_bc(3,2)

       case ('--material_properties_1_1')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) material_properties(1,1)
       case ('--material_properties_2_1')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) material_properties(2,1)
       case ('--material_properties_3_1')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) material_properties(3,1)

       case ('--material_properties_1_2')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) material_properties(1,2)
       case ('--material_properties_2_2')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) material_properties(2,2)
       case ('--material_properties_3_2')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) material_properties(3,2)

       case ('--material_properties_1_3')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) material_properties(1,3)
       case ('--material_properties_2_3')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) material_properties(2,3)
       case ('--material_properties_3_3')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) material_properties(3,3)

       case ('--boussinesq_beta')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) boussinesq_beta

       case ('--barodiffusion_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) barodiffusion_type

       case ('--analyze_binary')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) analyze_binary

       case ('--plot_stag')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) plot_stag

       case ('--')
          farg = farg + 1
          exit

       case default

       end select

       farg = farg + 1
    end do
    
  end subroutine probin_binarylm_init

end module probin_binarylm_module
