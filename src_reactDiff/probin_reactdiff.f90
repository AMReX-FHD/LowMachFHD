module probin_reactdiff_module

  use bl_types
  use probin_common_module, only: MAX_SPECIES
 
  implicit none

  ! Control of algorithm
  !----------------------
  integer, save :: temporal_integrator = 0          ! 0=D + R (first-order splitting)
                                                    ! 1=(1/2)R + D + (1/2)R (Strang option 1)
                                                    ! 2=(1/2)D + R + (1/2)D (Strang option 2)
                                                    ! -1=unsplit forward Euler
                                                    ! -2=unsplit explicit midpoint 
                                                    ! -3=unsplit multinomial diffusion
                                                    ! -4=unsplit implicit midpoint

  integer, save :: diffusion_type = 0               ! only used for split schemes (temporal_integrator>=0)
                                                    ! 0=explicit trapezoidal predictor/corrector
                                                    ! 1=Crank-Nicolson semi-implicit
                                                    ! 2=explicit midpoint
                                                    ! 3=multinomial diffusion
                                                    ! 4=forward Euler  

  integer, save :: reaction_type = 0                ! only used for split schemes (temporal_integrator>=0)
                                                    ! 0=first-order (deterministic, tau leaping, CLE, or SSA)
                                                    ! 1=second-order (determinisitc, tau leaping, or CLE only)

  integer, save :: midpoint_stoch_flux_type = 1     ! only used for midpoint diffusion schemes (split as well as unsplit)
                                                    ! corrector formulation of noise
                                                    ! 1 = K(nold) * W1 + K(nold)         * W2
                                                    ! 2 = K(nold) * W1 + K(npred)        * W2
                                                    ! 3 = K(nold) * W1 + K(2*npred-nold) * W2

  integer, save :: avg_type = 1                     ! how to compute n on faces for stochastic weighting
                                                    ! 1=arithmetic (with C0-Heaviside), 2=geometric, 3=harmonic
                                                    ! 10=arithmetic average with discontinuous Heaviside function
                                                    ! 11=arithmetic average with C1-smoothed Heaviside function
                                                    ! 12=arithmetic average with C2-smoothed Heaviside function

  logical, save :: inhomogeneous_bc_fix = .false.   ! use the Einkemmer boundary condition fix (split schemes only)

  ! Initial and boundary conditions
  !----------------------
  real(kind=dp_t), save    :: n_init_in(2,MAX_SPECIES) = 1.d0  ! initial values to be used in init_n.f90

  integer, save            :: model_file_init = 0              ! initialize from model files:
                                                               ! 0=no, 1=usual order (Fortran), -1=transpose order (C)

  character(len=128), save :: model_file(MAX_SPECIES)          ! one model file for each species
  
  logical, save :: integer_populations = .false.               ! initialize with all number of molecules strictly integer

  real(kind=dp_t), save    :: n_bc(3,2,MAX_SPECIES) = 0.d0     ! n_i boundary conditions (dir,lohi,species)

  ! Diffusion     
  !----------------------                          
  real(kind=dp_t), save :: D_Fick(MAX_SPECIES) = 1.d0          ! Fickian diffusion coeffs

  integer, save         :: diffusion_stencil_order = 1         ! diffusion boundary stencil order

  integer, save         :: mg_verbose = 0                      ! implicit diffusion solve verbosity

  integer, save         :: cg_verbose = 0                      ! implicit diffusion solve bottom solver verbosity

  real(kind=dp_t), save :: implicit_diffusion_rel_eps = 1.d-10 ! relative eps for implicit diffusion solve

  real(kind=dp_t), save :: implicit_diffusion_abs_eps = -1.d0  ! absolute eps for implicit diffusion solve
  
  ! Controlling output
  !----------------------                          
  integer, save :: n_steps_write_avg = 0 ! If non-zero, its absolute value tells how many steps before writing total densites
                                         ! If positive, it writes average number densities in the system
                                         ! If negative, it writes the total number of molecules in the system
  
  namelist /probin_reactdiff/ temporal_integrator, diffusion_type, reaction_type
  namelist /probin_reactdiff/ midpoint_stoch_flux_type, avg_type, inhomogeneous_bc_fix
  namelist /probin_reactdiff/ n_init_in, model_file_init, model_file, integer_populations, n_bc
  namelist /probin_reactdiff/ D_Fick, diffusion_stencil_order, mg_verbose, cg_verbose
  namelist /probin_reactdiff/ implicit_diffusion_rel_eps, implicit_diffusion_abs_eps
  namelist /probin_reactdiff/ n_steps_write_avg

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
          read(unit = un, nml = probin_reactdiff)
          close(unit = un)
          need_inputs = .false.
       end if
    end if
    
    ! also can be read in from the command line by appending 
    do while ( farg <= narg )
       call get_command_argument(farg, value = fname)
       select case (fname)

       case ('--temporal_integrator')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) temporal_integrator 

       case ('--diffusion_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) diffusion_type

       case ('--reaction_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) reaction_type

       case ('--midpoint_stoch_flux_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) midpoint_stoch_flux_type

       case ('--avg_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) avg_type

       case ('--inhomogeneous_bc_fix')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) inhomogeneous_bc_fix

       case ('--model_file_init')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) model_file_init

       case ('--integer_populations')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) integer_populations 

       case ('--D_Fick_1')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) D_Fick(1)

       case ('--D_Fick_2')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) D_Fick(2)

       case ('--D_Fick_3')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) D_Fick(3)

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

       case ('--n_steps_write_avg')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) n_steps_write_avg 

       case default
          if (parallel_IOProcessor() ) then
             print*,'probin_reactdiff: command-line input ',trim(fname),' not read'
          end if

       end select

       farg = farg + 1
    end do
    
    if (inhomogeneous_bc_fix .and. temporal_integrator .lt. 0) then
       call bl_error("inhomogeneous_bc_fix only appropriate for split schemes")
    end if

  end subroutine probin_reactdiff_init

end module probin_reactdiff_module
