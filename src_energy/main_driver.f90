subroutine main_driver()

  use bl_IO_module
  use init_lowmach_module
  use compute_mixture_properties_module
  use initial_projection_module
  use write_plotfile_module
  use advance_timestep_inertial_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use multifab_physbc_stag_module
  use analysis_module
  use analyze_spectra_module
  use div_and_grad_module
  use eos_check_module
  use estdt_module
  use stag_mg_layout_module
  use macproject_module
  use stochastic_mass_fluxdiv_module
  use stochastic_m_fluxdiv_module
  use fill_rho_ghost_cells_module
  use ParallelRNGs 
  use bl_rng_module
  use bl_random_module
  use mass_flux_utilities_module
  use compute_HSE_pres_module
  use convert_m_to_umac_module
  use sum_momenta_module
  use restart_module
  use checkpoint_module
  use reservoir_bc_fill_module
  use probin_common_module, only: prob_lo, prob_hi, n_cells, dim_in, hydro_grid_int, &
                                  max_grid_size, n_steps_save_stats, n_steps_skip, &
                                  plot_int, chk_int, seed, stats_int, bc_lo, bc_hi, restart, &
                                  probin_common_init, print_int, nspecies, &
                                  advection_type, fixed_dt, max_step, cfl, &
                                  algorithm_type, variance_coef_mom, initial_variance, &
                                  variance_coef_mass, barodiffusion_type, use_bl_rng
  use probin_multispecies_module, only: Dbar, start_time, probin_multispecies_init
  use probin_gmres_module, only: probin_gmres_init
  use probin_charged_module, only: probin_charged_init, use_charged_fluid, dielectric_const, &
                                   dielectric_type
  use probin_chemistry_module, only: probin_chemistry_init, nreactions
  use probin_energy_module, only: probin_energy_init

  implicit none

  ! quantities will be allocated with dm components
  integer, allocatable :: lo(:), hi(:)

  ! quantities will be allocated with (nlevs,dm) components
  real(kind=dp_t), allocatable :: dx(:,:)
  real(kind=dp_t)              :: dt,time,runtime1,runtime2,Dbar_max,dt_diffusive
  integer                      :: n,nlevs,i,dm,istep,ng_s,init_step,n_Dbar
  type(box)                    :: bx
  type(ml_boxarray)            :: mba
  type(ml_layout)              :: mla
  type(bc_tower)               :: the_bc_tower
  logical, allocatable         :: pmask(:)
  
  ! will be allocated on nlevels
  type(multifab), allocatable  :: rho_old(:)
  type(multifab), allocatable  :: rhotot_old(:)
  type(multifab), allocatable  :: rho_new(:)
  type(multifab), allocatable  :: rhotot_new(:)
  type(multifab), allocatable  :: Temp(:)
  type(multifab), allocatable  :: Temp_ed(:,:)
  type(multifab), allocatable  :: diff_mass_fluxdiv(:)
  type(multifab), allocatable  :: stoch_mass_fluxdiv(:)
  type(multifab), allocatable  :: stoch_mass_flux(:,:)
  type(multifab), allocatable  :: umac(:,:)
  type(multifab), allocatable  :: mtemp(:,:)
  type(multifab), allocatable  :: rhotot_fc(:,:)
  type(multifab), allocatable  :: gradp_baro(:,:)
  type(multifab), allocatable  :: pi(:)
  type(multifab), allocatable  :: eta(:)
  type(multifab), allocatable  :: eta_ed(:,:)
  type(multifab), allocatable  :: kappa(:)

  real(kind=dp_t)              :: total_charge
  type(multifab), allocatable  :: Epot_mass_fluxdiv(:)
  type(multifab), allocatable  :: charge_old(:)
  type(multifab), allocatable  :: charge_new(:)
  type(multifab), allocatable  :: permittivity(:)
  type(multifab), allocatable  :: grad_Epot_old(:,:)
  type(multifab), allocatable  :: grad_Epot_new(:,:)
  type(multifab), allocatable  :: Epot(:)
  type(multifab), allocatable  :: gradPhiApprox(:,:)

  type(multifab), allocatable  :: chem_rate(:)

  ! For HydroGrid
  integer :: narg, farg, un, n_rngs
  character(len=128) :: fname
  logical :: lexist
  logical :: nodal_temp(3)

  ! energy-specific data
  type(multifab), allocatable  :: rhoh_old(:)
  type(multifab), allocatable  :: rhoh_new(:)
  real(kind=dp_t) :: p0_old, p0_new

  !==============================================================
  ! Initialization
  !==============================================================

  call probin_common_init()
  call probin_multispecies_init() 
  call probin_gmres_init()
  call probin_charged_init() 
  call probin_chemistry_init()
  call probin_energy_init()
  call eos_model_init()

   if (use_bl_rng) then
      ! Build the random number engine and give initial distributions for the
      ! F_BaseLib/bl_random RNG module
      call rng_init()
   else
      ! Initialize random numbers *after* the global (root) seed has been set:
      ! This is for the RNG module that sits in Hydrogrid
      call SeedParallelRNG(seed)
   end if

  ! in this example we fix nlevs to be 1
  ! for adaptive simulations where the grids change, cells at finer
  ! resolution don't necessarily exist depending on your tagging criteria,
  ! so max_levs isn't necessary equal to nlevs
  nlevs = 1

  ! dimensionality is set in inputs file
  dm = dim_in
 
  ! now that we have dm, we can allocate these
  allocate(lo(dm),hi(dm),pmask(dm))
  allocate(rho_old(nlevs),rhotot_old(nlevs),pi(nlevs))
  allocate(rho_new(nlevs),rhotot_new(nlevs))
  allocate(Temp(nlevs))
  allocate(diff_mass_fluxdiv(nlevs),stoch_mass_fluxdiv(nlevs))
  allocate(stoch_mass_flux(nlevs,dm))
  allocate(umac(nlevs,dm),mtemp(nlevs,dm),rhotot_fc(nlevs,dm),gradp_baro(nlevs,dm))
  allocate(eta(nlevs),kappa(nlevs))

  ! 1 component in 2D, 3 components in 3D
  allocate(eta_ed(nlevs,2*dm-3))
  allocate(Temp_ed(nlevs,2*dm-3))

  allocate(Epot_mass_fluxdiv(nlevs))
  allocate(charge_old(nlevs),charge_new(nlevs))
  allocate(permittivity(nlevs))
  allocate(grad_Epot_old(nlevs,dm),grad_Epot_new(nlevs,dm))
  allocate(Epot(nlevs))
  allocate(gradPhiApprox(nlevs,dm))

  allocate(chem_rate(nlevs))

  ! energy specific
  allocate(rhoh_old(nlevs),rhoh_new(nlevs))

  ! build pmask
  pmask = .false.
  do i=1,dm
     if (bc_lo(i) .eq. PERIODIC .and. bc_hi(i) .eq. PERIODIC) then
        pmask(i) = .true.
     end if
  end do

  if (advection_type .eq. 0) then
     ng_s = 2 ! centered advection
  else if (advection_type .le. 3) then
     ng_s = 3 ! bilinear bds or unlimited quadratic bds
  else if (advection_type .eq. 4) then
     ng_s = 4 ! limited quadratic bds
  end if

  if (restart .ge. 0) then

     init_step = restart + 1
     call bl_error('Error: restart function not supported for energy code')

  else

     init_step = 1
     time = start_time
     
     ! tell mba how many levels and dimensionality of problem
     call ml_boxarray_build_n(mba,nlevs,dm)

     ! tell mba about the ref_ratio between levels
     ! mba%rr(n-1,i) is the refinement ratio between levels n-1 and n in direction i
     ! we use refinement ratio of 2 in every direction between all levels
     do n=2,nlevs
        mba%rr(n-1,:) = 2
     end do

     ! create a box from (0,0) to (n_cells-1,n_cells-1)
     lo(1:dm) = 0
     hi(1:dm) = n_cells(1:dm)-1
     bx = make_box(lo,hi)

     ! tell mba about the problem domain at every level
     mba%pd(1) = bx
     do n=2,nlevs
        mba%pd(n) = refine(mba%pd(n-1),mba%rr((n-1),:))
     end do

     ! initialize the boxarray at level 1 to be one single box
     call boxarray_build_bx(mba%bas(1),bx)

     ! overwrite the boxarray at level 1 to respect max_grid_size
     call boxarray_maxsize(mba%bas(1),max_grid_size)

     ! now build the boxarray at other levels
     if (nlevs .ge. 2) then
        call bl_error("Need to build boxarray for n>1")
     end if

     ! build the ml_layout, mla
     call ml_layout_build(mla,mba,pmask)

     ! don't need this anymore - free up memory
     call destroy(mba)

     do n=1,nlevs
        call multifab_build(rho_old(n)   ,mla%la(n),nspecies,ng_s)
        call multifab_build(rhotot_old(n),mla%la(n),1       ,ng_s)
        call multifab_build(rhoh_old(n)  ,mla%la(n),1       ,ng_s)
        call multifab_build(Temp(n)      ,mla%la(n),1       ,ng_s)
        call multifab_build(pi(n)        ,mla%la(n),1       ,1)
        do i=1,dm
           call multifab_build_edge(umac(n,i),mla%la(n),1,1,i)
        end do
     end do

     if (use_charged_fluid) then
        do n=1,nlevs        
           call multifab_build(Epot_mass_fluxdiv(n), mla%la(n),nspecies,0) 
        end do
     end if

  end if

  ! data structures to help with reservoirs
  call build_bc_multifabs(mla)

  ! set grid spacing at each level
  allocate(dx(nlevs,MAX_SPACEDIM))
  dx(1,1:MAX_SPACEDIM) = (prob_hi(1:MAX_SPACEDIM)-prob_lo(1:MAX_SPACEDIM)) &
       / n_cells(1:MAX_SPACEDIM)
  ! check that the grid spacing is the same in each direction for the first dm dimensions
  select case (dm) 
  case(2)
     if (dx(1,1) .ne. dx(1,2)) then
        call bl_error('ERROR: main_driver.f90, in 2D we only support dx=dy')
     end if
  case(3)
     if ((dx(1,1) .ne. dx(1,2)) .or. (dx(1,1) .ne. dx(1,3))) then
        call bl_error('ERROR: main_driver.f90, in 3D we only support dx=dy=dz')
     end if
  case default
     call bl_error('ERROR: main_driver.f90, dimension should be only equal to 2 or 3')
  end select

  ! use refined dx for next level
  ! assume refinement ratio is the same in each direction
  ! we do this because dx is allocated over MAX_SPACEDIM dimensions,
  ! whereas mba is only allocated over dm dimensions
  do n=2,nlevs
     dx(n,:) = dx(n-1,:) / mba%rr(n-1,1)
  end do

  !=======================================================
  ! Setup boundary condition bc_tower
  !=======================================================
 
  ! bc_tower structure in memory
  ! 1:dm = velocity
  ! dm+1 = pressure
  !
  ! next, for scalars, scal_bc_comp=dm+2
  ! there are 2*nspecies+2 "scalars"
  ! scal_bc_comp = rhotot
  ! scal_bc_comp+1 = c_i
  ! scal_bc_comp+nspecies+1 = molfrac or massfrac (dimensionless fractions)
  ! scal_bc_comp+2*nspecies+1 = temp_bc_comp = temperature
  ! scal_bc_comp+2*nspecies+2 = Epot_bc_comp = electric potential
  ! scal_bc_comp+2*nspecies+3 = h_bc_comp = enthalpy
  !
  ! next, for transport coefficients, tran_bc_comp = scal_bc_comp+2*nspecies+4
  ! we say there is "one" transport coefficient
  ! It may be better if each transport coefficient has its own BC code?
  ! I think the only place this is used is average_cc_to_node/face/edge
  ! I cannot right now foresee a case where different values would be used in different places
  ! so it is OK to keep num_tran_bc_in=1. But note the same code applies to eta,kappa and chi's
  ! tran_bc_comp = diffusion coefficients (eta,kappa,chi)
  call initialize_bc(the_bc_tower,nlevs,dm,mla%pmask, &
                     num_scal_bc_in=2*nspecies+4, &
                     num_tran_bc_in=1)

  do n=1,nlevs
     ! define level n of the_bc_tower
     call bc_tower_level_build(the_bc_tower,n,mla%la(n))
  end do

  ! build layouts for staggered multigrid solver and macproject within preconditioner
  call stag_mg_layout_build(mla)
  call mgt_macproj_precon_build(mla,dx,the_bc_tower)

  if (restart .lt. 0) then

     ! initialize rho and umac in valid region only
     ! also initialize rhoh, Temp, p0
     call init_lowmach(mla,umac,rhotot_old,rho_old,rhoh_old,Temp,p0_old)

     ! initialize pi, including ghost cells
     do n=1,nlevs
        call multifab_setval(pi(n),0.d0,all=.true.)
     end do

  end if

  ! fill temperature ghost cells
  call fill_Temp_ghost_cells(mla,Temp,dx,the_bc_tower)

  ! compute rhotot from rho in VALID REGION
  call compute_rhotot(mla,rho_old,rhotot_old)

  ! fill rho and rhotot ghost cells
  ! FIXME - different EOS; won't work for physical boundaries
  ! we need to write a different version that computes rhotot from the EOS
  ! which means this needs p0, Temp, and c
  call fill_rho_rhotot_ghost(mla,rho_old,rhotot_old,dx,the_bc_tower)

  ! fill rhoh ghost cells - note that rhotot needs filled ghost cells
  call fill_rhoh_ghost_cells(mla,rhoh_old,rhotot_old,dx,the_bc_tower)

  !=======================================================
  ! Build multifabs for all the variables
  !=======================================================

  do n=1,nlevs 
     call multifab_build(rho_new(n)          ,mla%la(n),nspecies,ng_s)
     call multifab_build(rhotot_new(n)       ,mla%la(n),1       ,ng_s) 
     call multifab_build(rhoth_new(n)        ,mla%la(n),1       ,ng_s) 
     call multifab_build(eta(n)              ,mla%la(n),1       ,1)
     call multifab_build(kappa(n)            ,mla%la(n),1       ,1)
     call multifab_build(diff_mass_fluxdiv(n),mla%la(n),nspecies,0) 
     do i=1,dm
        call multifab_build_edge     (mtemp(n,i),mla%la(n),1,0,i)
        call multifab_build_edge( rhotot_fc(n,i),mla%la(n),1,0,i)
        call multifab_build_edge(gradp_baro(n,i),mla%la(n),1,0,i)
     end do
  end do

  do n=1,nlevs
     ! eta and Temp on nodes (2d) or edges (3d)
     if (dm .eq. 2) then
        call multifab_build_nodal(eta_ed(n,1),mla%la(n),1,0)
        call multifab_build_nodal(Temp_ed(n,1),mla%la(n),1,0)
     else
        nodal_temp(1) = .true.
        nodal_temp(2) = .true.
        nodal_temp(3) = .false.
        call multifab_build(eta_ed(n,1),mla%la(n),1,0,nodal_temp)
        call multifab_build(Temp_ed(n,1),mla%la(n),1,0,nodal_temp)
        nodal_temp(1) = .true.
        nodal_temp(2) = .false.
        nodal_temp(3) = .true.
        call multifab_build(eta_ed(n,2),mla%la(n),1,0,nodal_temp)
        call multifab_build(Temp_ed(n,2),mla%la(n),1,0,nodal_temp)
        nodal_temp(1) = .false.
        nodal_temp(2) = .true.
        nodal_temp(3) = .true.
        call multifab_build(eta_ed(n,3),mla%la(n),1,0,nodal_temp)
        call multifab_build(Temp_ed(n,3),mla%la(n),1,0,nodal_temp)
     end if
  end do

  if (variance_coef_mass .ne. 0.d0) then
     do n=1,nlevs
        call multifab_build(stoch_mass_fluxdiv(n),mla%la(n),nspecies,0) 
        do i=1,dm
           call multifab_build_edge(stoch_mass_flux(n,i),mla%la(n),nspecies,0,i)
        end do
     end do
  end if

  if (use_charged_fluid) then
     do n=1,nlevs
        call multifab_build(charge_old(n)  ,mla%la(n),1,1)
        call multifab_build(charge_new(n)  ,mla%la(n),1,1)
        call multifab_build(permittivity(n),mla%la(n),1,1)
        call multifab_build(Epot(n)        ,mla%la(n),1,1)
        do i=1,dm
           call multifab_build_edge(grad_Epot_old(n,i),mla%la(n),1,1,i)
           call multifab_build_edge(grad_Epot_new(n,i),mla%la(n),1,1,i)
           call multifab_build_edge(gradPhiApprox(n,i),mla%la(n),1,0,i)
        end do
     end do
  end if

  if (nreactions .gt. 0) then
     do n=1,nlevs
        call multifab_build(chem_rate(n),mla%la(n),nspecies,0)
     end do
  end if

  ! allocate and build multifabs that will contain random numbers
  if (algorithm_type .eq. 2 .or. algorithm_type .eq. 5 ) then
     n_rngs = 2
  else
     n_rngs = 1
  end if
  call init_mass_stochastic(mla,n_rngs)
  call init_m_stochastic(mla,n_rngs)

  if (use_bl_rng) then
     ! save random state for writing checkpoint
     call bl_rng_copy_engine(rng_eng_diffusion_chk,rng_eng_diffusion)
  end if

  !=====================================================================
  ! Initialize values
  !=====================================================================

  if (use_charged_fluid) then
     call bl_error("use_charged_fluid does not work in energy code")
  end if

  if (barodiffusion_type .gt. 0) then

     ! this computes an initial guess at p using HSE
     call compute_HSE_pres(mla,rhotot_old,pi,dx,the_bc_tower)

     ! compute grad p for barodiffusion
     call compute_grad(mla,pi,gradp_baro,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

  end if

  ! initialize eta and kappa
  call compute_eta_kappa(mla,eta,eta_ed,kappa,rho_old,rhotot_old,Temp,dx, &
                         the_bc_tower%bc_tower_array)

  ! now that we have eta, we can initialize the inhomogeneous velocity bcs
  ! set inhomogeneous bc conditions
  call set_inhomogeneous_vel_bcs(mla,vel_bc_n,vel_bc_t,eta_ed,dx,time, &
                                 the_bc_tower%bc_tower_array)

  do n=1,nlevs
     do i=1,dm
        ! set normal velocity on physical domain boundaries
        call multifab_physbc_domainvel(umac(n,i),vel_bc_comp+i-1, &
                                       the_bc_tower%bc_tower_array(n), &
                                       dx(n,:),vel_bc_n(n,:))
        ! set transverse velocity behind physical boundaries
        call multifab_physbc_macvel(umac(n,i),vel_bc_comp+i-1, &
                                    the_bc_tower%bc_tower_array(n), &
                                    dx(n,:),vel_bc_t(n,:))
        ! fill periodic and interior ghost cells
        call multifab_fill_boundary(umac(n,i))
     end do
  end do

  if (restart .lt. 0) then

     if (print_int .gt. 0) then
        if (parallel_IOProcessor()) write(*,*) "Initial state:"  
        call sum_mass(rho_old, 0) ! print out the total mass to check conservation
        ! compute rhotot on faces
        call average_cc_to_face(nlevs,rhotot_old,rhotot_fc,1,scal_bc_comp,1, &
                                the_bc_tower%bc_tower_array)
        ! compute momentum
        call convert_m_to_umac(mla,rhotot_fc,mtemp,umac,.false.)
        call sum_momenta(mla,mtemp)
        call eos_check(mla,rho_old)
     end if

     ! add initial momentum fluctuations
     ! do not call for overdamped codes since the steady Stokes solver will 
     ! wipe out the initial condition to solver tolerance
     if (algorithm_type .ne. 2 .and. &
         variance_coef_mom .ne. 0.d0 .and. &
         initial_variance .ne. 0.d0) then
        call add_m_fluctuations(mla,dx,initial_variance*variance_coef_mom, &
                                umac,rhotot_old,Temp,the_bc_tower)

        if (print_int .gt. 0) then
           if (parallel_IOProcessor()) write(*,*) "After adding momentum fluctuations:"
           call sum_mass(rho_old, 0) ! print out the total mass to check conservation
           ! compute rhotot on faces
           call average_cc_to_face(nlevs,rhotot_old,rhotot_fc,1,scal_bc_comp,1, &
                                   the_bc_tower%bc_tower_array)
           ! compute momentum
           call convert_m_to_umac(mla,rhotot_fc,mtemp,umac,.false.)
           call sum_momenta(mla,mtemp)
           call eos_check(mla,rho_old)
        end if

     end if

     if (fixed_dt .gt. 0.d0) then
        dt = fixed_dt
     else
        call estdt(mla,umac,dx,dt)
        n_Dbar = nspecies*(nspecies-1)/2
        Dbar_max = maxval(Dbar(1:n_Dbar))
        dt_diffusive = cfl*dx(1,1)**2/(2*dm*Dbar_max)
        dt = min(dt,dt_diffusive)
     end if
     
  end if

  !=====================================================================
  ! Initialize HydroGrid for analysis
  !=====================================================================
  if((hydro_grid_int>0) .or. (stats_int>0)) then
     narg = command_argument_count()
     farg = 1
     if (narg >= 1) then
        call get_command_argument(farg, value = fname)
        inquire(file = fname, exist = lexist )
        if ( lexist ) then
           un = unit_new()
           open(unit=un, file = fname, status = 'old', action = 'read')
           
           ! We will also pass temperature here but no additional scalars
           call initialize_hydro_grid(mla,rho_old,dt,dx,namelist_file=un, &
                                      nspecies_in=nspecies, &
                                      nscal_in=0, &
                                      exclude_last_species_in=.false., &
                                      analyze_velocity=.true., &
                                      analyze_density=.true., &
                                      analyze_temperature=.true.) 
           
           close(unit=un)
        end if
     end if
  end if

  ! this routine is only called for all inertial simulations (both restart and non-restart)
  ! it does the following:
  ! 1. fill mass random numbers
  ! 2. computes mass fluxes and flux divergences
  ! if restarting, the subroutine ends; otherwise
  ! 3. perform an initial projection
  !
  ! overdamped schemes need to do 1. and 2. within the advance_timestep routine
  ! in principle, performing an initial projection for overdamped will change
  ! the reference state for the GMRES solver
  ! For overdamped the first ever solve cannot have a good reference state
  ! so in general there is the danger it will be less accurate than subsequent solves
  ! but I do not see how one can avoid that
  ! From this perspective it may be useful to keep initial_projection even in overdamped
  ! because different gmres tolerances may be needed in the first step than in the rest
  if (algorithm_type .ne. 2) then
     call initial_projection(mla,umac,rho_old,rhotot_old,gradp_baro, &
                             Epot_mass_fluxdiv,diff_mass_fluxdiv, &
                             stoch_mass_fluxdiv,stoch_mass_flux,chem_rate, &
                             Temp,eta,eta_ed,dt,dx,the_bc_tower, &
                             charge_old,grad_Epot_old,Epot,permittivity)
  end if

  if (restart .lt. 0) then

     !=====================================================================
     ! Process initial conditions (non-restart runs)
     !=====================================================================

     if (print_int .gt. 0) then
        if (parallel_IOProcessor()) write(*,*) "After initial projection:"  
        call sum_mass(rho_old,0) ! print out the total mass to check conservation
        ! compute rhotot on faces
        call average_cc_to_face(nlevs,rhotot_old,rhotot_fc,1,scal_bc_comp,1, &
                                the_bc_tower%bc_tower_array)
        ! compute momentum
        call convert_m_to_umac(mla,rhotot_fc,mtemp,umac,.false.)
        call sum_momenta(mla,mtemp)
        ! FIXME - need to write a custom eos_check using compute_p() and subtracting p0
     end if

     !=====================================================================
     ! Hydrogrid analysis and output for initial data
     !=====================================================================

     ! write initial plotfile
     if (plot_int .gt. 0) then
        if (parallel_IOProcessor()) then
           write(*,*), 'writing initial plotfile 0'
        end if
        call write_plotfile(mla,"plt",rho_old,rhotot_old,rhoh_old,Temp, &
                            umac,pi,p0_old,0,dx,time)
     end if

     ! write initial checkpoint
     if (chk_int .gt. 0) then
        call bl_error('Error: checkpoint not supported for energy code')
     end if
     
     if (stats_int .gt. 0) then
        ! write initial vertical and horizontal averages (hstat and vstat files)   
        call print_stats(mla,dx,0,time,umac=umac,rho=rho_old,temperature=Temp)
     end if

     ! We do the analysis first so we include the initial condition in the files if n_steps_skip=0
     if (n_steps_skip .eq. 0) then

        ! Add this snapshot to the average in HydroGrid
        if (hydro_grid_int > 0) then
           call analyze_hydro_grid(mla,dt,dx,istep,umac=umac,rho=rho_old,temperature=Temp)
        end if

        if (hydro_grid_int > 0 .and. n_steps_save_stats > 0) then
           call save_hydro_grid(id=0, step=0)
        end if

     end if    

  end if

  !=======================================================
  ! Begin time stepping loop
  !=======================================================

  do istep=init_step,max_step

     if (fixed_dt .le. 0.d0) then
        call estdt(mla,umac,dx,dt)
        n_Dbar = nspecies*(nspecies-1)/2
        Dbar_max = maxval(Dbar(1:n_Dbar))
        dt_diffusive = cfl*dx(1,1)**2/(2*dm*Dbar_max)
        dt = min(dt,dt_diffusive)
     end if

     if (parallel_IOProcessor()) then
        if ( (print_int .gt. 0 .and. mod(istep,print_int) .eq. 0) ) &
             print*,"Begin Advance; istep =",istep,"dt =",dt,"time =",time
     end if

     runtime1 = parallel_wtime()

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! advance the solution by dt
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! notes: eta, eta_ed, and kappa could be built and initialized within the advance routines
      ! but for now we pass them around (it does save a few flops)
      ! diff/stoch_mass_fluxdiv could be built locally within the overdamped
      ! routine, but since we have them around anyway for inertial we pass them in
     if (algorithm_type .eq. 1) then
        call advance_timestep_inertial(mla,umac,rho_old,rho_new, &
                                       rhotot_old,rhotot_new,rhoh_old,rhoh_new, &
                                       p0_old,p0_new,gradp_baro,Temp, &
                                       pi,dx,dt,the_bc_tower)
     else
        call bl_error("Error: invalid algorithm_type")
     end if

     time = time + dt

     if (parallel_IOProcessor()) then
        if ( (print_int .gt. 0 .and. mod(istep,print_int) .eq. 0) ) &
             print*,"End Advance; istep =",istep,"DT =",dt,"TIME =",time
     end if

      runtime2 = parallel_wtime() - runtime1
      call parallel_reduce(runtime1, runtime2, MPI_MAX, proc=parallel_IOProcessorNode())
      if (parallel_IOProcessor()) then
        if ( (print_int .gt. 0 .and. mod(istep,print_int) .eq. 0) ) &
           print*,'Time to advance timestep: ',runtime1,' seconds'
      end if
      
      if ( (print_int .gt. 0 .and. mod(istep,print_int) .eq. 0) &
          .or. &
          (istep .eq. max_step) ) then
          if (parallel_IOProcessor()) write(*,*) "At time step ", istep, " t=", time           
          call sum_mass(rho_new, istep) ! print out the total mass to check conservation
          ! compute rhotot on faces
          call average_cc_to_face(nlevs,rhotot_new,rhotot_fc,1,scal_bc_comp,1, &
                                  the_bc_tower%bc_tower_array)
          ! compute momentum
          call convert_m_to_umac(mla,rhotot_fc,mtemp,umac,.false.)
          call sum_momenta(mla,mtemp)
      end if

      ! We do the analysis first so we include the initial condition in the files if n_steps_skip=0
      if (istep >= n_steps_skip) then

         ! write plotfile at specific intervals
         if ((plot_int.gt.0 .and. mod(istep,plot_int).eq.0) .or. (istep.eq.max_step)) then
            if (parallel_IOProcessor()) then
               write(*,*), 'writing plotfiles at timestep =', istep 
            end if
            call write_plotfile(mla,"plt",rho_new,rhotot_new,rhoh_new,Temp, &
                                umac,pi,p0_new,istep,dx,time)
         end if

         ! write checkpoint at specific intervals
         if ((chk_int.gt.0 .and. mod(istep,chk_int).eq.0)) then
            call bl_error('Error: checkpoint function not supported for energy code')
         end if

         ! print out projection (average) and variance
         if ( (stats_int > 0) .and. &
               (mod(istep,stats_int) .eq. 0) ) then
            ! Compute vertical and horizontal averages (hstat and vstat files)   
            call print_stats(mla,dx,istep,time,umac=umac,rho=rho_new,temperature=Temp)    
         end if

         ! Add this snapshot to the average in HydroGrid
         if ( (hydro_grid_int > 0) .and. &
              ( mod(istep,hydro_grid_int) .eq. 0 ) ) then
            call analyze_hydro_grid(mla,dt,dx,istep,umac=umac,rho=rho_new,temperature=Temp)
         end if

         if ( (hydro_grid_int > 0) .and. &
              (n_steps_save_stats > 0) .and. &
              ( mod(istep,n_steps_save_stats) .eq. 0 ) ) then
              call save_hydro_grid(id=istep/n_steps_save_stats, step=istep)            
         end if

      end if

      ! set old state to new state

      p0_old = p0_new
      do n=1,nlevs
         call multifab_copy_c(rho_old(n)   ,1,rho_new(n)   ,1,nspecies,rho_old(n)%ng)
         call multifab_copy_c(rhotot_old(n),1,rhotot_new(n),1,1       ,rhotot_old(n)%ng)
         call multifab_copy_c(rhoh_old(n)  ,1,rhoh_new(n)  ,1,1       ,rhoh_old(n)%ng)
      end do

  end do

  !=======================================================
  ! Destroy multifabs and layouts
  !=======================================================

  if((hydro_grid_int>0) .or. (stats_int>0)) then
     call finalize_hydro_grid()
  end if

  call destroy_bc_multifabs(mla)
  call destroy_mass_stochastic(mla)
  call destroy_m_stochastic(mla)

  do n=1,nlevs
     call multifab_destroy(rho_old(n))
     call multifab_destroy(rhotot_old(n))
     call multifab_destroy(rhoh_old(n))
     call multifab_destroy(pi(n))
     do i=1,dm
        call multifab_destroy(umac(n,i))
     end do
  end do

  if (use_charged_fluid) then
     do n=1,nlevs
        call multifab_destroy(Epot_mass_fluxdiv(n))
     end do
  end if

  do n=1,nlevs
     call multifab_destroy(rho_new(n))
     call multifab_destroy(rhotot_new(n))
     call multifab_destroy(rhoh_new(n))
     call multifab_destroy(Temp(n))
     call multifab_destroy(eta(n))
     call multifab_destroy(kappa(n))
     call multifab_destroy(diff_mass_fluxdiv(n))
     do i=1,dm
        call multifab_destroy(mtemp(n,i))
        call multifab_destroy(rhotot_fc(n,i))
        call multifab_destroy(gradp_baro(n,i))
     end do
     do i=1,size(eta_ed,dim=2)
        call multifab_destroy(eta_ed(n,i))
        call multifab_destroy(Temp_ed(n,i))
     end do
  end do

  if (variance_coef_mass .ne. 0.d0) then
     do n=1,nlevs
        call multifab_destroy(stoch_mass_fluxdiv(n))
        do i=1,dm
           call multifab_destroy(stoch_mass_flux(n,i))
        end do
     end do
  end if

  if (use_charged_fluid) then
     do n=1,nlevs
        call multifab_destroy(charge_old(n))
        call multifab_destroy(charge_new(n))
        call multifab_destroy(permittivity(n))
        call multifab_destroy(Epot(n))
        do i=1,dm
           call multifab_destroy(grad_Epot_old(n,i))
           call multifab_destroy(grad_Epot_new(n,i))
           call multifab_destroy(gradPhiApprox(n,i))
        end do
     end do
  end if
  
  if (nreactions .gt. 0) then
     do n=1,nlevs
        call multifab_destroy(chem_rate(n))
     end do
  end if

  deallocate(lo,hi,pmask)
  deallocate(rho_old,rhotot_old,pi)
  deallocate(rho_new,rhotot_new)
  deallocate(Temp)
  deallocate(diff_mass_fluxdiv,stoch_mass_fluxdiv)
  deallocate(stoch_mass_flux)
  deallocate(umac,mtemp,rhotot_fc,gradp_baro)
  deallocate(eta,kappa)
  deallocate(eta_ed)
  deallocate(Temp_ed)

  deallocate(Epot_mass_fluxdiv)
  deallocate(charge_old,charge_new)
  deallocate(permittivity)
  deallocate(grad_Epot_old,grad_Epot_new)
  deallocate(Epot)
  deallocate(gradPhiApprox)

  deallocate(chem_rate)

  deallocate(rhoh_old,rhoh_new)

  deallocate(dx)

  if (use_bl_rng) then
     call rng_destroy()
  end if

  call stag_mg_layout_destroy()
  call mgt_macproj_precon_destroy()
  call destroy(mla)
  call bc_tower_destroy(the_bc_tower)

end subroutine main_driver
