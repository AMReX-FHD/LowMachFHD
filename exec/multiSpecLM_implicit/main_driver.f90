subroutine main_driver()

  use boxlib
  use multifab_module
  use bl_IO_module
  use layout_module
  use init_module
  use compute_mixture_properties_module
  use initial_projection_module
  use write_plotfileLM_module
  use advance_diffusion_module
  use advance_timestep_overdamped_module
  use define_bc_module
  use bc_module
  use analysis_module
  use analyze_spectra_module
  use convert_m_to_umac_module
  use eos_check_module
  use stochastic_mass_fluxdiv_module
  use stochastic_m_fluxdiv_module
  use ParallelRNGs 
  use mass_flux_utilities_module
  use convert_stag_module
  use probin_common_module, only: prob_lo, prob_hi, n_cells, dim_in, hydro_grid_int, &
                                  max_grid_size, n_steps_save_stats, n_steps_skip, &
                                  plot_int, seed, stats_int, bc_lo, bc_hi, probin_common_init, &
                                  advection_type, fixed_dt, max_step, &
                                  algorithm_type, variance_coef_mom, initial_variance
  use probin_multispecies_module, only: nspecies, mol_frac_bc_comp, &
                                        rho_part_bc_comp, start_time, temp_bc_comp, &
                                        probin_multispecies_init
  use probin_gmres_module, only: probin_gmres_init

  implicit none

  ! quantities will be allocated with dm components
  integer, allocatable :: lo(:), hi(:)

  ! quantities will be allocated with (nlevs,dm) components
  real(kind=dp_t), allocatable :: dx(:,:)
  real(kind=dp_t)              :: dt,time
  integer                      :: n,nlevs,i,dm,istep,ng_s
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
  type(multifab), allocatable  :: umac(:,:)
  type(multifab), allocatable  :: pres(:)
  type(multifab), allocatable  :: eta(:)
  type(multifab), allocatable  :: eta_ed(:,:)
  type(multifab), allocatable  :: kappa(:)

  ! For HydroGrid
  integer :: narg, farg, un, n_rngs
  character(len=128) :: fname
  logical :: lexist
  logical :: nodal_temp(3)
  
  !==============================================================
  ! Initialization
  !==============================================================

  call probin_common_init()
  call probin_multispecies_init() 
  call probin_gmres_init()
  
  ! Initialize random numbers *after* the global (root) seed has been set:
  call SeedParallelRNG(seed)

  ! in this example we fix nlevs to be 1
  ! for adaptive simulations where the grids change, cells at finer
  ! resolution don't necessarily exist depending on your tagging criteria,
  ! so max_levs isn't necessary equal to nlevs
  nlevs = 1

  ! dimensionality is set in inputs file
  dm = dim_in
 
  ! now that we have dm, we can allocate these
  allocate(lo(dm),hi(dm))
  allocate(dx(nlevs,dm))
  allocate(rho_old(nlevs),rhotot_old(nlevs))
  allocate(rho_new(nlevs),rhotot_new(nlevs))
  allocate(Temp(nlevs),diff_mass_fluxdiv(nlevs),stoch_mass_fluxdiv(nlevs))
  allocate(umac(nlevs,dm),pres(nlevs))
  allocate(eta(nlevs),kappa(nlevs))
  if (dm .eq. 2) then
     allocate(eta_ed(nlevs,1))
     allocate(Temp_ed(nlevs,1))
  else if (dm .eq. 3) then
     allocate(eta_ed(nlevs,3))
     allocate(Temp_ed(nlevs,3))
  end if

  !==============================================================
  ! Setup parallelization: Create boxes and layouts for multifabs
  !==============================================================
  
  ! tell mba how many levels and dimensionality of problem
  call ml_boxarray_build_n(mba,nlevs,dm)

  ! tell mba about the ref_ratio between levels
  ! mba%rr(n-1,i) is the refinement ratio between levels n-1 and n in direction i
  ! we use refinement ratio of 2 in every direction between all levels
  do n=2,nlevs
     mba%rr(n-1,:) = 2
  end do

  ! set grid spacing at each level; presently grid spacing is same in all direction
  dx(1,1:dm) = (prob_hi(1:dm)-prob_lo(1:dm)) / n_cells(1:dm)
  
  ! check whether dimensionality & grid spacing passed correct 
  select case (dm) 
    case(2)
      if (dx(1,1) .ne. dx(1,2)) then
        call bl_error('ERROR: main_driver.f90, we only support dx=dy')
      end if    
    case(3)
      if ((dx(1,1) .ne. dx(1,2)) .or. (dx(1,1) .ne. dx(1,3))) then
        call bl_error('ERROR: main_driver.f90, we only support dx=dy=dz')
      end if    
    case default
      call bl_error('ERROR: main_driver.f90, dimension should be only equal to 2 or 3')
  end select

  ! use refined dx for next level
  do n=2,nlevs
     dx(n,:) = dx(n-1,:) / mba%rr(n-1,:)
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

  ! build pmask
  allocate(pmask(dm))
  pmask = .false.
  do i=1,dm
     if (bc_lo(i) .eq. PERIODIC .and. bc_hi(i) .eq. PERIODIC) then
        pmask(i) = .true.
     end if
  end do

  ! build the ml_layout, mla
  call ml_layout_build(mla,mba,pmask)
  deallocate(pmask)

  ! don't need this anymore - free up memory
  call destroy(mba)

  !=======================================================
  ! Setup boundary condition bc_tower
  !=======================================================
 
  ! bc_tower structure in memory
  ! 1:dm = velocity
  ! dm+1 = pressure
  ! dm+2 = scal_bc_comp = rhotot
  ! scal_bc_comp+1 = rho_i
  ! scal_bc_comp+nspecies+1 = molfrac or massfrac (dimensionless fractions)
  ! scal_bc_comp+2*nspecies+1 = temp_bc_comp = temperature
  ! scal_bc_comp+2*nspecies+2 = tran_bc_comp = diffusion coefficients (eta,kappa,chi)
  ! It may be better if each transport coefficient has its own BC code?
  ! I think the only place this is used is average_cc_to_node/face/edge
  ! I cannot right now foresee a case where different values would be used in different places
  ! so it is OK to keep num_tran_bc_in=1. But note the same code applies to eta,kappa and chi's
  call initialize_bc(the_bc_tower,nlevs,dm,mla%pmask, &
                     num_scal_bc_in=2*nspecies+2,num_tran_bc_in=1)

  do n=1,nlevs
     ! define level n of the_bc_tower
     call bc_tower_level_build(the_bc_tower,n,mla%la(n))
  end do

  ! these quantities are populated here and defined in probin_multispecies 
  rho_part_bc_comp   = scal_bc_comp + 1
  mol_frac_bc_comp   = scal_bc_comp + nspecies + 1
  temp_bc_comp       = scal_bc_comp + 2*nspecies + 1

  !=======================================================
  ! Build multifabs for all the variables
  !=======================================================

  if (advection_type .eq. 0) then
     ng_s = 1 ! centered advection
  else
     ng_s = 3 ! bds advection
  end if

  ! build multifab with nspecies component and one ghost cell
  do n=1,nlevs

     call multifab_build(rho_old(n),           mla%la(n),nspecies,ng_s)
     call multifab_build(rhotot_old(n),        mla%la(n),1,       ng_s) 
     call multifab_build(rho_new(n),           mla%la(n),nspecies,ng_s)
     call multifab_build(rhotot_new(n),        mla%la(n),1,       ng_s) 
     call multifab_build(Temp(n),              mla%la(n),1,       ng_s)
     call multifab_build(diff_mass_fluxdiv(n), mla%la(n),nspecies,0) 
     call multifab_build(stoch_mass_fluxdiv(n),mla%la(n),nspecies,0) 
     ! pressure - need 1 ghost cell since we calculate its gradient
     call multifab_build(pres(n),mla%la(n),1,1)
     call multifab_build(eta(n)  ,mla%la(n),1,1)
     call multifab_build(kappa(n),mla%la(n),1,1)
     do i=1,dm
        call multifab_build_edge(umac(n,i),mla%la(n),1,1,i)
     end do

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

  ! allocate and build multifabs that will contain random numbers
  if (algorithm_type .eq. 0 .or. algorithm_type .eq. 1) then
     n_rngs = 1
  else if (algorithm_type .eq. 2) then
     n_rngs = 2
  end if
  call init_mass_stochastic(mla,n_rngs)
  call init_m_stochastic(mla,n_rngs)

  ! fill random flux multifabs with new random numbers
  call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)
  call fill_m_stochastic(mla)

  !=====================================================================
  ! Initialize values
  !=====================================================================

  ! initialize the time 
  time = start_time

  ! initialize rho
  call init_rho(rho_old,dx,time,the_bc_tower%bc_tower_array)
  call eos_check(mla,rho_old)
  call compute_rhotot(mla,rho_old,rhotot_old)

  ! initialize Temp
  call init_Temp(Temp,dx,time,the_bc_tower%bc_tower_array)
  if (dm .eq. 2) then
     call average_cc_to_node(nlevs,Temp,Temp_ed(:,1),1,tran_bc_comp,1,the_bc_tower%bc_tower_array)
  else if (dm .eq. 3) then
     call average_cc_to_edge(nlevs,Temp,Temp_ed,1,tran_bc_comp,1,the_bc_tower%bc_tower_array)
  end if

  ! initialize kappa
  do n=1,nlevs
     call multifab_setval(kappa(n), 1.d0, all=.true.)
  end do

  ! initialize eta
  call compute_eta(mla,eta,eta_ed,rho_old,rhotot_old,Temp,pres,dx,the_bc_tower%bc_tower_array)

  ! initialize pressure and velocity
  do n=1,nlevs
     call multifab_setval(pres(n), 0.d0, all=.true.)
     do i=1,dm
        call multifab_setval(umac(n,i), 0.d0, all=.true.)
     end do
  end do

  ! add initial momentum fluctuations - only call in inertial code for now
  ! Note, for overdamped code, the steady Stokes solver will wipe out the initial condition
  if (algorithm_type .eq. 0 .and. initial_variance .ne. 0.d0) then
     call add_m_fluctuations(mla,dx,initial_variance*variance_coef_mom, &
                             umac,rhotot_old,Temp,the_bc_tower)
     
  end if

  dt = fixed_dt

  !=====================================================================
  ! Initialize HydroGrid for analysis
  !=====================================================================
  if((abs(hydro_grid_int)>0) .or. (stats_int>0)) then
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

  ! initial projection - only truly needed for inertial algorithm
  ! for the overdamped algorithm, this only changes the reference state for the first
  ! gmres solve in the first time step
  ! Yes, I think in the purely overdamped version this can be removed
  ! In either case the first ever solve cannot have a good reference state
  ! so in general there is the danger it will be less accurate than subsequent solves
  ! but I do not see how one can avoid that
  ! From this perspective it may be useful to keep initial_projection even in overdamped
  ! because different gmres tolerances may be needed in the first step than in the rest
  if (algorithm_type .eq. 0) then
     call initial_projection(mla,umac,rho_old,rhotot_old,diff_mass_fluxdiv, &
                             stoch_mass_fluxdiv,Temp,dt,dx,n_rngs,the_bc_tower)
  end if

  !=======================================================
  ! Begin time stepping loop
  !=======================================================

  ! free up memory counters 
  istep      = 0

  do while(istep<=max_step)

      if (parallel_IOProcessor()) then
         print*,"Begin Advance; istep =",istep,"dt =",dt,"time =",time
      end if

      ! We do the analysis first so we include the initial condition in the files if n_steps_skip=0
      if (istep >= n_steps_skip) then

         ! print out projection (average) and variance
         if ( (stats_int > 0) .and. &
               (mod(istep-n_steps_skip,stats_int) .eq. 0) ) then
            ! Compute vertical and horizontal averages (hstat and vstat files)   
            call print_stats(mla,dx,istep-n_steps_skip,time,rho=rho_old,temperature=Temp)            
         end if

         ! Add this snapshot to the average in HydroGrid
         if ( (hydro_grid_int > 0) .and. &
              ( mod(istep-n_steps_skip,hydro_grid_int) .eq. 0 ) ) then
            call analyze_hydro_grid(mla,dt,dx,istep-n_steps_skip,rho=rho_old,temperature=Temp)           
         end if

         if ( (hydro_grid_int > 0) .and. &
              (n_steps_save_stats > 0) .and. &
              ( mod(istep-n_steps_skip,n_steps_save_stats) .eq. 0 ) ) then
              call save_hydro_grid(id=(istep-n_steps_skip)/n_steps_save_stats, step=istep)            
         end if

      end if

      ! write plotfile at specific intervals
      if ((plot_int.gt.0 .and. mod(istep,plot_int).eq.0) .or. (istep.eq.max_step)) then
         if (parallel_IOProcessor()) then
            write(*,*), 'writing plotfiles at timestep =', istep 
         end if
         call write_plotfileLM(mla,"plt",rho_old,rhotot_old,Temp,umac,pres,istep,dx,time)
      end if

      ! advance the solution by dt
      if (algorithm_type .eq. 0) then
         call bl_error("main_driver: inertial algorithm not written yet")
      else if (algorithm_type .eq. 1 .or. algorithm_type .eq. 2) then
         ! It appears to me there is no need to be passing eta, eta_ed etc. around here
         ! They are only used locally inside the routine when computing updates
         ! and do not appear to be needed here
         ! I think they should be local temps inside advance_timestep_overdamped
         ! Similar comment applies to diff_mass_fluxdiv,stoch_mass_fluxdiv
         ! I know they are also passed to initial_projection but it seems to me all this can be done locally with temps
         ! Also note that for advance_timestep_overdamped it is not necessary to call initial_projection
         ! In fact, if one makes this main_driver only do overdamped it can be simplified greatly
         ! One option is do this by converting this main_driver.f90 into overdamped_driver.f90 and simplifying it.
         call advance_timestep_overdamped(mla,umac,rho_old,rho_new,rhotot_old,rhotot_new, &
                                          pres,eta,eta_ed,kappa,Temp,Temp_ed, &
                                          diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                          dx,dt,time,the_bc_tower,n_rngs)
      end if

      ! increment simulation time
      istep = istep + 1
      time = time + dt

      ! set old state to new state
     do n=1,nlevs
        call multifab_copy_c(rho_old(n)   ,1,rho_new(n)   ,1,nspecies,rho_old(n)%ng)
        call multifab_copy_c(rhotot_old(n),1,rhotot_new(n),1,1       ,rhotot_old(n)%ng)
     end do

      ! print out the total mass to check conservation
     call sum_mass(rho_old, istep)

  end do

  ! print out the total mass to check conservation
  call sum_mass(rho_old, istep)

  !=======================================================
  ! Destroy multifabs and layouts
  !=======================================================

  if((abs(hydro_grid_int)>0) .or. (stats_int>0)) then
     call finalize_hydro_grid()
  end if

  call destroy_mass_stochastic(mla)
  call destroy_m_stochastic(mla)

  do n=1,nlevs
     call multifab_destroy(rho_old(n))
     call multifab_destroy(rhotot_old(n))
     call multifab_destroy(rho_new(n))
     call multifab_destroy(rhotot_new(n))
     call multifab_destroy(Temp(n))
     call multifab_destroy(diff_mass_fluxdiv(n))
     call multifab_destroy(stoch_mass_fluxdiv(n))
     call multifab_destroy(pres(n))
     call multifab_destroy(eta(n))
     call multifab_destroy(kappa(n))
     do i=1,dm
        call multifab_destroy(umac(n,i))
     end do
     do i=1,size(eta_ed,dim=2)
        call multifab_destroy(eta_ed(n,i))
        call multifab_destroy(Temp_ed(n,i))
     end do
  end do
  deallocate(lo,hi,dx)
  deallocate(rho_old,rhotot_old,Temp)
  deallocate(diff_mass_fluxdiv,stoch_mass_fluxdiv,umac)
  call destroy(mla)
  call bc_tower_destroy(the_bc_tower)

end subroutine main_driver
