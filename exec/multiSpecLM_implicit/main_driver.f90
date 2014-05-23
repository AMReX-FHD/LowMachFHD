subroutine main_driver()

  use boxlib
  use multifab_module
  use bl_IO_module
  use layout_module
  use init_module
  use initial_projection_module
  use write_plotfileLM_module
  use advance_diffusion_module
  use advance_timestep_overdamped_module
  use define_bc_module
  use bc_module
  use analysis_module
  use analyze_spectra_module
  use eos_check_module
  use stochastic_mass_fluxdiv_module
  use stochastic_m_fluxdiv_module
  use ParallelRNGs 
  use convert_mass_variables_module
  use convert_stag_module
  use probin_common_module, only: prob_lo, prob_hi, n_cells, dim_in, hydro_grid_int, &
       k_B, max_grid_size, n_steps_save_stats, n_steps_skip, plot_int, seed, stats_int, &
       bc_lo, bc_hi, probin_common_init, advection_type, fixed_dt
  use probin_multispecies_module, only: nspecies, rho_init, rho_bc, chi, &
       max_step, mol_frac_bc_comp, print_error_norms, rho_part_bc_comp, &
       start_time, molmass, temp_bc_comp, timeinteg_type, use_stoch, variance_coef_mass, &
       probin_multispecies_init
  use probin_gmres_module, only: probin_gmres_init
 
  implicit none

  ! quantities will be allocated with dm components
  integer, allocatable :: lo(:), hi(:)

  ! quantities will be allocated with (nlevs,dm) components
  real(kind=dp_t), allocatable :: dx(:,:)
  real(kind=dp_t)              :: dt,time,loc_param
  integer                      :: n,nlevs,i,j,dm,istep,step_count,ng_s
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
  real(kind=dp_t),allocatable  :: covW(:,:) 
  real(kind=dp_t),allocatable  :: covW_theo(:,:) 
  real(kind=dp_t),allocatable  :: wiwjt(:,:) 
  real(kind=dp_t),allocatable  :: wit(:) 

  ! For HydroGrid
  integer :: narg, farg, un, n_cell, rng, n_rngs
  character(len=128) :: fname
  logical :: lexist
  logical :: nodal_temp(3)
  
  !==============================================================
  ! Initialization
  !==============================================================

  call probin_common_init()
  call probin_multispecies_init() 
  call probin_gmres_init()
  
  if(.true.) then ! Confirm that gcc read the input file correctly
    write(*,*) "rho_init=", rho_init(1:2,1:nspecies)
    write(*,*) "rho_bc=", rho_bc(1:dim_in,1:2,1:nspecies)
  end if

  ! for time being, we fix nlevs to be 1. for adaptive simulations where the grids 
  ! change, cells at finer resolution don't necessarily exist depending on your 
  ! tagging criteria, so max_levs isn't necessary equal to nlevs
  nlevs = 1
  dm = dim_in

  n_rngs = 1
 
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
  allocate(covW(nspecies,nspecies))
  allocate(covW_theo(nspecies,nspecies))
  allocate(wiwjt(nspecies,nspecies))
  allocate(wit(nspecies))

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
  ! scal_bc_comp+nspecies+1 = mol_frac
  ! scal_bc_comp+2*nspecies+1 = temp_bc_comp = temperature
  ! scal_bc_comp+2*nspecies+2 = tran_bc_comp = diff_coef
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

  ! Initialize random numbers *after* the global (root) seed has been set:
  if(use_stoch) call SeedParallelRNG(seed)

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
  do n=1,nlevs
     call multifab_setval(eta(n), 1.d-2, all=.true.)
  end do
  if (dm .eq. 2) then
     call average_cc_to_node(nlevs,eta,eta_ed(:,1),1,tran_bc_comp,1,the_bc_tower%bc_tower_array)
  else if (dm .eq. 3) then
     call average_cc_to_edge(nlevs,eta,eta_ed,1,tran_bc_comp,1,the_bc_tower%bc_tower_array)
  end if

  ! initialize pressure and velocity
  do n=1,nlevs
     call multifab_setval(pres(n), 0.d0, all=.true.)
     do i=1,dm
        call multifab_setval(umac(n,i), 0.d0, all=.true.)
     end do
  end do

  ! initialize multifabs that hold random fluxes
  call init_mass_stochastic(mla,n_rngs)
  call init_m_stochastic(mla,n_rngs)

  ! fill random flux multifabs
  call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)
  call fill_m_stochastic(mla)

  ! choice of time step with a diffusive CFL of 0.1; CFL=minimum[dx^2/(2*chi)]; 
  ! chi is the largest eigenvalue of diffusion matrix to be input for n-species
  dt = fixed_dt

  if (dt .gt. dx(1,1)**2/(chi*2.d0*dm)) then
     call bl_error("time step violates diffusive mass cfl")
  end if
  
  if (parallel_IOProcessor()) then
     if(timeinteg_type .eq. 1) write(*,*) "Using Euler method"
     if(timeinteg_type .eq. 2) write(*,*) "Using Predictor-corrector method"
     if(timeinteg_type .eq. 3) write(*,*) "Using Midpoint method"
     if(timeinteg_type .eq. 4) write(*,*) "Using Runge-Kutta 3 method"
     write(*,*) "Using time step dt =", dt
     if(use_stoch) write(*,*), "Using noise variance =", sqrt(2.d0*k_B*&
                               variance_coef_mass/(product(dx(1,1:dm))*dt))
  end if

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
           call initialize_hydro_grid(mla,rho_old,dt,dx, namelist_file=un, &
                                      nspecies_in=nspecies, &
                                      nscal_in=0, &
                                      exclude_last_species_in=.false., &
                                      analyze_velocity=.false., &
                                      analyze_density=.true., &
                                      analyze_temperature=.true.) 
           
           close(unit=un)
        end if
     end if
  end if

  ! initial projection - only needed for inertial algorithm
!  call initial_projection(mla,umac,rho_old,rhotot_old,diff_mass_fluxdiv,stoch_mass_fluxdiv, &
!                          Temp,dt,dx,n_rngs,the_bc_tower)

  !=======================================================
  ! Begin time stepping loop
  !=======================================================

  ! free up memory counters 
  step_count = 0.d0
  covW       = 0.d0 
  covW_theo  = 0.d0 
  wit        = 0.d0 
  wiwjt      = 0.d0 
  istep      = 0

  do while(istep<=max_step)

      if (parallel_IOProcessor()) then
         print*,"Begin Advance; istep =",istep,"dt =",dt,"time =",time
      end if

      ! We do the analysis first so we include the initial condition in the files if n_steps_skip=0
      if (istep >= n_steps_skip) then
         ! Compute covariances manually for initial testing (HydroGrid now does the same)
         call compute_cov(mla,rho_old,wit,wiwjt)    
         step_count = step_count + 1 

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
         write(*,*), 'writing plotfiles at timestep =', istep 
         call write_plotfileLM(mla,"plt",rho_old,rhotot_old,Temp,umac,pres,istep,dx,time)
      end if

      ! advance the solution by dt
      call advance_timestep_overdamped(mla,umac,rho_old,rho_new,rhotot_old,rhotot_new, &
                                       pres,eta,eta_ed,kappa,Temp,Temp_ed, &
                                       diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                       dx,dt,time,the_bc_tower,n_rngs)

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

  ! print out the standard deviation
  if (parallel_IOProcessor()) then
     if(use_stoch .and. variance_coef_mass .ne. 0.d0 .and. step_count .ne. 0.d0) then
     
     write(*,*), ''
     write(1,*), 'Normalized numeric cov of W'
     do i=1,nspecies
        do j=1,nspecies
           covW(i,j) = (wiwjt(i,j)/real(step_count) - wit(i)*wit(j)/real(step_count)**2)/variance_coef_mass
        end do
        write(1,*), covW(i,:)
     end do

     if(nspecies .eq. 2) then 
     
        covW_theo(1,1) = (molmass(1)*rho_init(1,2) + molmass(2)*rho_init(1,1))*rho_init(1,1)*rho_init(1,2)/(&
                          product(dx(1,1:dm))*(rho_init(1,1)+rho_init(1,2))**4) 
        covW_theo(1,2) = -covW_theo(1,1) 
        covW_theo(2,1) = covW_theo(1,2) 
        covW_theo(2,2) = covW_theo(1,1) 

     else if(nspecies .eq. 3) then
        
        covW_theo(1,1) = (rho_init(1,1)*(molmass(2)*rho_init(1,1)*rho_init(1,2) + molmass(3)*rho_init(1,1)*&
                         rho_init(1,3) + molmass(1)*(rho_init(1,2) + rho_init(1,3))**2))/(product(dx(1,1:dm))*&
                         (rho_init(1,1) + rho_init(1,2) + rho_init(1,3))**4) 
        covW_theo(1,2) = -((rho_init(1,1)*rho_init(1,2)*(-(molmass(3)*rho_init(1,3)) + molmass(2)*(rho_init(1,1) +& 
                         rho_init(1,3)) + molmass(1)*(rho_init(1,2) + rho_init(1,3))))/(product(dx(1,1:dm))*(&
                         rho_init(1,1) + rho_init(1,2) + rho_init(1,3))**4)) 
        covW_theo(1,3) = -((rho_init(1,1)*rho_init(1,3)*(-(molmass(2)*rho_init(1,2)) + molmass(3)*(rho_init(1,1) +& 
                         rho_init(1,2)) + molmass(1)*(rho_init(1,2) + rho_init(1,3))))/(product(dx(1,1:dm))*(&
                         rho_init(1,1) + rho_init(1,2) + rho_init(1,3))**4))
        covW_theo(2,1) = covW_theo(1,2) 
        covW_theo(2,2) = (rho_init(1,2)*(molmass(2)*(rho_init(1,1) + rho_init(1,3))**2 + rho_init(1,2)*(molmass(1)*&
                         rho_init(1,1) + molmass(3)*rho_init(1,3))))/(product(dx(1,1:dm))*(rho_init(1,1) +& 
                         rho_init(1,2) + rho_init(1,3))**4)
        covW_theo(2,3) = -((rho_init(1,2)*rho_init(1,3)*(-(molmass(1)*rho_init(1,1)) + molmass(3)*(rho_init(1,1) +& 
                         rho_init(1,2)) + molmass(2)*(rho_init(1,1) + rho_init(1,3))))/(product(dx(1,1:dm))*(&
                         rho_init(1,1) + rho_init(1,2) + rho_init(1,3))**4)) 
        covW_theo(3,1) = covW_theo(1,3) 
        covW_theo(3,2) = covW_theo(2,3) 
        covW_theo(3,3) = (rho_init(1,3)*(molmass(3)*(rho_init(1,1) + rho_init(1,2))**2 + (molmass(1)*rho_init(1,1) +& 
                         molmass(2)*rho_init(1,2))*rho_init(1,3)))/(product(dx(1,1:dm))*(rho_init(1,1) +& 
                         rho_init(1,2) + rho_init(1,3))**4)
 
      else if(nspecies .eq. 4) then
        
        loc_param = 1.0d0/(product(dx(1,1:dm))*(rho_init(1,1) + rho_init(1,2) + rho_init(1,3) +& 
                    rho_init(1,4)))
        
        covW_theo(1,1) = 0.174d0*loc_param 
        covW_theo(1,2) = -0.024d0*loc_param
        covW_theo(1,3) = -0.018d0*loc_param
        covW_theo(1,4) = -0.132d0*loc_param 
 
        covW_theo(2,1) = covW_theo(1,2) 
        covW_theo(2,2) = 0.2115d0*loc_param
        covW_theo(2,3) = -0.0195d0*loc_param
        covW_theo(2,4) = -0.168d0*loc_param

        covW_theo(3,1) = covW_theo(1,3)
        covW_theo(3,2) = covW_theo(2,3)
        covW_theo(3,3) = 0.1135d0*loc_param
        covW_theo(3,4) = -0.076d0*loc_param 

        covW_theo(4,1) = covW_theo(1,4)
        covW_theo(4,2) = covW_theo(2,4)
        covW_theo(4,3) = covW_theo(3,4)
        covW_theo(4,4) = 0.376d0*loc_param

      else  
        write(*,*), 'analytic covariance of W for nspecies > 4 is not coded'
     end if

     ! correction made for infinite to periodic approximation of covariance
     n_cell = multifab_volume(rhotot_old(1))
     covW_theo = (1 - 1/n_cell)*covW_theo
     
     write(2,*), 'analytic cov of W' 
     do i=1,nspecies
        write(2,*) covW_theo(i,:)
     end do

     write(*,*), 'linf-norm of cov of W for',nspecies,' species is =', maxval(abs(covW_theo - covW))
     write(*,*), 'l1-norm of cov of W for',  nspecies,' species is =', sum(abs(covW_theo - covW))
     write(*,*), 'l2-norm of cov of W for',  nspecies,' species is =', sqrt(sum((covW_theo - covW)**2))
 
     end if
  end if

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
  deallocate(rho_old,rhotot_old,Temp,diff_mass_fluxdiv,stoch_mass_fluxdiv,umac)
  deallocate(covW,covW_theo,wiwjt,wit)
  call destroy(mla)
  call bc_tower_destroy(the_bc_tower)

end subroutine main_driver
