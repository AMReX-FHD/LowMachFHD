subroutine main_driver()

  use bl_types
  use bl_IO_module
  use ml_boxarray_module
  use ml_layout_module
  use bl_error_module
  use ParallelRNGs
  use bc_module
  use define_bc_module
  use init_module
  use init_pres_module
  use initial_projection_module
  use write_plotfile_module
  use advance_timestep_module
  use advance_timestep_overdamped_module
  use convert_variables_module
  use div_and_grad_module
  use analysis_module
  use sum_momenta_module
  use stochastic_m_fluxdiv_module
  use stochastic_rhoc_fluxdiv_module
  use project_onto_eos_module
  use multifab_physbc_module
  use multifab_physbc_stag_module
  use fill_rho_ghost_cells_module
  use analyze_spectra_binary_module
  use analyze_spectra_module
  use estdt_module
  use convert_stag_module
  use probin_binarylm_module, only: probin_binarylm_init, max_step, print_int, &
                                   project_eos_int, initial_variance, analyze_binary, &
                                   conc_scal, barodiffusion_type, algorithm_type
  use probin_common_module , only: probin_common_init, seed, dim_in, n_cells, &
                                   prob_lo, prob_hi, max_grid_size, &
                                   hydro_grid_int, n_steps_save_stats, n_steps_skip, &
                                   stats_int, variance_coef, &
                                   bc_lo, bc_hi, fixed_dt, plot_int, advection_type
  use probin_gmres_module  , only: probin_gmres_init

  implicit none
  
  ! will be allocated with dm components
  integer, allocatable :: lo(:), hi(:)

  ! will be allocated with (nlevs,dm) components
  real(dp_t), allocatable :: dx(:,:)

  integer :: n,nlevs,i,dm,istep

  real(kind=dp_t) :: dt,time,runtime1,runtime2,max_vel,av_mass

  type(box)         :: bx
  type(ml_boxarray) :: mba
  type(ml_layout)   :: mla

  type(bc_tower) :: the_bc_tower

  logical, allocatable :: pmask(:)
  logical :: nodal_temp(3)

  type(multifab), allocatable :: mold(:,:)         ! face-based
  type(multifab), allocatable :: mnew(:,:)         ! face-based
  type(multifab), allocatable :: umac(:,:)         ! face-based
  type(multifab), allocatable :: sold(:)           ! cell-centered
  type(multifab), allocatable :: snew(:)           ! cell-centered
  type(multifab), allocatable :: s_fc(:,:)         ! face-centered
  type(multifab), allocatable :: gp_fc(:,:)       ! face-centered
  type(multifab), allocatable :: prim(:)           ! cell-centered
  type(multifab), allocatable :: pold(:)           ! cell-centered
  type(multifab), allocatable :: pnew(:)           ! cell-centered
  type(multifab), allocatable :: chi(:)            ! cell-centered
  type(multifab), allocatable :: chi_fc(:,:)       ! face-centered
  type(multifab), allocatable :: eta(:)            ! cell-centered
  type(multifab), allocatable :: eta_ed(:,:)       ! nodal (2d); edge-centered (3d)
  type(multifab), allocatable :: kappa(:)          ! cell-centered
  type(multifab), allocatable :: rhoc_d_fluxdiv(:) ! cell-centered
  type(multifab), allocatable :: rhoc_s_fluxdiv(:) ! cell-centered
  type(multifab), allocatable :: rhoc_b_fluxdiv(:) ! cell-centered

  ! special inhomogeneous boundary condition multifab
  ! vel_bc_n(nlevs,dm) are the normal velocities
  ! in 2D, vel_bc_t(nlevs,2) respresents
  !   1. y-velocity bc on x-faces (nodal)
  !   2. x-velocity bc on y-faces (nodal)
  ! in 3D, vel_bc_t(nlevs,6) represents
  !   1. y-velocity bc on x-faces (nodal in y and x)
  !   2. z-velocity bc on x-faces (nodal in z and x)
  !   3. x-velocity bc on y-faces (nodal in x and y)
  !   4. z-velocity bc on y-faces (nodal in z and y)
  !   5. x-velocity bc on z-faces (nodal in x and z)
  !   6. y-velocity bc on z-faces (nodal in y and z)
  type(multifab), allocatable :: vel_bc_n(:,:)
  type(multifab), allocatable :: vel_bc_t(:,:)

  integer :: narg, farg, un
  character(len=128) :: fname
  logical :: lexist

  ! uncomment this once lowMach_implicit/probin.f90 is written
  call probin_binarylm_init()
  call probin_common_init()
  call probin_gmres_init()

  ! Initialize random numbers *after* the global (root) seed has been set:
  call SeedParallelRNG(seed)

  ! in this example we fix nlevs to be 1
  ! for adaptive simulations where the grids change, cells at finer
  ! resolution don't necessarily exist depending on your tagging criteria,
  ! so max_levs isn't necessary equal to nlevs
  nlevs = 1

  dm = dim_in

  ! now that we have dm, we can allocate these
  allocate(lo(dm),hi(dm))

  ! now that we have nlevs and dm, we can allocate these
  allocate(dx(nlevs,dm))
  allocate(mold(nlevs,dm),mnew(nlevs,dm),umac(nlevs,dm),vel_bc_n(nlevs,dm))
  allocate(sold(nlevs),snew(nlevs),prim(nlevs),pold(nlevs),pnew(nlevs))
  allocate(chi(nlevs),eta(nlevs),kappa(nlevs))
  allocate(rhoc_d_fluxdiv(nlevs),rhoc_s_fluxdiv(nlevs),rhoc_b_fluxdiv(nlevs))
  allocate(chi_fc(nlevs,dm),s_fc(nlevs,dm),gp_fc(nlevs,dm))
  if (dm .eq. 2) then
     allocate(eta_ed(nlevs,1))
     allocate(vel_bc_t(nlevs,2))
  else if (dm .eq. 3) then
     allocate(eta_ed(nlevs,3))
     allocate(vel_bc_t(nlevs,6))
  end if

  ! tell mba how many levels and dmensionality of problem
  call ml_boxarray_build_n(mba,nlevs,dm)

  ! tell mba about the ref_ratio between levels
  ! mba%rr(n-1,i) is the refinement ratio between levels n-1 and n in direction i
  ! we use refinement ratio of 2 in every direction between all levels
  do n=2,nlevs
     mba%rr(n-1,:) = 2
  enddo

  ! set grid spacing at each level
  ! the grid spacing is the same in each direction
  dx(1,1:dm) = (prob_hi(1:dm)-prob_lo(1:dm)) / n_cells(1:dm)
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
  do n=2,nlevs
     dx(n,:) = dx(n-1,:) / mba%rr(n-1,:)
  end do

  ! create a box from (0,0) to (n_cells-1,n_cells-1)
  lo(1:dm) = 0
  hi(1:dm) = n_cells(1:dm)-1
  bx = make_box(lo,hi)

  ! tell mba about the problem domain at each level
  mba%pd(1) = bx
  do n=2,nlevs
     mba%pd(n) = refine(mba%pd(n-1),mba%rr((n-1),:))
  enddo

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

  ! tell the_bc_tower about max_levs, dm, and domain_phys_bc
  call initialize_bc(the_bc_tower,nlevs,dm,mla%pmask,2,1)
  do n=1,nlevs
     ! define level n of the_bc_tower
     call bc_tower_level_build(the_bc_tower,n,mla%la(n))
  end do

  do n=1,nlevs
     do i=1,dm
        ! edge-momentum and velocity; 1 ghost cell
        call multifab_build_edge(mold(n,i),mla%la(n)  ,1,1,i)
        call multifab_build_edge(mnew(n,i),mla%la(n)  ,1,1,i)
        call multifab_build_edge(umac(n,i),mla%la(n)  ,1,1,i)
     end do
     ! conservative variables; 2 components (rho,rho1)
     ! need 2 ghost cells to average to ghost faces used in 
     ! converting m to umac in m ghost cells
     ! if using advection_type .ge. 1 (bds), need 3 ghost cells
     if (advection_type .ge. 1) then
        call multifab_build(sold(n) ,mla%la(n),2,3)
        call multifab_build(snew(n) ,mla%la(n),2,3)
        call multifab_build(prim(n) ,mla%la(n),2,3)
     else
        call multifab_build(sold(n) ,mla%la(n),2,2)
        call multifab_build(snew(n) ,mla%la(n),2,2)
        call multifab_build(prim(n) ,mla%la(n),2,2)
     end if

     ! s on faces, gp on faces
     do i=1,dm
        call multifab_build_edge(  s_fc(n,i),mla%la(n),2,1,i)
        call multifab_build_edge(gp_fc(n,i),mla%la(n),1    ,0,i)
     end do

     ! pressure
     ! need 1 ghost cell since we calculate its gradient
     call multifab_build(pold(n),mla%la(n),1,1)
     call multifab_build(pnew(n),mla%la(n),1,1)

     ! transport coefficients
     call multifab_build(chi(n)  ,mla%la(n),1,1)
     call multifab_build(eta(n)  ,mla%la(n),1,1)
     call multifab_build(kappa(n),mla%la(n),1,1)

     ! chi on faces
     do i=1,dm
        call multifab_build_edge(chi_fc(n,i),mla%la(n),1,0,i)
     end do

     ! eta on nodes (2d) or edges (3d)
     if (dm .eq. 2) then
        call multifab_build_nodal(eta_ed(n,1),mla%la(n),1,0)
     else
        nodal_temp(1) = .true.
        nodal_temp(2) = .true.
        nodal_temp(3) = .false.
        call multifab_build(eta_ed(n,1),mla%la(n),1,0,nodal_temp)
        nodal_temp(1) = .true.
        nodal_temp(2) = .false.
        nodal_temp(3) = .true.
        call multifab_build(eta_ed(n,2),mla%la(n),1,0,nodal_temp)
        nodal_temp(1) = .false.
        nodal_temp(2) = .true.
        nodal_temp(3) = .true.
        call multifab_build(eta_ed(n,3),mla%la(n),1,0,nodal_temp)
     end if

     ! this stores divergence of stochastic and diffusive fluxes for rhoc
     call multifab_build(rhoc_d_fluxdiv(n),mla%la(n),1,0)
     call multifab_build(rhoc_s_fluxdiv(n),mla%la(n),1,0)
     call multifab_build(rhoc_b_fluxdiv(n),mla%la(n),1,0)
     call multifab_setval(rhoc_d_fluxdiv(n),0.d0,all=.true.)
     call multifab_setval(rhoc_s_fluxdiv(n),0.d0,all=.true.)
     call multifab_setval(rhoc_b_fluxdiv(n),0.d0,all=.true.)

     ! boundary conditions
     do i=1,dm
        call multifab_build_edge(vel_bc_n(n,i),mla%la(n),1,0,i)
     end do
     if (dm .eq. 2) then
        ! y-velocity bc on x-faces (nodal)
        call multifab_build_nodal(vel_bc_t(n,1),mla%la(n),1,0)
        ! x-velocity bc on y-faces (nodal)
        call multifab_build_nodal(vel_bc_t(n,2),mla%la(n),1,0)
     else
        ! y-velocity bc on x-faces (nodal in y and x)
        nodal_temp(1) = .true.
        nodal_temp(2) = .true.
        nodal_temp(3) = .false.
        call multifab_build(vel_bc_t(n,1),mla%la(n),1,0,nodal_temp)
        ! z-velocity bc on x-faces (nodal in z and x)
        nodal_temp(1) = .true.
        nodal_temp(2) = .false.
        nodal_temp(3) = .true.
        call multifab_build(vel_bc_t(n,2),mla%la(n),1,0,nodal_temp)
        ! x-velocity bc on y-faces (nodal in x and y)
        nodal_temp(1) = .true.
        nodal_temp(2) = .true.
        nodal_temp(3) = .false.
        call multifab_build(vel_bc_t(n,3),mla%la(n),1,0,nodal_temp)
        ! z-velocity bc on y-faces (nodal in z and y)
        nodal_temp(1) = .false.
        nodal_temp(2) = .true.
        nodal_temp(3) = .true.
        call multifab_build(vel_bc_t(n,4),mla%la(n),1,0,nodal_temp)
        ! x-velocity bc on z-faces (nodal in x and z)
        nodal_temp(1) = .true.
        nodal_temp(2) = .false.
        nodal_temp(3) = .true.
        call multifab_build(vel_bc_t(n,5),mla%la(n),1,0,nodal_temp)
        ! y-velocity bc on z-faces (nodal in y and z)
        nodal_temp(1) = .false.
        nodal_temp(2) = .true.
        nodal_temp(3) = .true.
        call multifab_build(vel_bc_t(n,6),mla%la(n),1,0,nodal_temp)
     end if

  end do

  time = 0.d0

  ! initialize sold = s^0 and mold = m^0
  call init(mold,sold,pold,dx,mla,time)

  if (initial_variance .ne. 0.d0) then
     call average_cc_to_face(nlevs,sold,s_fc,1,scal_bc_comp,1,the_bc_tower%bc_tower_array)
     call add_m_fluctuations(mla,dx,initial_variance*variance_coef,sold,s_fc,mold,umac)
  end if

  if (barodiffusion_type .gt. 0) then
     ! this computes an initial guess at p using HSE
     call init_pres(mla,sold,pold,dx,the_bc_tower)
  end if

  ! compute grad p
  call compute_grad(mla,pold,gp_fc,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

  if (print_int .gt. 0) then
     call eos_check(mla,sold)
     call sum_momenta(mla,mold)
     do i=1,2
        av_mass = multifab_sum_c(sold(1),i,1,all=.false.)/product(n_cells(1:dm))
        if (parallel_IOProcessor()) then
           write(*,"(A,100G17.9)") "CONSERVE: <rho_i>=", av_mass
        end if
     end do
  end if

  ! convert cons to prim in valid region
  call convert_cons_to_prim(mla,sold,prim,.true.)

  if (initial_variance .ne. 0.d0) then
     call add_c_fluctuations(mla,dx,initial_variance*variance_coef*conc_scal,prim,sold)
  end if

  ! fill ghost cells for prim
  do n=1,nlevs
     call multifab_fill_boundary(prim(n))
     call multifab_physbc(prim(n),2,scal_bc_comp+1,1,the_bc_tower%bc_tower_array(n),dx(n,:))
     call fill_rho_ghost_cells(prim(n),the_bc_tower%bc_tower_array(n))
  end do

  ! convert prim to cons in valid and ghost region
  ! now cons has properly filled ghost cells
  call convert_cons_to_prim(mla,sold,prim,.false.)

  ! initialize chi, eta, and kappa
  call compute_chi(mla,chi,chi_fc,prim,dx,the_bc_tower%bc_tower_array)
  call compute_eta(mla,eta,eta_ed,prim,dx,the_bc_tower%bc_tower_array)
  call compute_kappa(mla,kappa,prim,dx)

  ! set inhomogeneous bc condition for velocities
  call set_inhomogeneous_vel_bcs(mla,vel_bc_n,vel_bc_t,eta_ed,dx, &
                                 the_bc_tower%bc_tower_array)

  ! allocate and build multifabs that will contain random numbers
  if (algorithm_type .eq. 0 .or. algorithm_type .eq. 1) then
     call init_m_stochastic(mla,1)
     call init_rhoc_stochastic(mla,1)
  else if (algorithm_type .eq. 2) then
     call init_m_stochastic(mla,2)
     call init_rhoc_stochastic(mla,2)
  end if

  ! fill the stochastic multifabs with a new set of random numbers
  call fill_m_stochastic(mla)  
  call fill_rhoc_stochastic(mla)  

  ! set the initial time step
  if (fixed_dt .gt. 0.d0) then
     dt = fixed_dt
  else
     call average_cc_to_face(nlevs,sold,s_fc,1,scal_bc_comp,1,the_bc_tower%bc_tower_array)
     call convert_m_to_umac(mla,s_fc,mold,umac,.true.)
     call estdt(mla,umac,dx,dt)
  end if

  if(abs(hydro_grid_int)>0 .or. stats_int>0) then
     narg = command_argument_count()
     farg = 1
     if (narg >= 1) then
        call get_command_argument(farg, value = fname)
        inquire(file = fname, exist = lexist )
        if ( lexist ) then
           un = unit_new()
           open(unit=un, file = fname, status = 'old', action = 'read')
           if(analyze_binary) then
              call initialize_hydro_grid_bin(mla,sold,mold,dt,dx,un,2)
           else
              call initialize_hydro_grid(mla,sold,dt,dx, &
                      namelist_file=un, nspecies_in=2, nscal_in=0, exclude_last_species_in=.true., &
                      analyze_velocity=.true., analyze_density=.true., analyze_temperature=.false.)
           end if 
           close(unit=un)
        end if
     end if
  end if

  ! need to do an initial projection to get an initial velocity field
  call initial_projection(mla,mold,umac,sold,s_fc,prim,chi_fc,gp_fc,rhoc_d_fluxdiv, &
                          rhoc_s_fluxdiv,rhoc_b_fluxdiv,dx,dt, &
                          the_bc_tower,vel_bc_n,vel_bc_t)

  if (print_int .gt. 0) then
     call sum_momenta(mla,mold)
     do i=1,2
        av_mass = multifab_sum_c(sold(1),i,1,all=.false.)/product(n_cells(1:dm))
        if (parallel_IOProcessor()) then
           write(*,"(A,100G17.9)") "CONSERVE: <rho_i>=", av_mass
        end if
     end do
  end if

  ! write initial plotfile
  if (plot_int .gt. 0) then
     call write_plotfile(mla,mold,umac,sold,pold,dx,time,0)
  end if
  ! print out projection (average) and variance)
  if (stats_int .gt. 0) then
     if(analyze_binary) then   
        call print_stats_bin(mla,sold,mold,umac,prim,dx,0,time)
     else
        call print_stats(mla,dx,0,time,umac=umac,rho=sold)
     end if
  end if
  
  do istep=1,max_step

     runtime1 = parallel_wtime()

     if (fixed_dt .le. 0.d0) then
        call average_cc_to_face(nlevs,sold,s_fc,1,scal_bc_comp,1,the_bc_tower%bc_tower_array)
        call convert_m_to_umac(mla,s_fc,mold,umac,.true.)
        call estdt(mla,umac,dx,dt)
     end if

     if (parallel_IOProcessor()) then
        print*,"Begin Advance; istep =",istep,"DT =",dt,"TIME =",time
     end if

     if (fixed_dt .gt. 0.d0) then
        call convert_m_to_umac(mla,s_fc,mold,umac,.true.)
        max_vel = 0.d0
        max_vel = max(max_vel,multifab_norm_inf_c(umac(1,1),1,1,all=.false.))
        max_vel = max(max_vel,multifab_norm_inf_c(umac(1,2),1,1,all=.false.))
        if (parallel_IOProcessor()) then
           print*,'effective advective CFL=',fixed_dt*max_vel/dx(1,1)
        end if
     end if

     ! advance the solution by dt
     if (algorithm_type .eq. 0) then
        call advance_timestep(mla,mold,mnew,umac,sold,snew,s_fc,prim,pold,pnew,chi,chi_fc, &
                              eta,eta_ed,kappa,rhoc_d_fluxdiv,rhoc_s_fluxdiv,rhoc_b_fluxdiv, &
                              gp_fc,dx,dt,time,the_bc_tower,vel_bc_n,vel_bc_t)
     else if (algorithm_type .eq. 1 .or. algorithm_type .eq. 2) then
        call advance_timestep_overdamped(mla,mnew,umac,sold,snew,s_fc,prim,pold,pnew, &
                                         chi,chi_fc,eta,eta_ed,kappa,dx,dt,time,the_bc_tower, &
                                         vel_bc_n,vel_bc_t)
     end if

     ! increment simulation time
     time = time + dt

    if (parallel_IOProcessor()) then
        print*,"End Advance; istep =",istep,"DT =",dt,"TIME =",time
     end if

     runtime2 = parallel_wtime() - runtime1
     call parallel_reduce(runtime1, runtime2, MPI_MAX, proc=parallel_IOProcessorNode())
     if (parallel_IOProcessor()) then
        print*,'Time to advance timestep: ',runtime1,' seconds'
     end if

     ! project rho and rho1 back onto EOS
     if ( project_eos_int .gt. 0 .and. mod(istep,project_eos_int) .eq. 0) then
        call project_onto_eos(mla,snew)
     end if

     if ( (print_int .gt. 0 .and. mod(istep,print_int) .eq. 0) &
          .or. &
          (istep .eq. max_step) ) then
        call eos_check(mla,snew)
        call sum_momenta(mla,mnew)
        do i=1,2
           av_mass = multifab_sum_c(sold(1),i,1,all=.false.)/product(n_cells(1:dm))
           if (parallel_IOProcessor()) then
              write(*,"(A,100G17.9)") "CONSERVE: <rho_i>=", av_mass
           end if
        end do
     end if

      if (istep > n_steps_skip) then

         ! write plotfile
         if ( (plot_int > 0) .and. &
              ( mod(istep-n_steps_skip,plot_int) .eq. 0) ) then
            call write_plotfile(mla,mnew,umac,snew,pnew,dx,time,istep-n_steps_skip)
         end if

         ! print out projection (average) and variance
         if ( (stats_int > 0) .and. &
               (mod(istep-n_steps_skip,stats_int) .eq. 0) ) then
            if(analyze_binary) then   
               call print_stats_bin(mla,snew,mnew,umac,prim,dx,istep-n_steps_skip,time)
            else
               call print_stats(mla,dx,istep-n_steps_skip,time,umac=umac,rho=snew)            
            end if   
            if (analyze_binary.and.(hydro_grid_int<0)) then
               ! Do some specialized analysis for low Mach binary mixing studies
               call analyze_hydro_grid_bin(mla,snew,mnew,umac,prim,dt,dx, &
                                       istep-n_steps_skip,custom_analysis=.true.)
            end if   
         end if

         ! Add this snapshot to the average in HydroGrid
         if ( (hydro_grid_int > 0) .and. &
              ( mod(istep-n_steps_skip,hydro_grid_int) .eq. 0 ) ) then
            if(analyze_binary) then  
               call analyze_hydro_grid_bin(mla,snew,mnew,umac,prim,dt,dx, &
                                    istep-n_steps_skip,custom_analysis=.false.)
            else
               call analyze_hydro_grid(mla,dt,dx,istep-n_steps_skip,umac=umac,rho=snew)           
            end if                                    
         end if

         if ( (hydro_grid_int > 0) .and. &
              (n_steps_save_stats > 0) .and. &
              ( mod(istep-n_steps_skip,n_steps_save_stats) .eq. 0 ) ) then
            if(analyze_binary) then  
               call save_hydro_grid_bin(id=(istep-n_steps_skip)/n_steps_save_stats, step=istep)
            else
               call save_hydro_grid(id=(istep-n_steps_skip)/n_steps_save_stats, step=istep)            
            end if
         end if

      end if


     ! set old state to new state
     do n=1,nlevs
        call multifab_copy_c(pold(n),1,pnew(n),1,    1,pold(n)%ng)
        call multifab_copy_c(sold(n),1,snew(n),1,2,sold(n)%ng)
        do i=1,dm
           call multifab_copy_c(mold(n,i),1,mnew(n,i),1,1,mold(n,i)%ng)
        end do
     end do
        
  end do

  !!!!!!!!!!!! convergence testing
!  call init(mold,sold,pold,dx,mla,time)
!  do n=1,nlevs
!     call multifab_sub_sub_c(sold(n),1,snew(n),1,2,0)
!  end do
!  linf = multifab_norm_inf_c(sold(1),2,1,all=.false.)
!  l1 = multifab_norm_l1_c(sold(1),2,1,all=.false.)
!  l2 = multifab_norm_l2_c(sold(1),2,1,all=.false.)
!  n_cell = multifab_volume(sold(1)) / 2
!  if (parallel_IOProcessor()) then
!     print*,'linf error in rho*c',linf
!     print*,'l1   error in rho*c',l1 / dble(n_cell)
!     print*,'l2   error in rho*c',l2 / sqrt(dble(n_cell))
!  end if
  !!!!!!!!!!!! convergence testing

  ! destroy and deallocate multifabs that contain random numbers
  call destroy_m_stochastic(mla)
  call destroy_rhoc_stochastic(mla)

  if(abs(hydro_grid_int)>0 .or. stats_int>0) then
     if(analyze_binary) then
        call finalize_hydro_grid_bin()
     else
        call finalize_hydro_grid()
     end if
  end if

  do n=1,nlevs
     call multifab_destroy(sold(n))
     call multifab_destroy(snew(n))
     call multifab_destroy(prim(n))
     call multifab_destroy(pold(n))
     call multifab_destroy(pnew(n))
     call multifab_destroy(chi(n))
     call multifab_destroy(eta(n))
     call multifab_destroy(kappa(n))
     call multifab_destroy(rhoc_d_fluxdiv(n))
     call multifab_destroy(rhoc_s_fluxdiv(n))
     call multifab_destroy(rhoc_b_fluxdiv(n))
     do i=1,dm
        call multifab_destroy(mold(n,i))
        call multifab_destroy(mnew(n,i))
        call multifab_destroy(umac(n,i))
        call multifab_destroy(vel_bc_n(n,i))
        call multifab_destroy(chi_fc(n,i))
        call multifab_destroy(s_fc(n,i))
        call multifab_destroy(gp_fc(n,i))
     end do
     do i=1,size(eta_ed,dim=2)
        call multifab_destroy(eta_ed(n,i))
     end do
     do i=1,size(vel_bc_t,dim=2)
        call multifab_destroy(vel_bc_t(n,i))
     end do

  end do

  call destroy(mla)
  call bc_tower_destroy(the_bc_tower)

end subroutine main_driver
