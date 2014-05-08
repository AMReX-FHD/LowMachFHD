subroutine main_driver()

  use bl_IO_module
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
  use restart_module
  use checkpoint_module
  use estdt_module
  use convert_stag_module
  use probin_binarylm_module, only: probin_binarylm_init, max_step, print_int, &
                                   project_eos_int, initial_variance, analyze_binary, &
                                   conc_scal, barodiffusion_type, algorithm_type, restart
  use probin_common_module , only: probin_common_init, seed, dim_in, n_cells, &
                                   prob_lo, prob_hi, max_grid_size, &
                                   hydro_grid_int, n_steps_save_stats, n_steps_skip, &
                                   stats_int, variance_coef, chk_int, &
                                   bc_lo, bc_hi, fixed_dt, plot_int, advection_type
  use probin_gmres_module  , only: probin_gmres_init

  implicit none
  
  ! will be allocated with dm components
  integer, allocatable :: lo(:), hi(:)

  ! will be allocated with (nlevs,dm) components
  real(dp_t), allocatable :: dx(:,:)

  integer :: n,nlevs,i,dm,istep,n_cell

  real(kind=dp_t) :: dt,time,runtime1,runtime2,max_vel,av_mass
  real(kind=dp_t) :: l1, l2, linf

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
  type(multifab), allocatable :: gp_fc(:,:)        ! face-centered
  type(multifab), allocatable :: prim(:)           ! cell-centered
  type(multifab), allocatable :: pres(:)           ! cell-centered
  type(multifab), allocatable :: chi(:)            ! cell-centered
  type(multifab), allocatable :: chi_fc(:,:)       ! face-centered
  type(multifab), allocatable :: eta(:)            ! cell-centered
  type(multifab), allocatable :: eta_ed(:,:)       ! nodal (2d); edge-centered (3d)
  type(multifab), allocatable :: kappa(:)          ! cell-centered
  type(multifab), allocatable :: rhoc_fluxdiv(:)   ! cell-centered

  integer :: narg, farg, un, init_step, n_rngs
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

  ! dimensionality is set in inputs file
  dm = dim_in

  allocate(lo(dm),hi(dm))
  allocate(mold(nlevs,dm),mnew(nlevs,dm),umac(nlevs,dm))
  allocate(sold(nlevs),snew(nlevs),prim(nlevs),pres(nlevs))
  allocate(chi(nlevs),eta(nlevs),kappa(nlevs))
  allocate(rhoc_fluxdiv(nlevs))
  allocate(chi_fc(nlevs,dm),s_fc(nlevs,dm),gp_fc(nlevs,dm))
  if (dm .eq. 2) then
     allocate(eta_ed(nlevs,1))
  else if (dm .eq. 3) then
     allocate(eta_ed(nlevs,3))
  end if

  ! set grid spacing at each level
  ! the grid spacing is the same in each direction
  allocate(dx(nlevs,dm))
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
  if (nlevs .ge. 2) then
     call bl_error('ERROR: main_driver.f90, need to define dx for n>1')
  end if

  ! build pmask
  allocate(pmask(dm))
  pmask = .false.
  do i=1,dm
     if (bc_lo(i) .eq. PERIODIC .and. bc_hi(i) .eq. PERIODIC) then
        pmask(i) = .true.
     end if
  end do

  if (restart .ge. 0) then

     init_step = restart + 1

     ! build the ml_layout
     ! read in time and dt from checkpoint
     ! build and fill mold, sold, and pres
     ! build the mass flux divergence multifabs, and fill them for the inertial algorithm
     call initialize_from_restart(mla,time,dt,mold,sold,pres,rhoc_fluxdiv,pmask)

  else

     ! non-restart initialization
     ! need to build the ml_layout
     ! set time to zero
     ! build and fill mold, sold, and pres
     init_step = 1
     time = 0.d0

     ! tell mba how many levels and dmensionality of problem
     call ml_boxarray_build_n(mba,nlevs,dm)

     ! tell mba about the ref_ratio between levels
     ! mba%rr(n-1,i) is the refinement ratio between levels n-1 and n in direction i
     ! we use refinement ratio of 2 in every direction between all levels
     do n=2,nlevs
        mba%rr(n-1,:) = 2
     enddo

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
     
     ! build the ml_layout, mla
     call ml_layout_build(mla,mba,pmask)

     ! don't need this anymore - free up memory
     call destroy(mba)

     do n=1,nlevs
        do i=1,dm
           ! edge-momentum and velocity; 1 ghost cell
           call multifab_build_edge(mold(n,i),mla%la(n),1,1,i)
        end do

        ! conservative variables; 2 components (rho,rho1)
        ! need 2 ghost cells to average to ghost faces used in 
        ! converting m to umac in m ghost cells
        ! if using advection_type .ge. 1 (bds), need 3 ghost cells
        if (advection_type .ge. 1) then
           call multifab_build(sold(n),mla%la(n),2,3)
        else
           call multifab_build(sold(n),mla%la(n),2,2)
        end if

        ! pressure - need 1 ghost cell since we calculate its gradient
        call multifab_build(pres(n),mla%la(n),1,1)

        ! this stores divergence of stochastic and diffusive fluxes for rhoc
        call multifab_build(rhoc_fluxdiv(n),mla%la(n),1,0)
        call multifab_setval(rhoc_fluxdiv(n),0.d0,all=.true.)
     end do

     call init(mold,sold,pres,dx,mla,time)

  end if

  deallocate(pmask)

  ! tell the_bc_tower about max_levs, dm, and domain_phys_bc
  call initialize_bc(the_bc_tower,nlevs,dm,mla%pmask,2,1)
  do n=1,nlevs
     ! define level n of the_bc_tower
     call bc_tower_level_build(the_bc_tower,n,mla%la(n))
  end do

  do n=1,nlevs
     do i=1,dm
        ! edge-momentum and velocity; 1 ghost cell
        call multifab_build_edge(mnew(n,i),mla%la(n),1,1,i)
        call multifab_build_edge(umac(n,i),mla%la(n),1,1,i)
     end do
     ! conservative variables; 2 components (rho,rho1)
     ! need 2 ghost cells to average to ghost faces used in 
     ! converting m to umac in m ghost cells
     ! if using advection_type .ge. 1 (bds), need 3 ghost cells
     if (advection_type .ge. 1) then
        call multifab_build(snew(n),mla%la(n),2,3)
        call multifab_build(prim(n),mla%la(n),2,3)
     else
        call multifab_build(snew(n),mla%la(n),2,2)
        call multifab_build(prim(n),mla%la(n),2,2)
     end if

     ! s on faces, gp on faces
     do i=1,dm
        call multifab_build_edge( s_fc(n,i),mla%la(n),2,1,i)
        call multifab_build_edge(gp_fc(n,i),mla%la(n),1,0,i)
     end do

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

  end do

  ! convert cons to prim in valid region
  call convert_cons_to_prim(mla,sold,prim,.true.)

  if (restart .le. 0 .and. initial_variance .ne. 0.d0) then
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

  call average_cc_to_face(nlevs,sold,s_fc,1,scal_bc_comp,2,the_bc_tower%bc_tower_array)
  call convert_m_to_umac(mla,s_fc,mold,umac,.true.)

  do n=1,nlevs
     ! presure ghost cells
     call multifab_fill_boundary(pres(n))
     call multifab_physbc(pres(n),1,pres_bc_comp,1,the_bc_tower%bc_tower_array(n),dx(n,:))
  end do

  if (restart .le. 0) then

     if (initial_variance .ne. 0.d0) then
        ! umac is passed in as a temporary
        call add_m_fluctuations(mla,dx,initial_variance*variance_coef,sold,s_fc,mold,umac)
     end if
     
     if (barodiffusion_type .gt. 0) then
        ! this computes an initial guess at p using HSE
        call init_pres(mla,sold,pres,dx,the_bc_tower)
     end if

     ! set the initial time step
     if (fixed_dt .gt. 0.d0) then
        dt = fixed_dt
     else
        call estdt(mla,umac,dx,dt)
     end if

  end if

  ! compute grad p
  call compute_grad(mla,pres,gp_fc,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

  if (print_int .gt. 0) then
     if (parallel_IOProcessor()) write(*,*) "Initial conditions before initial projection:"
     call eos_check(mla,sold)
     call sum_momenta(mla,mold)
     do i=1,2
        av_mass = multifab_sum_c(sold(1),i,1,all=.false.)/product(n_cells(1:dm))
        if (parallel_IOProcessor()) then
           write(*,"(A,100G17.9)") "Conservation: rho=", i, av_mass
        end if
     end do
  end if

  ! initialize chi, eta, and kappa
  call compute_chi(mla,chi,chi_fc,prim,dx,the_bc_tower%bc_tower_array)
  call compute_eta(mla,eta,eta_ed,prim,dx,the_bc_tower%bc_tower_array)
  call compute_kappa(mla,kappa,prim,dx)

  ! allocate and build multifabs that will contain random numbers
  if (algorithm_type .eq. 0 .or. algorithm_type .eq. 1) then
     n_rngs = 1
  else if (algorithm_type .eq. 2) then
     n_rngs = 2
  end if
  call init_m_stochastic   (mla,n_rngs)
  call init_rhoc_stochastic(mla,n_rngs)

  ! fill the stochastic multifabs with a new set of random numbers
  call fill_m_stochastic(mla)  
  call fill_rhoc_stochastic(mla)  

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
              call initialize_hydro_grid_bin(mla,sold,dt,dx,un,2)
           else
              call initialize_hydro_grid(mla,sold,dt,dx, &
                      namelist_file=un, nspecies_in=2, nscal_in=0, exclude_last_species_in=.true., &
                      analyze_velocity=.true., analyze_density=.true., analyze_temperature=.false.)
           end if 
           close(unit=un)
        end if
     end if
  end if

  if (restart .gt. 0) then

     ! fill ghost cells for umac (but leave boundary values untouched)
     call fill_ghost_umac(mla,umac,eta_ed,dx,the_bc_tower)

     ! compute m that has all ghost values properly filled
     call convert_m_to_umac(mla,s_fc,mold,umac,.false.)

  else

     ! need to do an initial projection to get an initial velocity field
     call initial_projection(mla,mold,umac,sold,s_fc,prim,eta_ed,chi_fc,gp_fc, &
                             rhoc_fluxdiv,dx,dt, &
                             the_bc_tower)

     if (print_int .gt. 0) then
        if (parallel_IOProcessor()) write(*,*) "After initial projection:"
        call sum_momenta(mla,mold)
        do i=1,2
           av_mass = multifab_sum_c(sold(1),i,1,all=.false.)/product(n_cells(1:dm))
           if (parallel_IOProcessor()) then
              write(*,"(A,100G17.9)") "Conservation: rho=", i, av_mass
           end if
        end do
     end if

     ! write initial plotfile
     if (plot_int .gt. 0) then
        call write_plotfile(mla,mold,umac,sold,pres,dx,time,0)
     end if
     ! print out projection (average) and variance)
     if (stats_int .gt. 0) then
        if(analyze_binary) then   
           call print_stats_bin(mla,sold,mold,umac,prim,dx,0,time)
        else
           call print_stats(mla,dx,0,time,umac=umac,rho=sold)
        end if
     end if
  end if
  
  do istep=init_step,max_step

     runtime1 = parallel_wtime()

     if (fixed_dt .le. 0.d0) then
        call estdt(mla,umac,dx,dt)
     end if

     if ( (print_int .gt. 0 .and. mod(istep,print_int) .eq. 0) &
          .or. &
          (istep .eq. max_step) ) then
        if (parallel_IOProcessor()) then
           print*,"Begin Advance; istep =",istep,"DT =",dt,"TIME =",time
        end if
      
        if (fixed_dt .gt. 0.d0) then
           max_vel = 0.d0
           max_vel = max(max_vel,multifab_norm_inf_c(umac(1,1),1,1,all=.false.))
           max_vel = max(max_vel,multifab_norm_inf_c(umac(1,2),1,1,all=.false.))
           if (parallel_IOProcessor()) then
              print*,'effective advective CFL=',fixed_dt*max_vel/dx(1,1)
           end if
        end if
     end if

     ! advance the solution by dt
     if (algorithm_type .eq. 0) then
        call advance_timestep(mla,mold,mnew,umac,sold,snew,s_fc,prim,pres,chi,chi_fc, &
                              eta,eta_ed,kappa,rhoc_fluxdiv, &
                              gp_fc,dx,dt,time,the_bc_tower)
     else if (algorithm_type .eq. 1 .or. algorithm_type .eq. 2) then
        call advance_timestep_overdamped(mla,mnew,umac,sold,snew,s_fc,prim,pres, &
                                         chi,chi_fc,eta,eta_ed,kappa,dx,dt,time,the_bc_tower)
     end if

     ! increment simulation time
     time = time + dt

     if ( (print_int .gt. 0 .and. mod(istep,print_int) .eq. 0) &
          .or. &
          (istep .eq. max_step) ) then
       if (parallel_IOProcessor()) then
           print*,"End Advance; istep =",istep,"DT =",dt,"TIME =",time
        end if

        runtime2 = parallel_wtime() - runtime1
        call parallel_reduce(runtime1, runtime2, MPI_MAX, proc=parallel_IOProcessorNode())
        if (parallel_IOProcessor()) then
           print*,'Time to advance timestep: ',runtime1,' seconds'
        end if
     end if
      
     if ( (print_int .gt. 0 .and. mod(istep,print_int) .eq. 0) &
          .or. &
          (istep .eq. max_step) ) then
        if (parallel_IOProcessor()) write(*,*) "At time step ", istep, " t=", time  
        call eos_check(mla,snew)
        call sum_momenta(mla,mnew)
        do i=1,2
           av_mass = multifab_sum_c(sold(1),i,1,all=.false.)/product(n_cells(1:dm))
           if (parallel_IOProcessor()) then
              write(*,"(A,100G17.9)") "Conservation: rho=", i, av_mass
           end if
        end do
     end if

     ! project rho and rho1 back onto EOS
     if ( project_eos_int .gt. 0 .and. mod(istep,project_eos_int) .eq. 0) then
        call project_onto_eos(mla,snew)
     end if

      if (istep >= n_steps_skip) then

         if ( (plot_int > 0) .and. &
              ( mod(istep,plot_int) .eq. 0) ) then
            call write_plotfile(mla,mnew,umac,snew,pres,dx,time,istep)
         end if

         ! write checkpoint
         if ( (chk_int > 0) .and. &
              ( mod(istep,chk_int) .eq. 0) ) then
            call checkpoint_write(mla,snew,mnew,pres,rhoc_fluxdiv,time,dt,istep)
         end if

         ! print out projection (average) and variance
         if ( (stats_int > 0) .and. &
               (mod(istep,stats_int) .eq. 0) ) then
            if(analyze_binary) then   
               call print_stats_bin(mla,snew,mnew,umac,prim,dx,istep,time)
            else
               call print_stats(mla,dx,istep,time,umac=umac,rho=snew)            
            end if   
            if (analyze_binary.and.(hydro_grid_int<0)) then
               ! Do some specialized analysis for low Mach binary mixing studies
               call analyze_hydro_grid_bin(mla,snew,mnew,umac,prim,dt,dx, &
                                       istep,custom_analysis=.true.)
            end if   
         end if

         ! Add this snapshot to the average in HydroGrid
         if ( (hydro_grid_int > 0) .and. &
              ( mod(istep,hydro_grid_int) .eq. 0 ) ) then
            if(analyze_binary) then  
               call analyze_hydro_grid_bin(mla,snew,mnew,umac,prim,dt,dx, &
                                    istep,custom_analysis=.false.)
            else
               call analyze_hydro_grid(mla,dt,dx,istep,umac=umac,rho=snew)           
            end if                                    
         end if

         if ( (hydro_grid_int > 0) .and. &
              (n_steps_save_stats > 0) .and. &
              ( mod(istep,n_steps_save_stats) .eq. 0 ) ) then
            if(analyze_binary) then  
               call save_hydro_grid_bin(id=(istep)/n_steps_save_stats, step=istep)
            else
               call save_hydro_grid(id=(istep)/n_steps_save_stats, step=istep)            
            end if
         end if

      end if

     ! set old state to new state
     do n=1,nlevs
        call multifab_copy_c(sold(n),1,snew(n),1,2,sold(n)%ng)
        do i=1,dm
           call multifab_copy_c(mold(n,i),1,mnew(n,i),1,1,mold(n,i)%ng)
        end do
     end do
        
  end do
  
  !!!!!!!!!!!! convergence testing
  if(.false.) then
     call init(mold,sold,pres,dx,mla,time)
     do n=1,nlevs
        call multifab_sub_sub_c(sold(n),1,snew(n),1,2,0)
     end do
     linf = multifab_norm_inf_c(sold(1),2,1,all=.false.)
     l1 = multifab_norm_l1_c(sold(1),2,1,all=.false.)
     l2 = multifab_norm_l2_c(sold(1),2,1,all=.false.)
     n_cell = multifab_volume(sold(1)) / 2
     if (parallel_IOProcessor()) then
        print*,'linf error in rho*c',linf
        print*,'l1   error in rho*c',l1 / dble(n_cell)
        print*,'l2   error in rho*c',l2 / sqrt(dble(n_cell))
     end if
  end if   
  !!!!!!!!!!!! convergence testing

  ! destroy and deallocate multifabs that contain random numbers
  call destroy_m_stochastic(mla)
  call destroy_rhoc_stochastic(mla)

  if(abs(hydro_grid_int)>0 .or. stats_int>0) then
     ! Note that these will also write out statistics if n_steps_save_stats<=0
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
     call multifab_destroy(pres(n))
     call multifab_destroy(chi(n))
     call multifab_destroy(eta(n))
     call multifab_destroy(kappa(n))
     call multifab_destroy(rhoc_fluxdiv(n))
     do i=1,dm
        call multifab_destroy(mold(n,i))
        call multifab_destroy(mnew(n,i))
        call multifab_destroy(umac(n,i))
        call multifab_destroy(chi_fc(n,i))
        call multifab_destroy(s_fc(n,i))
        call multifab_destroy(gp_fc(n,i))
     end do
     do i=1,size(eta_ed,dim=2)
        call multifab_destroy(eta_ed(n,i))
     end do

  end do

  call destroy(mla)
  call bc_tower_destroy(the_bc_tower)

end subroutine main_driver
