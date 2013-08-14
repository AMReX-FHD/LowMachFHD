subroutine main_driver()

  use bl_types
  use ml_boxarray_module
  use ml_layout_module
  use bl_error_module
  use ParallelRNGs
  use bc_module
  use define_bc_module
  use init_module
  use initial_projection_module
  use write_plotfile_module
  use advance_timestep_module
  use convert_variables_module
  use analysis_module
  use mk_stochastic_fluxdiv_module
  use project_onto_eos_module
  use multifab_physbc_module
  use multifab_physbc_stag_module
  use probin_lowmach_module, only: probin_lowmach_init, max_step, nscal, print_int, &
                                   project_eos_int, visc_coef
  use probin_common_module , only: probin_common_init, seed, dim_in, n_cells, &
                                   prob_lo, prob_hi, max_grid_size, &
                                   bc_lo, bc_hi, fixed_dt, plot_int, visc_type
  use probin_gmres_module  , only: probin_gmres_init

  implicit none

  ! will be allocated with dm components
  integer, allocatable :: lo(:), hi(:)

  ! will be allocated with (nlevs,dm) components
  real(dp_t), allocatable :: dx(:,:)

  integer :: n,nlevs,i,dm,istep

  real(kind=dp_t) :: time

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
  type(multifab), allocatable :: gp0_fc(:,:)       ! face-centered
  type(multifab), allocatable :: prim(:)           ! cell-centered
  type(multifab), allocatable :: pres(:)           ! cell-centered
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

  ! uncomment this once lowMach_implicit/probin.f90 is written
  call probin_lowmach_init()
  call probin_common_init()
  call probin_gmres_init()

  ! inputs file error checking
  if (visc_coef .lt. 0.d0 .and. visc_type .gt. 0) then
     call bl_error("negative visc_coef requires negative visc_type")
  end if

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
  allocate(sold(nlevs),snew(nlevs),prim(nlevs),pres(nlevs))
  allocate(chi(nlevs),eta(nlevs),kappa(nlevs))
  allocate(rhoc_d_fluxdiv(nlevs),rhoc_s_fluxdiv(nlevs),rhoc_b_fluxdiv(nlevs))
  allocate(chi_fc(nlevs,dm),s_fc(nlevs,dm),gp0_fc(nlevs,dm))
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
  dx(1,1:dm) = (prob_hi(1)-prob_lo(1)) / n_cells(1:dm)
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
     call multifab_build(sold(n) ,mla%la(n),nscal,2)
     call multifab_build(snew(n) ,mla%la(n),nscal,2)
     call multifab_build(prim(n) ,mla%la(n),nscal,2)

     ! s on faces, gp0 on faces
     do i=1,dm
        call multifab_build_edge(  s_fc(n,i),mla%la(n),nscal,1,i)
        call multifab_build_edge(gp0_fc(n,i),mla%la(n),1    ,0,i)
     end do

     ! pressure
     ! need 1 ghost cell since we calculate its gradient
     call multifab_build(pres(n),mla%la(n),1,1)

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
  call init(mold,sold,pres,gp0_fc,dx,mla,time)

  if (print_int .gt. 0) then
     call eos_check(mla,sold)
     call sum_mass_momentum(mla,sold,mold)
  end if

  ! convert cons to prim in valid region
  call convert_cons_to_prim(mla,sold,prim,.true.)

  ! fill ghost cells for prim
  do n=1,nlevs
     call multifab_fill_boundary(prim(n))
     call multifab_physbc(prim(n),1,scal_bc_comp,2,the_bc_tower%bc_tower_array(n),dx(n,:))
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
  call init_stochastic(mla)

  ! fill the stochastic multifabs with a new set of random numbers
  call fill_stochastic(mla)  

  ! need to do an initial projection to get an initial velocity field
  call initial_projection(mla,mold,umac,sold,s_fc,prim,chi_fc,gp0_fc,rhoc_d_fluxdiv, &
                          rhoc_s_fluxdiv,rhoc_b_fluxdiv,dx,the_bc_tower,vel_bc_n,vel_bc_t)

  if (print_int .gt. 0) then
     call sum_mass_momentum(mla,sold,mold)
  end if

  ! write initial plotfile
  if (plot_int .gt. 0) then
     call write_plotfile(mla,mold,umac,sold,pres,dx,time,0)
  end if
  
  do istep=1,max_step

     if (parallel_IOProcessor()) then
        print*,"Begin Advance; istep =",istep,"DT =",fixed_dt,"TIME =",time
     end if

     ! advance the solution by dt
     call advance_timestep(mla,mold,mnew,umac,sold,snew,s_fc,prim,pres,chi,chi_fc, &
                           eta,eta_ed,kappa,rhoc_d_fluxdiv,rhoc_s_fluxdiv,rhoc_b_fluxdiv, &
                           gp0_fc,dx,the_bc_tower,vel_bc_n,vel_bc_t)

     ! increment simulation time
     time = time + fixed_dt

    if (parallel_IOProcessor()) then
        print*,"End Advance; istep =",istep,"DT =",fixed_dt,"TIME =",time
     end if

     ! project rho and rho1 back onto EOS
     if ( project_eos_int .gt. 0 .and. mod(istep,project_eos_int) .eq. 0) then
        call project_onto_eos(mla,snew)
     end if

     ! write a plotfile
     if ( (plot_int .gt. 0 .and. mod(istep,plot_int) .eq. 0) &
          .or. &
          (istep .eq. max_step) ) then
        call write_plotfile(mla,mnew,umac,snew,pres,dx,time,istep)
     end if
     
     if ( (print_int .gt. 0 .and. mod(istep,print_int) .eq. 0) &
          .or. &
          (istep .eq. max_step) ) then
        call eos_check(mla,snew)
        call sum_mass_momentum(mla,snew,mnew)
     end if

     ! set old state to new state
     do n=1,nlevs
        call multifab_copy_c(sold(n),1,snew(n),1,nscal,sold(n)%ng)
        do i=1,dm
           call multifab_copy_c(mold(n,i),1,mnew(n,i),1,1,mold(n,i)%ng)
        end do
     end do
        
  end do

  ! destroy and deallocate multifabs that contain random numbers
  call destroy_stochastic(mla)

  do n=1,nlevs
     call multifab_destroy(sold(n))
     call multifab_destroy(snew(n))
     call multifab_destroy(prim(n))
     call multifab_destroy(pres(n))
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
        call multifab_destroy(gp0_fc(n,i))
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
