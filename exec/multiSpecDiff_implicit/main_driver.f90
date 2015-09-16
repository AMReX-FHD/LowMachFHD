subroutine main_driver()

  use boxlib
  use bl_IO_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use bndry_reg_module
  use multifab_physbc_module
  use init_n_module
  use stochastic_n_fluxdiv_module
  use diffusive_n_fluxdiv_module
  use write_plotfile_n_module
  use ml_solve_module
  use ParallelRNGs 
  use probin_common_module, only: prob_lo, prob_hi, n_cells, dim_in, max_grid_size, &
                                  plot_int, chk_int, print_int, seed, bc_lo, bc_hi, restart, &
                                  probin_common_init, fixed_dt, max_step
  use probin_multispecies_module, only: nspecies, start_time, probin_multispecies_init
  use probin_gmres_module, only: probin_gmres_init, mg_verbose, cg_verbose

  use fabio_module

  implicit none

  ! quantities will be allocated with dm components
  integer, allocatable :: lo(:), hi(:)

  ! quantities will be allocated with (nlevs,dm) components
  real(kind=dp_t), allocatable :: dx(:,:)
  real(kind=dp_t)              :: dt,time,runtime1,runtime2
  integer                      :: n,nlevs,i,dm,istep,ng_s,init_step,comp
  type(box)                    :: bx
  type(ml_boxarray)            :: mba
  type(ml_layout)              :: mla
  type(bc_tower)               :: the_bc_tower
  logical, allocatable         :: pmask(:)
  
  type(multifab), allocatable :: n_old(:)
  type(multifab), allocatable :: n_new(:)
  type(multifab), allocatable :: diff_fluxdiv(:)
  type(multifab), allocatable :: stoch_fluxdiv(:)
  type(multifab), allocatable :: diff_coef_face(:,:)

  ! for multigrid solver; (alpha - div beta grad) phi = rhs
  type(multifab), allocatable :: alpha(:)
  type(multifab), allocatable :: rhs(:)
  type(multifab), allocatable :: phi(:)
  type(multifab), allocatable :: beta(:,:)

  ! for diffusion multigrid - not used but needs to be passed in
  type(bndry_reg), allocatable :: fine_flx(:)

  ! For HydroGrid
  integer :: n_rngs

  ! to test "conservation"
  real(kind=dp_t), allocatable :: n_sum(:)
  
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
  allocate(n_old(nlevs),n_new(nlevs))
  allocate(diff_fluxdiv(nlevs),stoch_fluxdiv(nlevs))
  allocate(diff_coef_face(nlevs,dm))
  allocate(alpha(nlevs),rhs(nlevs),phi(nlevs),beta(nlevs,dm))
  allocate(fine_flx(2:mla%nlevel))

  allocate(n_sum(nspecies))

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

  ! use refined dx for next level
  do n=2,nlevs
     dx(n,:) = dx(n-1,:) / mba%rr(n-1,:)
  end do

  ! build pmask
  allocate(pmask(dm))
  pmask = .false.
  do i=1,dm
     if (bc_lo(i) .eq. PERIODIC .and. bc_hi(i) .eq. PERIODIC) then
        pmask(i) = .true.
     end if
  end do

  ng_s = 2

  if (restart .ge. 0) then

     call bl_error("restart not supported yet")

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
        call multifab_build(n_new(n)        ,mla%la(n),nspecies,ng_s)
        call multifab_build(n_old(n)        ,mla%la(n),nspecies,ng_s) 
        call multifab_build(diff_fluxdiv(n) ,mla%la(n),nspecies,0) 
        call multifab_build(stoch_fluxdiv(n),mla%la(n),nspecies,0) 
        do i=1,dm
           call multifab_build_edge(diff_coef_face(n,i),mla%la(n),nspecies,0,i)
        end do
        ! for multigrid solver; (alpha - div beta grad) phi = rhs
        call multifab_build(alpha(n),mla%la(n),1,0)
        call multifab_build(rhs(n)  ,mla%la(n),1,0)
        call multifab_build(phi(n)  ,mla%la(n),1,1) 
        do i=1,dm
           call multifab_build_edge(beta(n,i),mla%la(n),1,0,i)
        end do
     end do

     ! stores beta*grad phi/dx_fine on coarse-fine interfaces
     ! this gets computed inside of ml_cc_solve
     ! we pass it back out because some algorithms (like projection methods) 
     ! use this information
     do n = 2,nlevs
        call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
     end do

  end if

  deallocate(pmask)

  ! allocate and build multifabs that will contain random numbers
  n_rngs = 1
  call init_mass_stochastic(mla,n_rngs)

  !=======================================================
  ! Setup boundary condition bc_tower
  !=======================================================
 
  ! bc_tower structure in memory
  ! 1:dm = velocity
  ! dm+1 = pressure
  ! dm+2 = scal_bc_comp (for n_i)
  ! scal_bc_comp+nspecies = tran_bc_comp = diffusion coefficients
  ! It may be better if each transport coefficient has its own BC code?
  ! I think the only place this is used is average_cc_to_node/face/edge
  ! I cannot right now foresee a case where different values would be used in different places
  ! so it is OK to keep num_tran_bc_in=1. But note the same code applies to eta,kappa and chi's
  call initialize_bc(the_bc_tower,nlevs,dm,mla%pmask, &
                     num_scal_bc_in=nspecies, &
                     num_tran_bc_in=1)

  do n=1,nlevs
     ! define level n of the_bc_tower
     call bc_tower_level_build(the_bc_tower,n,mla%la(n))
  end do

  !=====================================================================
  ! Initialize values
  !=====================================================================

  call init_n(mla,n_old,dx)
  do n=1,nlevs
     call multifab_fill_boundary(n_old(n))
  end do

  if (restart .lt. 0) then

     if (fixed_dt .gt. 0.d0) then
        dt = fixed_dt
     else
        call bl_error("Need to define what we mean by CFL time step")
     end if
     
  end if

  ! write a plotfile
  if (plot_int .gt. 0) then
     call write_plotfile_n(mla,n_old,dx,0.d0,0)
  end if

  !=======================================================
  ! Begin time stepping loop
  !=======================================================

  do istep=init_step,max_step

     runtime1 = parallel_wtime()

     if (parallel_IOProcessor()) then
        if ( (print_int .gt. 0 .and. mod(istep,print_int) .eq. 0) ) &
           print*,"Begin Advance; istep =",istep,"dt =",dt,"time =",time
     end if

     do comp=1,nspecies
        n_sum(comp) = multifab_sum_c(n_old(1),comp,1)
     end do
     if (parallel_IOProcessor()) then
        if ( (print_int .gt. 0 .and. mod(istep,print_int) .eq. 0) ) &
             print*,'sum of n',n_sum(:)
     end if


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! advance the solution by dt
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! fill random flux multifabs with new random numbers
     call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)

     ! compute the diffusion coefficients (for now just setting each to a different constant)
     do n=1,nlevs
        do i=1,dm
           do comp=1,nspecies
              call multifab_setval_c(diff_coef_face(n,i),dble(comp),comp,1,all=.true.)
           end do
        end do
     end do

     ! compute diffusive flux diverge
     call diffusive_n_fluxdiv(mla,n_old,diff_coef_face,diff_fluxdiv,dx,the_bc_tower)

     ! compute stochastic flux divergence
     call stochastic_n_fluxdiv(mla,n_old,diff_coef_face,stoch_fluxdiv,dx,the_bc_tower)

     ! n_k^{n+1} = n_k^n + (dt/2)(div D_k grad n_k)^n
     !                   + (dt/2)(div D_k grad n_k)^n+1
     !                   +  dt    div (sqrt(2 D_k n_k) Z)^n
     ! 
     ! in operator form
     !
     ! (I - (dt/2) div D_k grad)n_k^{n+1} = n_k^n + (dt/2)(div D_k grad n_k)^n
     !                                            +  dt    div (sqrt(2 D_k n_k) Z)^n
     !
     
     do comp=1,nspecies

        ! alpha = 1
        ! beta = (dt/2)*D_k
        do n=1,nlevs
           call multifab_setval(alpha(n),1.d0,all=.true.)
           do i=1,dm
              call multifab_copy_c(beta(n,i),1,diff_coef_face(n,i),comp,1,0)
              call multifab_mult_mult_s(beta(n,i),0.5d0*fixed_dt)
           end do
        end do

        ! rhs = n_k^n + (dt/2)(div D_k grad n_k)^n
        !             +  dt    div (sqrt(2 D_k n_k) Z)^n
        do n=1,nlevs
           call multifab_copy_c(rhs(n),1,diff_fluxdiv(n),comp,1,0)
           call multifab_mult_mult_s(rhs(n),0.5d0)
           call multifab_plus_plus_c(rhs(n),1,stoch_fluxdiv(n),comp,1,0)
           call multifab_mult_mult_s(rhs(n),fixed_dt)
           call multifab_plus_plus_c(rhs(n),1,n_old(n),comp,1,0)
        end do

        ! initial guess for phi is n_k^n
        do n=1,nlevs
           call multifab_copy_c(phi(n),1,n_old(n),comp,1,1)
        end do

        ! solve the implicit system
        call ml_cc_solve(mla,rhs,phi,fine_flx,alpha,beta,dx, &
                         the_bc_tower,scal_bc_comp+comp-1, &
                         verbose=mg_verbose, &
                         cg_verbose=cg_verbose)

        ! copy solution into n_new
        do n=1,nlevs
           call multifab_copy_c(n_new(n),comp,phi(n),1,1,0)
        end do

     end do

     do n=1,nlevs
        call multifab_fill_boundary(n_new(n))
     end do

     time = time + dt

      if (parallel_IOProcessor()) then
        if (print_int .gt. 0 .and. mod(istep,print_int) .eq. 0) then
           print*,"End Advance; istep =",istep,"DT =",dt,"TIME =",time
        end if
      end if

      runtime2 = parallel_wtime() - runtime1
      call parallel_reduce(runtime1, runtime2, MPI_MAX, proc=parallel_IOProcessorNode())
      if (parallel_IOProcessor()) then
        if (print_int .gt. 0 .and. mod(istep,print_int) .eq. 0) then
           print*,'Time to advance timestep: ',runtime1,' seconds'
        end if
      end if

      ! write a plotfile
      if (plot_int .gt. 0 .and. mod(istep,plot_int) .eq. 0) then
         call write_plotfile_n(mla,n_new,dx,time,istep)
      end if

      ! write a checkpoint
      if (chk_int .gt. 0 .and. mod(init_step,chk_int) .eq. 0) then

      end if

      ! set old state to new state
      do n=1,nlevs
         call multifab_copy_c(n_old(n),1,n_new(n),1,nspecies,n_old(n)%ng)
      end do

  end do

  !=======================================================
  ! Destroy multifabs and layouts
  !=======================================================

  call destroy_mass_stochastic(mla)

  do n=1,nlevs
     call multifab_destroy(n_new(n))
     call multifab_destroy(n_old(n))
     call multifab_destroy(diff_fluxdiv(n))
     call multifab_destroy(stoch_fluxdiv(n))
     do i=1,dm
        call multifab_destroy(diff_coef_face(n,i))
     end do
     call multifab_destroy(alpha(n))
     call multifab_destroy(rhs(n))
     call multifab_destroy(phi(n))
     do i=1,dm
        call multifab_destroy(beta(n,i))
     end do
  end do

  do n = 2,nlevs
     call bndry_reg_destroy(fine_flx(n))
  end do

  call destroy(mla)
  call bc_tower_destroy(the_bc_tower)

end subroutine main_driver
