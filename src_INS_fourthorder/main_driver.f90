subroutine main_driver()

  use ml_layout_module
  use define_bc_module
  use bc_module
  use init_module
  use stag_mg_layout_module
  use macproject_module
  use write_plotfile_module
  use probin_common_module, only: prob_lo, prob_hi, n_cells, dim_in, max_grid_size, print_int, &
                                  plot_int, bc_lo, bc_hi, restart, fixed_dt, max_step, &
                                  probin_common_init
  use probin_gmres_module, only: probin_gmres_init

  implicit none

  ! quantities will be allocated with dm components
  integer, allocatable :: lo(:), hi(:)

  ! quantities will be allocated with (nlevs,dm) components
  real(kind=dp_t), allocatable :: dx(:,:)
  real(kind=dp_t)              :: dt,time,runtime1,runtime2
  integer                      :: n,nlevs,i,dm,istep,init_step
  type(box)                    :: bx
  type(ml_boxarray)            :: mba
  type(ml_layout)              :: mla
  type(bc_tower)               :: the_bc_tower
  logical, allocatable         :: pmask(:)
  
  ! will be allocated on nlevels
  type(multifab), allocatable :: pres(:)
  type(multifab), allocatable :: umac(:,:)
  type(multifab), allocatable :: eta(:)
  type(multifab), allocatable :: eta_ed(:,:)
  type(multifab), allocatable :: kappa(:)

  logical :: nodal_temp(3)
  
  !==============================================================
  ! Initialization
  !==============================================================

  call probin_common_init()
  call probin_gmres_init()

  ! in this example we fix nlevs to be 1
  ! for adaptive simulations where the grids change, cells at finer
  ! resolution don't necessarily exist depending on your tagging criteria,
  ! so max_levs isn't necessary equal to nlevs
  nlevs = 1

  ! dimensionality is set in inputs file
  dm = dim_in
 
  ! now that we have dm, we can allocate these
  allocate(lo(dm),hi(dm))
  allocate(pres(nlevs))
  allocate(umac(nlevs,dm))
  allocate(eta(nlevs),kappa(nlevs))
  if (dm .eq. 2) then
     allocate(eta_ed(nlevs,1))
  else if (dm .eq. 3) then
     allocate(eta_ed(nlevs,3))
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

     ! initialize from restart

  else

     init_step = 1
     time = 0.d0
     
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
        call multifab_build(pres(n),mla%la(n),1,2)
        do i=1,dm
           call multifab_build_edge(umac(n,i),mla%la(n),1,2,i)
        end do
     end do

  end if

  deallocate(pmask)

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

  !=======================================================
  ! Setup boundary condition bc_tower
 
  ! bc_tower structure in memory
  ! 1:dm = velocity
  ! dm+1 = pressure
  call initialize_bc(the_bc_tower,nlevs,dm,mla%pmask, &
                     num_scal_bc_in=0, &
                     num_tran_bc_in=0)

  do n=1,nlevs
     ! define level n of the_bc_tower
     call bc_tower_level_build(the_bc_tower,n,mla%la(n))
  end do

  ! build layouts for staggered multigrid solver and macproject within preconditioner
  call stag_mg_layout_build(mla)
  call mgt_macproj_precon_build(mla,dx,the_bc_tower)

  ! initialize umac
  call init(mla,umac,dx,time)

  ! write initial plotfile
  call write_plotfile(mla,"plt",umac,pres,0,dx,0.d0)

  print*,'here'
  stop

  !=======================================================
  ! Build multifabs for all the variables
  !=======================================================

  ! build multifab with nspecies component and one ghost cell
  do n=1,nlevs 
     call multifab_build(eta(n)  ,mla%la(n),1,1)
     call multifab_build(kappa(n),mla%la(n),1,1)

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

  !=====================================================================
  ! Initialize values
  !=====================================================================

  ! initialize eta and kappa




  ! compute dt
  dt = fixed_dt

  ! write initial plotfile



  !=======================================================
  ! Begin time stepping loop
  !=======================================================

  do istep=init_step,max_step

     if (parallel_IOProcessor()) then
        if ( (print_int .gt. 0 .and. mod(istep,print_int) .eq. 0) ) &
           print*,"Begin Advance; istep =",istep,"dt =",dt,"time =",time
     end if

     runtime1 = parallel_wtime()

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! advance the solution by dt
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     
     
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

     ! write plotfile


     ! write checkpoint


  end do

  !=======================================================
  ! Destroy multifabs and layouts
  !=======================================================

  do n=1,nlevs
     call multifab_destroy(eta(n))
     call multifab_destroy(kappa(n))
     call multifab_destroy(pres(n))
     do i=1,dm
        call multifab_destroy(umac(n,i))
     end do
     do i=1,size(eta_ed,dim=2)
        call multifab_destroy(eta_ed(n,i))
     end do
  end do
  deallocate(lo,hi,dx)
  deallocate(umac,eta,eta_ed,kappa)

  call stag_mg_layout_destroy()
  call mgt_macproj_precon_destroy()

  call destroy(mla)
  call bc_tower_destroy(the_bc_tower)

end subroutine main_driver
