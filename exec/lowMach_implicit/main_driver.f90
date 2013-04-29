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
  use probin_lowmach_module, only: probin_lowmach_init, max_step, nscal
  use probin_common_module , only: probin_common_init, seed, dim_in, n_cells, &
                                   prob_lo, prob_hi, max_grid_size, &
                                   bc_lo, bc_hi, fixed_dt
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

  type(multifab), allocatable :: mold(:,:) ! face-based
  type(multifab), allocatable :: mnew(:,:) ! face-based
  type(multifab), allocatable :: umac(:,:) ! face-based
  type(multifab), allocatable :: sold(:)   ! cell-centered
  type(multifab), allocatable :: snew(:)   ! cell-centered
  type(multifab), allocatable :: chi(:)    ! cell-centered
  type(multifab), allocatable :: eta(:)    ! cell-centered

  ! uncomment this once lowMach_implicit/probin.f90 is written
  call probin_lowmach_init()
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
  allocate(mold(nlevs,dm),mnew(nlevs,dm),umac(nlevs,dm))
  allocate(sold(nlevs),snew(nlevs),chi(nlevs),eta(nlevs))

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
  call initialize_bc(the_bc_tower,nlevs,dm,mla%pmask)
  do n=1,nlevs
     ! define level n of the_bc_tower
     call bc_tower_level_build(the_bc_tower,n,mla%la(n))
  end do

  do n=1,nlevs
     do i=1,dm
        ! edge-momentum and velocity; 1 ghost cell
        call multifab_build_edge(mold(n,i),mla%la(n),1,1,i)
        call multifab_build_edge(mnew(n,i),mla%la(n),1,1,i)
        call multifab_build_edge(umac(n,i),mla%la(n),1,1,i)
     end do
     ! 2 components (rho,rho1) and 2 ghost cells
     call multifab_build(sold(n),mla%la(n),nscal,2)
     call multifab_build(snew(n),mla%la(n),nscal,2)
     call multifab_build(chi(n) ,mla%la(n),nscal,1)
     call multifab_build(eta(n) ,mla%la(n),nscal,1)
  end do

  time = 0.d0

  ! initialize sold and mold
  call init(mold,sold,dx,mla,time)

  ! initialize eta and chi - for now just use a setval
  do n=1,nlevs
     call setval(eta(n),1.d0,all=.true.)
     call setval(chi(n),1.d0,all=.true.)
  end do

  ! need to do an initial projection to get an initial velocity field
  call initial_projection(mla,mold,umac,sold,chi,dx,the_bc_tower)

  ! write initial plotfile
  call write_plotfile(mla,mold,umac,sold,dx,time,0)
  
  do istep=1,max_step


     time = time + fixed_dt

  end do


  do n=1,nlevs
     do i=1,dm
        call destroy(mold(n,i))
        call destroy(mnew(n,i))
        call destroy(umac(n,i))
     end do
     call destroy(sold(n))
     call destroy(snew(n))
     call destroy(chi(n))
     call destroy(eta(n))
  end do

  call destroy(mla)
  call bc_tower_destroy(the_bc_tower)
  deallocate(lo,hi,dx,mold,mnew,sold,snew)

end subroutine main_driver
