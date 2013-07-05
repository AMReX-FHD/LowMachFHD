subroutine main_driver()

  use boxlib
  use multifab_module
  use bl_IO_module
  use layout_module
  use init_module
  use write_plotfile_module
  use advance_module
  use define_bc_module
  use bc_module
  use probin_common_module
  use probin_multispecies_module

  implicit none

  ! will be allocated with dm components
  integer, allocatable :: lo(:), hi(:)

  ! will be allocated with (nlevs,dm) components
  real(dp_t), allocatable :: dx(:,:)
  real(kind=dp_t) :: time
  integer :: n,nlevs,i,dm,istep

  type(box)         :: bx
  type(ml_boxarray) :: mba
  type(ml_layout)   :: mla
  type(bc_tower),  intent(in ) :: the_bc_tower
  type(bc_level),  intent(in ) :: the_bc_level(:)
  logical, allocatable :: pmask(:)
  type(multifab), allocatable :: rho(:)

  !=======================================================
  ! Initialization
  !=======================================================

  call probin_common_init()
  call probin_multispecies_init() 

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
  allocate(rho(nlevs))

  !=======================================================
  ! Setup parallelization: Create boxes and layouts for multifabs
  !=======================================================
  
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

  !=======================================================
  ! Setup boundary condition bc_tower
  !=======================================================

  ! tell the_bc_tower about max_levs, dm, and domain_phys_bc
  ! Donev: Last argument to initialize_bc is the number of scalar variables
  ! nscal=nspecies ! Number of scalars (maybe add temperature later)
  !call initialize_bc(the_bc_tower,nlevs,dm,mla%pmask,nscal)
  call initialize_bc(the_bc_tower,nlevs,dm,mla%pmask,nspecies)

  do n=1,nlevs
     ! define level n of the_bc_tower
     call bc_tower_level_build(the_bc_tower,n,mla%la(n))
  end do

  !=======================================================
  ! Build multifabs for all the variables
  !=======================================================

  ! build multifab with nspecies component and one ghost cell
  do n=1,nlevs
     call multifab_build(rho(n),mla%la(n),nspecies,1)
  end do

  call init_rho(rho,dx,prob_lo,the_bc_tower)

  !=======================================================
  ! Begin time stepping loop
  !=======================================================

  istep = 0
  time = 0.d0

  ! write initial plotfile
  if (plot_int .gt. 0) then
     call write_plotfile(mla,rho,istep,dx,time,prob_lo,prob_hi)
  end if
  
  do istep=1,max_step

     if (parallel_IOProcessor()) then
        print*,"Begin Advance; istep =",istep,"DT =",fixed_dt,"TIME =",time
     end if

     ! advance the solution by dt
     call advance(rho,dx,fixed_dt,the_bc_level)

     ! increment simulation time
     time = time + fixed_dt

     ! write a plotfile
     if ( (plot_int .gt. 0 .and. mod(istep,plot_int) .eq. 0) &
          .or. &
          (istep .eq. max_step) ) then
        call write_plotfile(mla,rho,istep,dx,time,prob_lo,prob_hi)
     end if
        
  end do

  !=======================================================
  ! Destroy multifabs and layouts
  !=======================================================

  do n=1,nlevs
     call multifab_destroy(rho(n))
  end do

  call destroy(mla)
  call bc_tower_destroy(the_bc_tower)

end subroutine main_driver
