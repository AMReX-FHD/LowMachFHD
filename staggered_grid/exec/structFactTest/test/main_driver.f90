subroutine main_driver()

  use bl_IO_module
  use ml_layout_module
  use analyze_spectra_module
  use probin_common_module, only: probin_common_init, prob_lo, prob_hi, n_cells, dim_in, max_grid_size

  use fabio_module

  implicit none

  ! quantities will be allocated with dm components
  integer, allocatable :: lo(:), hi(:)

  ! quantities will be allocated with (nlevs,dm) components
  real(kind=dp_t), allocatable :: dx(:,:)
  real(kind=dp_t)              :: dt,time
  integer                      :: n,nlevs,i,dm,step,id
  type(box)                    :: bx
  type(ml_boxarray)            :: mba
  type(ml_layout)              :: mla
  logical, allocatable         :: pmask(:)
  
  ! will be allocated on nlevels
  type(multifab), allocatable  :: rho(:)
  type(multifab), allocatable  :: Temp(:)
  type(multifab), allocatable  :: umac(:,:)

  real(kind=dp_t), pointer :: dp(:,:,:,:)
  
  ! For HydroGrid
  integer :: narg, farg, un
  character(len=128) :: fname
  logical :: lexist
  
  !==============================================================
  ! Initialization
  !==============================================================

  call probin_common_init()

  ! dimensionality is set in inputs file
  dm = dim_in

  ! in this example we fix nlevs to be 1
  nlevs = 1

  allocate(pmask(dm))
  pmask(:) = .true.
       
  ! tell mba how many levels and dimensionality of problem
  call ml_boxarray_build_n(mba,nlevs,dm)

  ! tell mba about the ref_ratio between levels
  ! mba%rr(n-1,i) is the refinement ratio between levels n-1 and n in direction i
  ! we use refinement ratio of 2 in every direction between all levels
  do n=2,nlevs
     mba%rr(n-1,:) = 2
  end do

  ! create a box from (0,0) to (n_cells-1,n_cells-1)
  allocate(lo(dm))
  allocate(hi(dm))
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

  ! 2 species hard coded
  allocate(rho(nlevs))
  do n=1,nlevs
     call multifab_build(rho(n),mla%la(n),2,0)
  end do
  ! set grid spacing at each level
  allocate(dx(nlevs,dm))
  dx(1,1:dm) = (prob_hi(1:dm)-prob_lo(1:dm)) / n_cells(1:dm)
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

  !=====================================================================
  ! Initialize rho
  !=====================================================================

  do n=1,nlevs
     do i=1,nfabs(rho(n))
        dp => dataptr(rho(n),i)
        lo = lwb(get_box(rho(n),i))
        hi = upb(get_box(rho(n),i))
        select case (dm)
        case (2)
        case (3)
           dp => dataptr(rho(n),i)
           call init_rho_3d(dp(:,:,:,:),lo,hi)
        end select
     end do
  end do
  
  call fabio_ml_multifab_write_d(rho,mla%mba%rr(:,1),"a_rho")
  
  !=====================================================================
  ! Initialize HydroGrid for analysis
  !=====================================================================
  narg = command_argument_count()
  farg = 1
  if (narg >= 1) then
     call get_command_argument(farg, value = fname)
     inquire(file = fname, exist = lexist )
     if ( lexist ) then
        un = unit_new()
        open(unit=un, file = fname, status = 'old', action = 'read')
        
        ! We will also pass temperature
        call initialize_hydro_grid(mla,rho,dt,dx,namelist_file=un, & 
                                   nspecies_in=2, &
                                   nscal_in=0, &
                                   exclude_last_species_in=.false., &
                                   analyze_velocity=.false., &
                                   analyze_density=.true., &
                                   analyze_temperature=.false.)
        
        close(unit=un)

     else

        call bl_error('HydroGrid initialization requires a namelist in an input file')

     end if
  end if

  dt = 1.d0
  step = 0
  id = 0
  
  ! Add this snapshot to the average in HydroGrid
  call analyze_hydro_grid(mla,dt,dx,step,rho=rho)

  ! increment id each time you call this
  ! step is for the output filename
  call save_hydro_grid(id,step)
  id = id+1
  step = step+1

  ! now set rho to zero
!  do n=1,nlevs
!     call setval(rho(n),0.d0)
!  end do

  ! Add this snapshot to the average in HydroGrid
!  call analyze_hydro_grid(mla,dt,dx,step,rho=rho)

  ! increment id each time you call this
  ! step is for the output filename
!  call save_hydro_grid(id,step)
!  id = id+1
!  step = step+1

  ! cleanup

  do n=1,nlevs
     call multifab_destroy(rho(n))
  end do
  
  call finalize_hydro_grid()
  call destroy(mla)

contains

  subroutine init_rho_3d(rho,lo,hi)

    integer :: lo(3), hi(3)
    real(kind=dp_t) :: rho(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:2)

    integer :: i,j,k

    do k=lo(3),hi(3)
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)

       rho(i,j,k,1) = (i+j+k)
       rho(i,j,k,2) = 0.5*sqrt(dble(i+j+k))
       
    end do
    end do
    end do

  end subroutine init_rho_3d  
  
end subroutine main_driver
