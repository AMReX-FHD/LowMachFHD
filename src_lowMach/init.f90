module init_module

  use bl_types
  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: init

contains

  ! Important note: For periodic boundaries, the init routines should fill out
  ! *both* sides of the domain with values, even though this is duplicate information
  ! We ensure the two sides are bitwise identical, but low and/or high sides may win

  subroutine init(m,s,dx,mla,time)

    type(multifab) , intent(inout) :: m(:,:),s(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(ml_layout), intent(in   ) :: mla
    real(kind=dp_t), intent(in   ) :: time

    real(kind=dp_t), pointer :: mxp(:,:,:,:), myp(:,:,:,:), mzp(:,:,:,:)
    real(kind=dp_t), pointer :: sop(:,:,:,:)
    integer :: lo(mla%dim),hi(mla%dim)
    integer :: nlevs,n,i,ng_m,ng_s,dm

    nlevs = mla%nlevel
    dm = mla%dim
    
    ng_m = m(1,1)%ng
    ng_s = s(1)%ng

    do n=1,nlevs
       do i = 1, nfabs(s(n))
          mxp => dataptr(m(n,1),i)
          myp => dataptr(m(n,2),i)
          sop => dataptr(s(n),i)
          lo =  lwb(get_box(s(n),i))
          hi =  upb(get_box(s(n),i))
          select case (dm)
          case (2)
             call init_2d(mxp(:,:,1,1), myp(:,:,1,1), sop(:,:,1,:), &
                          lo, hi, ng_m, ng_s, dx(n,:), time)
          case (3)
             mzp => dataptr(m(n,3),i)
             call init_3d(mxp(:,:,:,1), myp(:,:,:,1), mzp(:,:,:,1), &
                          sop(:,:,:,:), &
                          lo, hi, ng_m, ng_s, dx(n,:), time)
          end select
       end do

       ! For periodic boundaries, ensure the low and high side are consistent:
       ! Note: multifab_internal_sync compares the box number of the two boxes 
       ! with overlapping values and the data on the box with lower number wins. 
       do i=1,dm
          call multifab_internal_sync(m(n,i))
          call multifab_fill_boundary(m(n,i))
       end do
       call multifab_fill_boundary(s(n))

    enddo

  end subroutine init

  subroutine init_2d(mx,my,s,lo,hi,ng_m,ng_s,dx,time)

    use exact_solutions_module, only: exact_2d

    integer, intent(in) :: lo(:), hi(:), ng_m, ng_s
    real (kind = dp_t), intent(out) :: mx(lo(1)-ng_m:,lo(2)-ng_m:)
    real (kind = dp_t), intent(out) :: my(lo(1)-ng_m:,lo(2)-ng_m:)
    real (kind = dp_t), intent(out) ::  s(lo(1)-ng_s:,lo(2)-ng_s:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)
    real (kind = dp_t), intent(in ) :: time

    call exact_2d(mx,my,s,lo,hi,ng_m,ng_s,dx,time)

  end subroutine init_2d

  subroutine init_3d(mx,my,mz,s,lo,hi,ng_m,ng_s,dx,time)

    use exact_solutions_module, only: exact_3d

    integer, intent(in) :: lo(:), hi(:), ng_m, ng_s
    real (kind = dp_t), intent(out) :: mx(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real (kind = dp_t), intent(out) :: my(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real (kind = dp_t), intent(out) :: mz(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real (kind = dp_t), intent(out) ::  s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind = dp_t), intent(in ) :: dx(:)
    real (kind = dp_t), intent(in ) :: time

    call exact_3d(mx,my,mz,s,lo,hi,ng_m,ng_s,dx,time)

  end subroutine init_3d

end module init_module
