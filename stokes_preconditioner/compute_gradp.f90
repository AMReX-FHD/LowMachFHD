module compute_gradp_module

  use multifab_module
  use ml_layout_module
  use macproject_module

  implicit none

  private

  public :: compute_gradp

contains

  subroutine compute_gradp(mla,pres,gradp,dx)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: pres(:)
    type(multifab) , intent(inout) :: gradp(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)

    ! local
    integer :: n,i,dm,nlevs,ng_p,ng_g
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: pp(:,:,:,:)
    real(kind=dp_t), pointer :: gpx(:,:,:,:)
    real(kind=dp_t), pointer :: gpy(:,:,:,:)
    real(kind=dp_t), pointer :: gpz(:,:,:,:)

    dm = mla%dim
    nlevs = mla%nlevel

    ng_p = pres(1)%ng
    ng_g = gradp(1,1)%ng

    do n=1,nlevs
       do i=1,nfabs(pres(n))
          pp  => dataptr(pres(n), i)
          gpx => dataptr(gradp(n,1), i)
          gpy => dataptr(gradp(n,2), i)
          lo = lwb(get_box(pres(n), i))
          hi = upb(get_box(pres(n), i))
          select case (dm)
          case (2)
             call compute_gradp_2d(pp(:,:,1,1), ng_p, &
                                   gpx(:,:,1,1), gpy(:,:,1,1), ng_g, &
                                   lo, hi, dx(n,:))
          case (3)
             gpz => dataptr(gradp(n,3), i)
             call compute_gradp_3d(pp(:,:,:,1), ng_p, &
                                   gpx(:,:,:,1), gpy(:,:,:,1), gpz(:,:,:,1), ng_g, &
                                   lo, hi, dx(n,:))
          end select
       end do
    end do

  contains
    
    subroutine compute_gradp_2d(pres,ng_p,gpx,gpy,ng_g,lo,hi,dx)

      integer        , intent(in   ) :: ng_p,ng_g,lo(:),hi(:)
      real(kind=dp_t), intent(in   ) :: pres(lo(1)-ng_p:,lo(2)-ng_p:)
      real(kind=dp_t), intent(inout) ::  gpx(lo(1)-ng_g:,lo(2)-ng_g:)
      real(kind=dp_t), intent(inout) ::  gpy(lo(1)-ng_g:,lo(2)-ng_g:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      ! local
      integer :: i,j

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
            gpx(i,j) = ( pres(i,j)-pres(i-1,j) ) / dx(1)
         end do
      end do

      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)
            gpy(i,j) = ( pres(i,j)-pres(i,j-1) ) / dx(2)
         end do
      end do

    end subroutine compute_gradp_2d

    subroutine compute_gradp_3d(pres,ng_p,gpx,gpy,gpz,ng_g,lo,hi,dx)

      integer        , intent(in   ) :: ng_p,ng_g,lo(:),hi(:)
      real(kind=dp_t), intent(in   ) :: pres(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
      real(kind=dp_t), intent(inout) ::  gpx(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:)
      real(kind=dp_t), intent(inout) ::  gpy(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:)
      real(kind=dp_t), intent(inout) ::  gpz(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      ! local
      integer :: i,j,k

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)+1
               gpx(i,j,k) = ( pres(i,j,k)-pres(i-1,j,k) ) / dx(1)
            end do
         end do
      end do

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)
               gpy(i,j,k) = ( pres(i,j,k)-pres(i,j-1,k) ) / dx(2)
            end do
         end do
      end do

      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               gpz(i,j,k) = ( pres(i,j,k)-pres(i,j,k-1) ) / dx(3)
            end do
         end do
      end do

    end subroutine compute_gradp_3d

  end subroutine compute_gradp

end module compute_gradp_module
