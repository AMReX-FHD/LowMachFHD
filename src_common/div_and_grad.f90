module div_and_grad_module

  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: compute_gradp, compute_divu

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

  subroutine compute_divu(mla,umac,rh,dx)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: umac(:,:)
    type(multifab) , intent(inout) :: rh(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)

    real(kind=dp_t), pointer :: ump(:,:,:,:) 
    real(kind=dp_t), pointer :: vmp(:,:,:,:) 
    real(kind=dp_t), pointer :: wmp(:,:,:,:) 
    real(kind=dp_t), pointer :: rhp(:,:,:,:) 
    real(kind=dp_t)          :: rhmax
    integer :: i,n,nlevs,dm,ng_u,ng_r,lo(mla%dim),hi(mla%dim)

    dm = mla%dim
    nlevs = mla%nlevel

    ng_u = umac(1,1)%ng
    ng_r = rh(1)%ng

    do n = 1,nlevs
       do i = 1, nfabs(rh(n))
          ump => dataptr(umac(n,1), i)
          vmp => dataptr(umac(n,2), i)
          rhp => dataptr(rh(n)  , i)
          lo =  lwb(get_box(rh(n), i))
          hi =  upb(get_box(rh(n), i))
          select case (dm)
          case (2)
             call compute_divu_2d(ump(:,:,1,1), vmp(:,:,1,1), ng_u, rhp(:,:,1,1), ng_r, &
                                  dx(n,:),lo,hi)
          case (3)
             wmp => dataptr(umac(n,3), i)
             call compute_divu_3d(ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_u, &
                                  rhp(:,:,:,1), ng_r, dx(n,:),lo,hi)
          end select
       end do
    end do

    rhmax = norm_inf(rh(nlevs))
    do n = nlevs,2,-1
       rhmax = max(rhmax,norm_inf(rh(n-1)))
    end do

  end subroutine compute_divu

  subroutine compute_divu_2d(umac,vmac,ng_u,rh,ng_r,dx,lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:),ng_u,ng_r
    real(kind=dp_t), intent(in   ) :: umac(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t), intent(in   ) :: vmac(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t), intent(inout) ::   rh(lo(1)-ng_r:,lo(2)-ng_r:)
    real(kind=dp_t), intent(in   ) ::   dx(:)

    integer :: i,j

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          rh(i,j) = &
               (umac(i+1,j) - umac(i,j)) / dx(1) + &
               (vmac(i,j+1) - vmac(i,j)) / dx(2)
       end do
    end do

  end subroutine compute_divu_2d

  subroutine compute_divu_3d(umac,vmac,wmac,ng_u,rh,ng_r,dx,lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:),ng_u,ng_r
    real(kind=dp_t), intent(in   ) :: umac(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(in   ) :: vmac(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(in   ) :: wmac(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(inout) ::   rh(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    integer :: i,j,k

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             rh(i,j,k) = &
                  (umac(i+1,j,k) - umac(i,j,k)) / dx(1) + &
                  (vmac(i,j+1,k) - vmac(i,j,k)) / dx(2) + &
                  (wmac(i,j,k+1) - wmac(i,j,k)) / dx(3)
          end do
       end do
    end do

  end subroutine compute_divu_3d

end module div_and_grad_module
