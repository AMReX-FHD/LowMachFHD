module mk_diffusive_fluxdiv_module

  use ml_layout_module
  use define_bc_module
  use bc_module

  use fabio_module

  implicit none

  private

  public :: mk_diffusive_rhoc_fluxdiv

contains

  subroutine mk_diffusive_rhoc_fluxdiv(mla,s_update,out_comp,prim,rho_face,chi_face, &
                                       dx,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: s_update(:)
    integer        , intent(in   ) :: out_comp      ! which component of s_update
    type(multifab) , intent(in   ) ::     prim(:)   ! rho and c
    type(multifab) , intent(in   ) :: rho_face(:,:) ! rho on faces
    type(multifab) , intent(in   ) :: chi_face(:,:) ! chi on faces
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local variables
    integer :: nlevs,dm,i,n,ng_u,ng_p,ng_s,ng_c
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: pp(:,:,:,:)
    real(kind=dp_t), pointer :: spx(:,:,:,:)
    real(kind=dp_t), pointer :: spy(:,:,:,:)
    real(kind=dp_t), pointer :: spz(:,:,:,:)
    real(kind=dp_t), pointer :: cpx(:,:,:,:)
    real(kind=dp_t), pointer :: cpy(:,:,:,:)
    real(kind=dp_t), pointer :: cpz(:,:,:,:)

    ng_u = s_update(1)%ng
    ng_p = prim(1)%ng
    ng_s = rho_face(1,1)%ng
    ng_c = chi_face(1,1)%ng

    nlevs = mla%nlevel
    dm    = mla%dim

    ! compute del dot (rhoD grad c) and add it to s_update
    do n=1,nlevs
       do i=1,nfabs(prim(n))
          up  => dataptr(s_update(n), i)
          pp  => dataptr(prim(n), i)
          spx => dataptr(rho_face(n,1), i)
          spy => dataptr(rho_face(n,2), i)
          cpx => dataptr(chi_face(n,1), i)
          cpy => dataptr(chi_face(n,2), i)
          lo = lwb(get_box(prim(n), i))
          hi = upb(get_box(prim(n), i))
          select case (dm)
          case (2)
             call mk_diffusive_rhoc_fluxdiv_2d(up(:,:,1,out_comp), ng_u, pp(:,:,1,2), ng_p, &
                                               spx(:,:,1,1), spy(:,:,1,1), ng_s, &
                                               cpx(:,:,1,1), cpy(:,:,1,1), ng_c, &
                                               lo, hi, dx(n,:), &
                                               the_bc_level(n)%adv_bc_level_array(i,:,:,:))
          case (3)
             spz => dataptr(rho_face(n,3), i)
             cpz => dataptr(chi_face(n,3), i)
             call mk_diffusive_rhoc_fluxdiv_3d(up(:,:,:,out_comp), ng_u, pp(:,:,:,2), ng_p, &
                                               spx(:,:,:,1), spy(:,:,:,1), spz(:,:,:,1), ng_s, &
                                               cpx(:,:,:,1), cpy(:,:,:,1), cpz(:,:,:,1), ng_c, &
                                               lo, hi, dx(n,:), &
                                               the_bc_level(n)%adv_bc_level_array(i,:,:,:))
          end select
       end do
    end do

  contains

    subroutine mk_diffusive_rhoc_fluxdiv_2d(s_update,ng_u,c,ng_p,rhox,rhoy,ng_s, &
                                            chix,chiy,ng_c,lo,hi,dx,adv_bc)

      integer        , intent(in   ) :: lo(:), hi(:), ng_u, ng_p, ng_s, ng_c
      real(kind=dp_t), intent(inout) :: s_update(lo(1)-ng_u:,lo(2)-ng_u:)
      real(kind=dp_t), intent(in   ) ::        c(lo(1)-ng_p:,lo(2)-ng_p:)
      real(kind=dp_t), intent(in   ) ::     rhox(lo(1)-ng_s:,lo(2)-ng_s:)
      real(kind=dp_t), intent(in   ) ::     rhoy(lo(1)-ng_s:,lo(2)-ng_s:)
      real(kind=dp_t), intent(in   ) ::     chix(lo(1)-ng_c:,lo(2)-ng_c:)
      real(kind=dp_t), intent(in   ) ::     chiy(lo(1)-ng_c:,lo(2)-ng_c:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: adv_bc(:,:,:)

      ! local
      integer :: i,j

      real(kind=dp_t) :: fluxx(lo(1):hi(1)+1,lo(2):hi(2))
      real(kind=dp_t) :: fluxy(lo(1):hi(1),lo(2):hi(2)+1)

      ! x-faces
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
            ! fluxx = rhoD * grad c
            fluxx(i,j) = rhox(i,j)*chix(i,j) * (c(i,j)-c(i-1,j)) / dx(1)
         end do
      end do

      ! alter grad(c) stencil at boundary since ghost value represents value at boundary
      if (adv_bc(1,1,4) .eq. FOEXTRAP .or. adv_bc(1,1,4) .eq. EXT_DIR) then
         i=lo(1)
         do j=lo(2),hi(2)
            fluxx(i,j) = rhox(i,j)*chix(i,j) * (c(i,j)-c(i-1,j)) / (0.5d0*dx(1))
         end do
      end if
      if (adv_bc(1,2,4) .eq. FOEXTRAP .or. adv_bc(1,2,4) .eq. EXT_DIR) then
         i=hi(1)+1
         do j=lo(2),hi(2)
            fluxx(i,j) = rhox(i,j)*chix(i,j) * (c(i,j)-c(i-1,j)) / (0.5d0*dx(1))
         end do
      end if

      ! y-faces
      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)
            ! fluxy = rhoD * grad c
            fluxy(i,j) = rhoy(i,j)*chiy(i,j) * (c(i,j)-c(i,j-1)) / dx(2)
         end do
      end do

      ! alter grad(c) stencil at boundary since ghost value represents value at boundary
      if (adv_bc(2,1,4) .eq. FOEXTRAP .or. adv_bc(2,1,4) .eq. EXT_DIR) then
         j=lo(2)
         do i=lo(1),hi(1)
            fluxy(i,j) = rhoy(i,j)*chiy(i,j) * (c(i,j)-c(i,j-1)) / (0.5d0*dx(2))
         end do
      end if
      if (adv_bc(2,2,4) .eq. FOEXTRAP .or. adv_bc(2,2,4) .eq. EXT_DIR) then
         j=hi(2)+1
         do i=lo(1),hi(1)
            fluxy(i,j) = rhoy(i,j)*chiy(i,j) * (c(i,j)-c(i,j-1)) / (0.5d0*dx(2))
         end do
      end if

      ! flux divergence
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            ! add div(rhoD * grad c) to s_update
            s_update(i,j) = s_update(i,j) + (fluxx(i+1,j) - fluxx(i,j)) / dx(1) &
                                          + (fluxy(i,j+1) - fluxy(i,j)) / dx(2)
         end do
      end do

    end subroutine mk_diffusive_rhoc_fluxdiv_2d

    subroutine mk_diffusive_rhoc_fluxdiv_3d(s_update,ng_u,c,ng_p,rhox,rhoy,rhoz,ng_s, &
                                            chix,chiy,chiz,ng_c,lo,hi,dx,adv_bc)

      integer        , intent(in   ) :: lo(:), hi(:), ng_u, ng_p, ng_s, ng_c
      real(kind=dp_t), intent(inout) :: s_update(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
      real(kind=dp_t), intent(in   ) ::        c(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
      real(kind=dp_t), intent(in   ) ::     rhox(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
      real(kind=dp_t), intent(in   ) ::     rhoy(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
      real(kind=dp_t), intent(in   ) ::     rhoz(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
      real(kind=dp_t), intent(in   ) ::     chix(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:)
      real(kind=dp_t), intent(in   ) ::     chiy(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:)
      real(kind=dp_t), intent(in   ) ::     chiz(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: adv_bc(:,:,:)

      ! local
      integer :: i,j,k

      real(kind=dp_t) :: fluxx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3))
      real(kind=dp_t) :: fluxy(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3))
      real(kind=dp_t) :: fluxz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1)

      ! x-faces
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)+1
               ! fluxx = rhoD * grad c
               fluxx(i,j,k) = rhox(i,j,k)*chix(i,j,k) * (c(i,j,k) - c(i-1,j,k)) / dx(1)
            end do
         end do
      end do

      ! alter grad(c) stencil at boundary since ghost value represents value at boundary
      if (adv_bc(1,1,5) .eq. FOEXTRAP .or. adv_bc(1,1,5) .eq. EXT_DIR) then
         i=lo(1)
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               fluxx(i,j,k) = rhox(i,j,k)*chix(i,j,k) * (c(i,j,k)-c(i-1,j,k)) / (0.5d0*dx(1))
            end do
         end do
      end if
      if (adv_bc(1,2,5) .eq. FOEXTRAP .or. adv_bc(1,2,5) .eq. EXT_DIR) then
         i=hi(1)+1
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               fluxx(i,j,k) = rhox(i,j,k)*chix(i,j,k) * (c(i,j,k)-c(i-1,j,k)) / (0.5d0*dx(1))
            end do
         end do
      end if

      ! y-faces
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)
               ! fluxy = rhoD * grad c
               fluxy(i,j,k) = rhoy(i,j,k)*chiy(i,j,k) * (c(i,j,k) - c(i,j-1,k)) / dx(2)
            end do
         end do
      end do

      ! alter grad(c) stencil at boundary since ghost value represents value at boundary
      if (adv_bc(2,1,5) .eq. FOEXTRAP .or. adv_bc(2,1,5) .eq. EXT_DIR) then
         j=lo(2)
         do k=lo(3),hi(3)
            do i=lo(1),hi(1)
               fluxy(i,j,k) = rhoy(i,j,k)*chiy(i,j,k) * (c(i,j,k)-c(i,j-1,k)) / (0.5d0*dx(2))
            end do
         end do
      end if
      if (adv_bc(2,2,5) .eq. FOEXTRAP .or. adv_bc(2,2,5) .eq. EXT_DIR) then
         j=hi(2)+1
         do k=lo(3),hi(3)
            do i=lo(1),hi(1)
               fluxy(i,j,k) = rhoy(i,j,k)*chiy(i,j,k) * (c(i,j,k)-c(i,j-1,k)) / (0.5d0*dx(2))
            end do
         end do
      end if

      ! z-faces
      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               ! fluxz = rhoD * grad c
               fluxz(i,j,k) = rhoz(i,j,k)*chiz(i,j,k) * (c(i,j,k) - c(i,j,k-1)) / dx(3)
            end do
         end do
      end do

      ! alter grad(c) stencil at boundary since ghost value represents value at boundary
      if (adv_bc(3,1,5) .eq. FOEXTRAP .or. adv_bc(3,1,5) .eq. EXT_DIR) then
         k=lo(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               fluxz(i,j,k) = rhoz(i,j,k)*chiz(i,j,k) * (c(i,j,k)-c(i,j,k-1)) / (0.5d0*dx(3))
            end do
         end do
      end if
      if (adv_bc(3,2,5) .eq. FOEXTRAP .or. adv_bc(3,2,5) .eq. EXT_DIR) then
         k=hi(3)+1
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               fluxz(i,j,k) = rhoz(i,j,k)*chiz(i,j,k) * (c(i,j,k)-c(i,j,k-1)) / (0.5d0*dx(3))
            end do
         end do
      end if

      ! flux divergence
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               ! add div(rhoD * grad c) to s_update
               s_update(i,j,k) = s_update(i,j,k) + (fluxx(i+1,j,k) - fluxx(i,j,k)) / dx(1) &
                                                 + (fluxy(i,j+1,k) - fluxy(i,j,k)) / dx(2) &
                                                 + (fluxz(i,j,k+1) - fluxz(i,j,k)) / dx(3)
            end do
         end do
      end do

    end subroutine mk_diffusive_rhoc_fluxdiv_3d

  end subroutine mk_diffusive_rhoc_fluxdiv

end module mk_diffusive_fluxdiv_module
