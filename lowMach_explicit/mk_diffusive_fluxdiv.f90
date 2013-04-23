module mk_diffusive_fluxdiv_module

  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: mk_diffusive_rhoc_fluxdiv, mk_diffusive_m_fluxdiv

contains

  subroutine mk_diffusive_rhoc_fluxdiv(mla,s_update,prim,s_face,chi_face, &
                                              umac,dx,the_bc_level)

    use bc_module
    use probin_module, only: rhobar

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: s_update(:)
    type(multifab) , intent(in   ) ::     prim(:)   ! rho and c
    type(multifab) , intent(in   ) ::   s_face(:,:) ! rho and rho1 on faces
    type(multifab) , intent(in   ) :: chi_face(:,:) ! chi on faces
    type(multifab) , intent(in   ) ::     umac(:,:) ! umac on faces
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local variables
    integer :: nlevs,dm,i,n,ng_u,ng_p,ng_s,ng_c,ng_m
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: pp(:,:,:,:)
    real(kind=dp_t), pointer :: spx(:,:,:,:)
    real(kind=dp_t), pointer :: spy(:,:,:,:)
    real(kind=dp_t), pointer :: spz(:,:,:,:)
    real(kind=dp_t), pointer :: cpx(:,:,:,:)
    real(kind=dp_t), pointer :: cpy(:,:,:,:)
    real(kind=dp_t), pointer :: cpz(:,:,:,:)
    real(kind=dp_t), pointer :: upx(:,:,:,:)
    real(kind=dp_t), pointer :: upy(:,:,:,:)
    real(kind=dp_t), pointer :: upz(:,:,:,:)

    ng_u = s_update(1)%ng
    ng_p = prim(1)%ng
    ng_s = s_face(1,1)%ng
    ng_c = chi_face(1,1)%ng
    ng_m = umac(1,1)%ng

    nlevs = mla%nlevel
    dm    = mla%dim

    ! compute del dot (rhoD grad c) and add it to s_update
    do n=1,nlevs
       do i=1,nfabs(prim(n))
          up  => dataptr(s_update(n), i)
          pp  => dataptr(prim(n), i)
          spx => dataptr(s_face(n,1), i)
          spy => dataptr(s_face(n,2), i)
          cpx => dataptr(chi_face(n,1), i)
          cpy => dataptr(chi_face(n,2), i)
          upx => dataptr(umac(n,1), i)
          upy => dataptr(umac(n,2), i)
          lo = lwb(get_box(prim(n), i))
          hi = upb(get_box(prim(n), i))
          select case (dm)
          case (2)
             call mk_diffusive_rhoc_fluxdiv_2d(up(:,:,1,2), ng_u, pp(:,:,1,2), ng_p, &
                                               spx(:,:,1,:), spy(:,:,1,:), ng_s, &
                                               cpx(:,:,1,1), cpy(:,:,1,1), ng_c, &
                                               upx(:,:,1,1), upy(:,:,1,1), ng_m, &
                                               lo, hi, dx(n,:), &
                                               the_bc_level(n)%adv_bc_level_array(i,:,:,:))
          case (3)
             spz => dataptr(s_face(n,3), i)
             cpz => dataptr(chi_face(n,3), i)
             upz => dataptr(umac(n,3), i)
             call mk_diffusive_rhoc_fluxdiv_3d(up(:,:,:,2), ng_u, pp(:,:,:,2), ng_p, &
                                               spx(:,:,:,1), spy(:,:,:,1), spz(:,:,:,1), ng_s, &
                                               cpx(:,:,:,1), cpy(:,:,:,1), cpz(:,:,:,1), ng_c, &
                                               upx(:,:,:,1), upy(:,:,:,1), upz(:,:,:,1), ng_m, &
                                               lo, hi, dx(n,:), &
                                               the_bc_level(n)%adv_bc_level_array(i,:,:,:))
          end select
       end do
    end do

  contains

    subroutine mk_diffusive_rhoc_fluxdiv_2d(s_update,ng_u,c,ng_p,sx,sy,ng_s, &
                                            chix,chiy,ng_c,umac,vmac,ng_m,lo,hi,dx,adv_bc)

      integer        , intent(in   ) :: lo(:), hi(:), ng_u, ng_p, ng_s, ng_c, ng_m
      real(kind=dp_t), intent(inout) :: s_update(lo(1)-ng_u:,lo(2)-ng_u:)
      real(kind=dp_t), intent(in   ) ::        c(lo(1)-ng_p:,lo(2)-ng_p:)
      real(kind=dp_t), intent(in   ) ::       sx(lo(1)-ng_s:,lo(2)-ng_s:,:)
      real(kind=dp_t), intent(in   ) ::       sy(lo(1)-ng_s:,lo(2)-ng_s:,:)
      real(kind=dp_t), intent(in   ) ::     chix(lo(1)-ng_c:,lo(2)-ng_c:)
      real(kind=dp_t), intent(in   ) ::     chiy(lo(1)-ng_c:,lo(2)-ng_c:)
      real(kind=dp_t), intent(inout) ::     umac(lo(1)-ng_m:,lo(2)-ng_m:)
      real(kind=dp_t), intent(inout) ::     vmac(lo(1)-ng_m:,lo(2)-ng_m:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: adv_bc(:,:,:)

      ! local
      integer :: i,j

      real(kind=dp_t) :: fluxx(lo(1):hi(1)+1,lo(2):hi(2))
      real(kind=dp_t) :: fluxy(lo(1):hi(1),lo(2):hi(2)+1)

      real(kind=dp_t) :: S_fac

      S_fac = 1.d0/rhobar(1)-1.d0/rhobar(2)

      ! x-faces
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
            ! fluxx = rhoD * grad c
            fluxx(i,j) = sx(i,j,1)*chix(i,j) * (c(i,j)-c(i-1,j)) / dx(1)
         end do
      end do

      ! alter grad(c) stencil at boundary since ghost value represents value at boundary
      if (adv_bc(1,1,4) .eq. FOEXTRAP .or. adv_bc(1,1,4) .eq. EXT_DIR) then
         i=lo(1)
         do j=lo(2),hi(2)
            fluxx(i,j) = sx(i,j,1)*chix(i,j) * (c(i,j)-c(i-1,j)) / (0.5d0*dx(1))
         end do
      end if
      if (adv_bc(1,2,4) .eq. FOEXTRAP .or. adv_bc(1,2,4) .eq. EXT_DIR) then
         i=hi(1)+1
         do j=lo(2),hi(2)
            fluxx(i,j) = sx(i,j,1)*chix(i,j) * (c(i,j)-c(i-1,j)) / (0.5d0*dx(1))
         end do
      end if

      ! update umac based on diffusive flux at boundary
      if (adv_bc(1,1,1) .eq. EXT_DIR .and. adv_bc(1,1,4) .eq. EXT_DIR) then
         umac(lo(1),lo(2):hi(2)) = umac(lo(1),lo(2):hi(2)) &
              + S_fac*fluxx(lo(1),lo(2):hi(2))
      end if
      if (adv_bc(1,2,1) .eq. EXT_DIR .and. adv_bc(1,2,4) .eq. EXT_DIR) then
         umac(hi(1)+1,lo(2):hi(2)) = umac(hi(1)+1,lo(2):hi(2)) &
              + S_fac*fluxx(hi(1)+1,lo(2):hi(2))
      end if

      ! y-faces
      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)
            ! fluxy = rhoD * grad c
            fluxy(i,j) = sy(i,j,1)*chiy(i,j) * (c(i,j)-c(i,j-1)) / dx(2)
         end do
      end do

      ! alter grad(c) stencil at boundary since ghost value represents value at boundary
      if (adv_bc(2,1,4) .eq. FOEXTRAP .or. adv_bc(2,1,4) .eq. EXT_DIR) then
         j=lo(2)
         do i=lo(1),hi(1)
            fluxy(i,j) = sy(i,j,1)*chiy(i,j) * (c(i,j)-c(i,j-1)) / (0.5d0*dx(2))
         end do
      end if
      if (adv_bc(2,2,4) .eq. FOEXTRAP .or. adv_bc(2,2,4) .eq. EXT_DIR) then
         j=hi(2)+1
         do i=lo(1),hi(1)
            fluxy(i,j) = sy(i,j,1)*chiy(i,j) * (c(i,j)-c(i,j-1)) / (0.5d0*dx(2))
         end do
      end if

      ! update vmac based on diffusive flux at boundary
      if (adv_bc(2,1,2) .eq. EXT_DIR .and. adv_bc(2,1,4) .eq. EXT_DIR) then
         vmac(lo(1):hi(1),lo(2)) = vmac(lo(1):hi(1),lo(2)) + &
              S_fac*fluxy(lo(1):hi(1),lo(2))
      end if
      if (adv_bc(2,2,2) .eq. EXT_DIR .and. adv_bc(2,2,4) .eq. EXT_DIR) then
         vmac(lo(1):hi(1),hi(2)+1) = vmac(lo(1):hi(1),hi(2)+1) + &
              S_fac*fluxy(lo(1):hi(1),hi(2)+1)
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
                                            chix,chiy,chiz,ng_c,umac,vmac,wmac,ng_m,lo,hi,dx,adv_bc)

      integer        , intent(in   ) :: lo(:), hi(:), ng_u, ng_p, ng_s, ng_c, ng_m
      real(kind=dp_t), intent(inout) :: s_update(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
      real(kind=dp_t), intent(in   ) ::        c(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
      real(kind=dp_t), intent(in   ) ::     rhox(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
      real(kind=dp_t), intent(in   ) ::     rhoy(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
      real(kind=dp_t), intent(in   ) ::     rhoz(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
      real(kind=dp_t), intent(in   ) ::     chix(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:)
      real(kind=dp_t), intent(in   ) ::     chiy(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:)
      real(kind=dp_t), intent(in   ) ::     chiz(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:)
      real(kind=dp_t), intent(inout) ::     umac(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
      real(kind=dp_t), intent(inout) ::     vmac(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
      real(kind=dp_t), intent(inout) ::     wmac(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: adv_bc(:,:,:)

      ! local
      integer :: i,j,k

      real(kind=dp_t) :: fluxx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3))
      real(kind=dp_t) :: fluxy(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3))
      real(kind=dp_t) :: fluxz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1)

      real(kind=dp_t) :: S_fac

      S_fac = 1.d0/rhobar(1)-1.d0/rhobar(2)

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

      ! update umac based on diffusive flux at boundary
      if (adv_bc(1,1,1) .eq. EXT_DIR .and. adv_bc(1,1,5) .eq. EXT_DIR) then
         umac(lo(1),lo(2):hi(2),lo(3):hi(3)) = umac(lo(1),lo(2):hi(2),lo(3):hi(3)) &
              + S_fac*fluxx(lo(1),lo(2):hi(2),lo(3):hi(3))
      end if
      if (adv_bc(1,2,1) .eq. EXT_DIR .and. adv_bc(1,2,5) .eq. EXT_DIR) then
         umac(hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = umac(hi(1)+1,lo(2):hi(2),lo(3):hi(3)) &
              + S_fac*fluxx(hi(1)+1,lo(2):hi(2),lo(3):hi(3))
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

      ! update vmac based on diffusive flux at boundary
      if (adv_bc(2,1,2) .eq. EXT_DIR .and. adv_bc(2,1,5) .eq. EXT_DIR) then
         vmac(lo(1):hi(1),lo(2),lo(3):hi(3)) = vmac(lo(1):hi(1),lo(2),lo(3):hi(3)) + &
              S_fac*fluxy(lo(1):hi(1),lo(2),lo(3):hi(3))
      end if
      if (adv_bc(2,2,2) .eq. EXT_DIR .and. adv_bc(2,2,5) .eq. EXT_DIR) then
         vmac(lo(1):hi(1),hi(2)+1,lo(3):hi(3)) = vmac(lo(1):hi(1),hi(2)+1,lo(3):hi(3)) + &
              S_fac*fluxy(lo(1):hi(1),hi(2)+1,lo(3):hi(3))
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

      ! update wmac based on diffusive flux at boundary
      if (adv_bc(3,1,3) .eq. EXT_DIR .and. adv_bc(3,1,5) .eq. EXT_DIR) then
         wmac(lo(1):hi(1),lo(2):hi(2),lo(3)) = wmac(lo(1):hi(1),lo(2):hi(2),lo(3)) + &
              S_fac*fluxz(lo(1):hi(1),lo(2):hi(2),lo(3))
      end if
      if (adv_bc(3,2,3) .eq. EXT_DIR .and. adv_bc(3,2,5) .eq. EXT_DIR) then
         wmac(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = wmac(lo(1):hi(1),lo(2):hi(2),hi(3)+1) + &
              S_fac*fluxz(lo(1):hi(1),lo(2):hi(2),hi(3)+1)
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

  subroutine mk_diffusive_m_fluxdiv(mla,m_update,umac,eta,eta_nodal,eta_edge, &
                                                  kappa,dx,the_bc_level)

    use stag_applyop_module, only: stag_applyop_2d, stag_applyop_3d

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) ::  m_update(:,:)
    type(multifab) , intent(in   ) ::      umac(:,:)
    type(multifab) , intent(in   ) ::       eta(:)
    type(multifab) , intent(in   ) :: eta_nodal(:)
    type(multifab) , intent(in   ) ::  eta_edge(:,:)
    type(multifab) , intent(in   ) ::     kappa(:)
    real(kind=dp_t), intent(in   ) ::  dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local variables
    integer :: nlevs,dm,i,n

    type(multifab) :: Lphi_fc(mla%nlevel,mla%dim)
    type(multifab) :: alpha_fc(mla%nlevel,mla%dim)

    nlevs = mla%nlevel
    dm    = mla%dim

    do n=1,nlevs
       do i=1,dm
          call multifab_build_edge(Lphi_fc(n,i),mla%la(n),1,0,i)
          call multifab_build_edge(alpha_fc(n,i),mla%la(n),1,0,i)
       end do
    end do

    ! set coefficients to explicitly compute
    ! (alpha - D) phi
    do n=1,nlevs
       do i=1,dm
          ! set alpha to zero
          call setval(alpha_fc(n,i),0.d0,all=.true.)
       end do
    end do

    do n=1,nlevs

       ! compute -L(phi)
       ! we could compute +L(phi) but then we'd have to multiply beta and kappa by -1
       if (dm .eq. 2) then
          call stag_applyop_2d(mla%la(n),the_bc_level(n),umac(n,:),Lphi_fc(n,:), &
                               alpha_fc(n,:),eta(n),eta_nodal(n),kappa(n),dx(n,:))
       else
          call stag_applyop_3d(mla%la(n),the_bc_level(n),umac(n,:),Lphi_fc(n,:), &
                               alpha_fc(n,:),eta(n),eta_edge(n,:),kappa(n),dx(n,:))
       end if

       ! subtract -L(phi) to m_update
       do i=1,dm
          call multifab_sub_sub_c(m_update(n,i),1,Lphi_fc(n,i),1,1,0)
       end do

    end do

    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(Lphi_fc(n,i))
          call multifab_destroy(alpha_fc(n,i))
       end do
    end do

  end subroutine mk_diffusive_m_fluxdiv

end module mk_diffusive_fluxdiv_module
