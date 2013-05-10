module mk_diffusive_fluxdiv_module

  use ml_layout_module
  use define_bc_module
  use bc_module
  use convert_stag_module
  use probin_lowmach_module, only: visc_coef, diff_coef

  implicit none

  private

  public :: mk_diffusive_rhoc_fluxdiv, mk_diffusive_m_fluxdiv

contains

  subroutine mk_diffusive_rhoc_fluxdiv(mla,s_update,out_comp,prim,rho_fc,chi, &
                                       dx,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: s_update(:)
    integer        , intent(in   ) :: out_comp    ! which component of s_update
    type(multifab) , intent(in   ) ::   prim(:)   ! rho and c
    type(multifab) , intent(in   ) :: rho_fc(:,:) ! rho on faces
    type(multifab) , intent(in   ) ::    chi(:)   ! chi
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

    type(multifab) :: chi_fc(mla%nlevel,mla%dim)

    ng_u = s_update(1)%ng
    ng_p = prim(1)%ng
    ng_s = rho_fc(1,1)%ng
    ng_c = chi_fc(1,1)%ng

    nlevs = mla%nlevel
    dm    = mla%dim

    do n=1,nlevs
       do i=1,dm
          call multifab_build_edge(chi_fc(n,i),mla%la(n),1,0,i)
       end do
    end do

    ! average chi to faces
    if (diff_coef < 0) then
       call average_cc_to_face(nlevs,chi,chi_fc,1,dm+2,1,the_bc_level)
    else
       do n=1,nlevs
          do i=1,dm
             call setval(chi_fc(n,i),diff_coef,all=.true.)
          end do
       end do
    end if

    ! compute del dot (rhoD grad c) and add it to s_update
    do n=1,nlevs
       do i=1,nfabs(prim(n))
          up  => dataptr(s_update(n), i)
          pp  => dataptr(prim(n), i)
          spx => dataptr(rho_fc(n,1), i)
          spy => dataptr(rho_fc(n,2), i)
          cpx => dataptr(chi_fc(n,1), i)
          cpy => dataptr(chi_fc(n,2), i)
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
             spz => dataptr(rho_fc(n,3), i)
             cpz => dataptr(chi_fc(n,3), i)
             call mk_diffusive_rhoc_fluxdiv_3d(up(:,:,:,out_comp), ng_u, pp(:,:,:,2), ng_p, &
                                               spx(:,:,:,1), spy(:,:,:,1), spz(:,:,:,1), ng_s, &
                                               cpx(:,:,:,1), cpy(:,:,:,1), cpz(:,:,:,1), ng_c, &
                                               lo, hi, dx(n,:), &
                                               the_bc_level(n)%adv_bc_level_array(i,:,:,:))
          end select
       end do
    end do

    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(chi_fc(n,i))
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

  subroutine mk_diffusive_m_fluxdiv(mla,m_update,umac,eta,kappa,dx,the_bc_level)

    use stag_applyop_module, only: stag_applyop_2d, stag_applyop_3d

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: m_update(:,:)
    type(multifab) , intent(in   ) ::     umac(:,:)
    type(multifab) , intent(in   ) ::      eta(:)
    type(multifab) , intent(in   ) ::    kappa(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local variables
    integer :: nlevs,dm,i,n

    type(multifab) :: Lphi_fc(mla%nlevel,mla%dim)
    type(multifab) :: alpha_fc(mla%nlevel,mla%dim)

    type(multifab), allocatable :: eta_ed(:,:)

    logical :: nodal_temp(mla%dim)

    nlevs = mla%nlevel
    dm    = mla%dim

    if (dm .eq. 2) then
       allocate(eta_ed(nlevs,1))  ! nodal
    else if (dm .eq. 3) then
       allocate(eta_ed(nlevs,3))  ! edge-based
    end if

    do n=1,nlevs
       do i=1,dm
          call multifab_build_edge(Lphi_fc(n,i),mla%la(n),1,0,i)
          call multifab_build_edge(alpha_fc(n,i),mla%la(n),1,0,i)
          ! set alpha to zero
          call setval(alpha_fc(n,i),0.d0,all=.true.)
       end do
    end do

    ! nodal (in 2D) and edge-based (in 3D) eta
    if (dm .eq. 2) then
       do n=1,nlevs
          call multifab_build_nodal(eta_ed(n,1),mla%la(n),1,0)
       end do
    else
       do n=1,nlevs
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
       end do
    end if

    ! compute eta on nodes (2D) or edges (3D)
    if (dm .eq. 2) then
       if (visc_coef < 0) then
          call average_cc_to_node(nlevs,eta,eta_ed(:,1),1,dm+2,1,the_bc_level)
       else
          do n=1,nlevs
             call setval(eta_ed(n,1),visc_coef,all=.true.)
          end do
       end if
    else if (dm .eq. 3) then
       if (visc_coef < 0) then
          call average_cc_to_edge(nlevs,eta,eta_ed,1,dm+2,1,the_bc_level)
       else
          do n=1,nlevs
             do i=1,dm
                call setval(eta_ed(n,i),visc_coef,all=.true.)
             end do
          end do
       end if
    end if

    do n=1,nlevs

       ! compute -L(phi)
       ! we could compute +L(phi) but then we'd have to multiply beta and kappa by -1
       if (dm .eq. 2) then
          call stag_applyop_2d(mla%la(n),the_bc_level(n),umac(n,:),Lphi_fc(n,:), &
                               alpha_fc(n,:),eta(n),eta_ed(n,1),kappa(n),dx(n,:))
       else
          call stag_applyop_3d(mla%la(n),the_bc_level(n),umac(n,:),Lphi_fc(n,:), &
                               alpha_fc(n,:),eta(n),eta_ed(n,:),kappa(n),dx(n,:))
       end if

       ! subtract -L(phi) to m_update
       do i=1,dm
          call multifab_sub_sub_c(m_update(n,i),1,Lphi_fc(n,i),1,1,0)
       end do

    end do

    do n=1,nlevs
       if (dm .eq. 2) then
          call multifab_destroy(eta_ed(n,1))
       else if (dm .eq. 3) then
          call multifab_destroy(eta_ed(n,1))
          call multifab_destroy(eta_ed(n,2))
          call multifab_destroy(eta_ed(n,3))
       end if
       do i=1,dm
          call multifab_destroy(Lphi_fc(n,i))
          call multifab_destroy(alpha_fc(n,i))
       end do
    end do

  end subroutine mk_diffusive_m_fluxdiv

end module mk_diffusive_fluxdiv_module
