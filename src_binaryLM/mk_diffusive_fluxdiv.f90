module mk_diffusive_fluxdiv_module

  use ml_layout_module
  use define_bc_module
  use bc_module
  use div_and_grad_module
  use probin_binarylm_module, only: rhobar

  implicit none

  private

  public :: mk_diffusive_rhoc_fluxdiv

contains

  subroutine mk_diffusive_rhoc_fluxdiv(mla,s_update,out_comp,prim,rho_fc,chi_fc, &
                                       dx,the_bc_level,vel_bc_n)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: s_update(:)
    integer        , intent(in   ) :: out_comp    ! which component of s_update
    type(multifab) , intent(in   ) ::   prim(:)   ! rho and c
    type(multifab) , intent(in   ) :: rho_fc(:,:) ! rho on faces
    type(multifab) , intent(in   ) :: chi_fc(:,:) ! chi on faces
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(multifab) , intent(inout) :: vel_bc_n(:,:) ! inhomogeneous normal vel value

    ! local variables
    integer :: nlevs,dm,i,n,ng_u,ng_g,ng_s,ng_c,ng_b
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: gcx(:,:,:,:)
    real(kind=dp_t), pointer :: gcy(:,:,:,:)
    real(kind=dp_t), pointer :: gcz(:,:,:,:)
    real(kind=dp_t), pointer :: spx(:,:,:,:)
    real(kind=dp_t), pointer :: spy(:,:,:,:)
    real(kind=dp_t), pointer :: spz(:,:,:,:)
    real(kind=dp_t), pointer :: cpx(:,:,:,:)
    real(kind=dp_t), pointer :: cpy(:,:,:,:)
    real(kind=dp_t), pointer :: cpz(:,:,:,:)
    real(kind=dp_t), pointer :: vpx(:,:,:,:)
    real(kind=dp_t), pointer :: vpy(:,:,:,:)
    real(kind=dp_t), pointer :: vpz(:,:,:,:)

    type(multifab) :: gradc(mla%nlevel,mla%dim)

    nlevs = mla%nlevel
    dm    = mla%dim

    do n=1,nlevs
       do i=1,dm
          call multifab_build_edge(gradc(n,i),mla%la(n),1,0,i)
       end do
    end do

    ng_u = s_update(1)%ng
    ng_s = rho_fc(1,1)%ng
    ng_c = chi_fc(1,1)%ng
    ng_b = vel_bc_n(1,1)%ng
    ng_g = gradc(1,1)%ng

    call compute_grad(mla,prim,gradc,dx,2,scal_bc_comp+1,1,1,the_bc_level)

    ! compute del dot (rhoD grad c) and add it to s_update
    do n=1,nlevs
       do i=1,nfabs(prim(n))
          up  => dataptr(s_update(n), i)
          gcx => dataptr(gradc(n,1), i)
          gcy => dataptr(gradc(n,2), i)
          spx => dataptr(rho_fc(n,1), i)
          spy => dataptr(rho_fc(n,2), i)
          cpx => dataptr(chi_fc(n,1), i)
          cpy => dataptr(chi_fc(n,2), i)
          vpx => dataptr(vel_bc_n(n,1), i)
          vpy => dataptr(vel_bc_n(n,2), i)
          lo = lwb(get_box(prim(n), i))
          hi = upb(get_box(prim(n), i))
          select case (dm)
          case (2)
             call mk_diffusive_rhoc_fluxdiv_2d(up(:,:,1,out_comp), ng_u, &
                                               gcx(:,:,1,1), gcy(:,:,1,1), ng_g, &
                                               spx(:,:,1,1), spy(:,:,1,1), ng_s, &
                                               cpx(:,:,1,1), cpy(:,:,1,1), ng_c, &
                                               vpx(:,:,1,1), vpy(:,:,1,1), ng_b, &
                                               lo, hi, dx(n,:), &
                                               the_bc_level(n)%phys_bc_level_array(i,:,:))
          case (3)
             gcz => dataptr(gradc(n,3), i)
             spz => dataptr(rho_fc(n,3), i)
             cpz => dataptr(chi_fc(n,3), i)
             vpz => dataptr(vel_bc_n(n,3), i)
             call mk_diffusive_rhoc_fluxdiv_3d(up(:,:,:,out_comp), ng_u, &
                                               gcx(:,:,:,1), gcy(:,:,:,1), gcz(:,:,:,1), ng_g, &
                                               spx(:,:,:,1), spy(:,:,:,1), spz(:,:,:,1), ng_s, &
                                               cpx(:,:,:,1), cpy(:,:,:,1), cpz(:,:,:,1), ng_c, &
                                               vpx(:,:,:,1), vpy(:,:,:,1), vpz(:,:,:,1), ng_b, &
                                               lo, hi, dx(n,:), &
                                               the_bc_level(n)%phys_bc_level_array(i,:,:))
          end select
       end do
    end do

    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(gradc(n,i))
       end do
    end do

  contains

    subroutine mk_diffusive_rhoc_fluxdiv_2d(s_update,ng_u,gradcx,gradcy,ng_g,rhox,rhoy,ng_s, &
                                            chix,chiy,ng_c,vel_bc_nx,vel_bc_ny,ng_b, &
                                            lo,hi,dx,bc)

      integer        , intent(in   ) :: lo(:), hi(:), ng_u, ng_g, ng_s, ng_c, ng_b
      real(kind=dp_t), intent(inout) ::  s_update(lo(1)-ng_u:,lo(2)-ng_u:)
      real(kind=dp_t), intent(in   ) ::    gradcx(lo(1)-ng_g:,lo(2)-ng_g:)
      real(kind=dp_t), intent(in   ) ::    gradcy(lo(1)-ng_g:,lo(2)-ng_g:)
      real(kind=dp_t), intent(in   ) ::      rhox(lo(1)-ng_s:,lo(2)-ng_s:)
      real(kind=dp_t), intent(in   ) ::      rhoy(lo(1)-ng_s:,lo(2)-ng_s:)
      real(kind=dp_t), intent(in   ) ::      chix(lo(1)-ng_c:,lo(2)-ng_c:)
      real(kind=dp_t), intent(in   ) ::      chiy(lo(1)-ng_c:,lo(2)-ng_c:)
      real(kind=dp_t), intent(inout) :: vel_bc_nx(lo(1)-ng_b:,lo(2)-ng_b:)
      real(kind=dp_t), intent(inout) :: vel_bc_ny(lo(1)-ng_b:,lo(2)-ng_b:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: bc(:,:)

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
            fluxx(i,j) = rhox(i,j)*chix(i,j) * gradcx(i,j)
         end do
      end do

      ! update umac bc based on diffusive flux at boundary
      if (bc(1,1) .eq. NO_SLIP_RESERVOIR .or. bc(1,1) .eq. SLIP_RESERVOIR) then
         vel_bc_nx(lo(1),lo(2):hi(2)) = vel_bc_nx(lo(1),lo(2):hi(2)) &
              + S_fac*fluxx(lo(1),lo(2):hi(2))
      end if
      if (bc(1,2) .eq. NO_SLIP_RESERVOIR .or. bc(1,2) .eq. SLIP_RESERVOIR) then
         vel_bc_nx(hi(1)+1,lo(2):hi(2)) = vel_bc_nx(hi(1)+1,lo(2):hi(2)) &
              + S_fac*fluxx(hi(1)+1,lo(2):hi(2))
      end if

      ! y-faces
      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)
            ! fluxy = rhoD * grad c
            fluxy(i,j) = rhoy(i,j)*chiy(i,j) * gradcy(i,j)
         end do
      end do

      ! update vmac bc based on diffusive flux at boundary
      if (bc(2,1) .eq. NO_SLIP_RESERVOIR .or. bc(2,1) .eq. SLIP_RESERVOIR) then
         vel_bc_ny(lo(1):hi(1),lo(2)) = vel_bc_ny(lo(1):hi(1),lo(2)) + &
              S_fac*fluxy(lo(1):hi(1),lo(2))
      end if
      if (bc(2,2) .eq. NO_SLIP_RESERVOIR .or. bc(2,2) .eq. SLIP_RESERVOIR) then
         vel_bc_ny(lo(1):hi(1),hi(2)+1) = vel_bc_ny(lo(1):hi(1),hi(2)+1) + &
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

    subroutine mk_diffusive_rhoc_fluxdiv_3d(s_update,ng_u,gradcx,gradcy,gradcz,ng_g, &
                                            rhox,rhoy,rhoz,ng_s, &
                                            chix,chiy,chiz,ng_c,vel_bc_nx, &
                                            vel_bc_ny,vel_bc_nz,ng_b,lo,hi,dx,bc)

      integer        , intent(in   ) :: lo(:), hi(:), ng_u, ng_g, ng_s, ng_c, ng_b
      real(kind=dp_t), intent(inout) ::  s_update(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
      real(kind=dp_t), intent(in   ) ::    gradcx(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:)
      real(kind=dp_t), intent(in   ) ::    gradcy(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:)
      real(kind=dp_t), intent(in   ) ::    gradcz(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:)
      real(kind=dp_t), intent(in   ) ::      rhox(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
      real(kind=dp_t), intent(in   ) ::      rhoy(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
      real(kind=dp_t), intent(in   ) ::      rhoz(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
      real(kind=dp_t), intent(in   ) ::      chix(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:)
      real(kind=dp_t), intent(in   ) ::      chiy(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:)
      real(kind=dp_t), intent(in   ) ::      chiz(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:)
      real(kind=dp_t), intent(inout) :: vel_bc_nx(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
      real(kind=dp_t), intent(inout) :: vel_bc_ny(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
      real(kind=dp_t), intent(inout) :: vel_bc_nz(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: bc(:,:)

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
               fluxx(i,j,k) = rhox(i,j,k)*chix(i,j,k) * gradcx(i,j,k)
            end do
         end do
      end do

      ! update umac bc based on diffusive flux at boundary
      if (bc(1,1) .eq. NO_SLIP_RESERVOIR .or. bc(1,1) .eq. SLIP_RESERVOIR) then
         vel_bc_nx(lo(1),lo(2):hi(2),lo(3):hi(3)) = vel_bc_nx(lo(1),lo(2):hi(2),lo(3):hi(3)) &
              + S_fac*fluxx(lo(1),lo(2):hi(2),lo(3):hi(3))
      end if
      if (bc(1,2) .eq. NO_SLIP_RESERVOIR .or. bc(1,2) .eq. SLIP_RESERVOIR) then
         vel_bc_nx(hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = vel_bc_nx(hi(1)+1,lo(2):hi(2),lo(3):hi(3)) &
              + S_fac*fluxx(hi(1)+1,lo(2):hi(2),lo(3):hi(3))
      end if

      ! y-faces
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)
               ! fluxy = rhoD * grad c
               fluxy(i,j,k) = rhoy(i,j,k)*chiy(i,j,k) * gradcy(i,j,k)
            end do
         end do
      end do

      ! update vmac bc based on diffusive flux at boundary
      if (bc(2,1) .eq. NO_SLIP_RESERVOIR .or. bc(2,1) .eq. SLIP_RESERVOIR) then
         vel_bc_ny(lo(1):hi(1),lo(2),lo(3):hi(3)) = vel_bc_ny(lo(1):hi(1),lo(2),lo(3):hi(3)) + &
              S_fac*fluxy(lo(1):hi(1),lo(2),lo(3):hi(3))
      end if
      if (bc(2,2) .eq. NO_SLIP_RESERVOIR .or. bc(2,2) .eq. SLIP_RESERVOIR) then
         vel_bc_ny(lo(1):hi(1),hi(2)+1,lo(3):hi(3)) = vel_bc_ny(lo(1):hi(1),hi(2)+1,lo(3):hi(3)) + &
              S_fac*fluxy(lo(1):hi(1),hi(2)+1,lo(3):hi(3))
      end if

      ! z-faces
      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               ! fluxz = rhoD * grad c
               fluxz(i,j,k) = rhoz(i,j,k)*chiz(i,j,k) * gradcz(i,j,k)
            end do
         end do
      end do

      ! update wmac bc based on diffusive flux at boundary
      if (bc(3,1) .eq. NO_SLIP_RESERVOIR .or. bc(3,1) .eq. SLIP_RESERVOIR) then
         vel_bc_nz(lo(1):hi(1),lo(2):hi(2),lo(3)) = vel_bc_nz(lo(1):hi(1),lo(2):hi(2),lo(3)) + &
              S_fac*fluxz(lo(1):hi(1),lo(2):hi(2),lo(3))
      end if
      if (bc(3,2) .eq. NO_SLIP_RESERVOIR .or. bc(3,2) .eq. SLIP_RESERVOIR) then
         vel_bc_nz(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = vel_bc_nz(lo(1):hi(1),lo(2):hi(2),hi(3)+1) + &
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

end module mk_diffusive_fluxdiv_module
