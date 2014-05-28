module mk_baro_fluxdiv_module

  use ml_layout_module
  use define_bc_module
  use bc_module
  use div_and_grad_module
  use probin_common_module, only: rhobar

  implicit none

  private

  public :: mk_baro_fluxdiv

contains

  subroutine mk_baro_fluxdiv(mla,s_update,out_comp,s_fc,chi_fc,gp0_fc, &
                             dx,the_bc_level,vel_bc_n)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: s_update(:)
    integer        , intent(in   ) :: out_comp    ! which component of s_update
    type(multifab) , intent(in   ) :: s_fc(:,:)   ! rho and rho*c on faces
    type(multifab) , intent(in   ) :: chi_fc(:,:) ! chi on faces
    type(multifab) , intent(in   ) :: gp0_fc(:,:) ! grad p0 on faces
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(multifab) , intent(inout) :: vel_bc_n(:,:) ! inhomogeneous normal vel value

    ! local variables
    integer :: nlevs,dm,i,n,ng_u,ng_g,ng_s,ng_c,ng_b
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: sup(:,:,:,:)
    real(kind=dp_t), pointer :: gpx(:,:,:,:)
    real(kind=dp_t), pointer :: gpy(:,:,:,:)
    real(kind=dp_t), pointer :: gpz(:,:,:,:)
    real(kind=dp_t), pointer :: spx(:,:,:,:)
    real(kind=dp_t), pointer :: spy(:,:,:,:)
    real(kind=dp_t), pointer :: spz(:,:,:,:)
    real(kind=dp_t), pointer :: cpx(:,:,:,:)
    real(kind=dp_t), pointer :: cpy(:,:,:,:)
    real(kind=dp_t), pointer :: cpz(:,:,:,:)
    real(kind=dp_t), pointer :: vpx(:,:,:,:)
    real(kind=dp_t), pointer :: vpy(:,:,:,:)
    real(kind=dp_t), pointer :: vpz(:,:,:,:)

    nlevs = mla%nlevel
    dm    = mla%dim

    ng_u = s_update(1)%ng
    ng_s = s_fc(1,1)%ng
    ng_c = chi_fc(1,1)%ng
    ng_b = vel_bc_n(1,1)%ng
    ng_g = gp0_fc(1,1)%ng

    ! compute del dot (rhoD grad c) and add it to s_update
    do n=1,nlevs
       do i=1,nfabs(s_update(n))
          sup => dataptr(s_update(n), i)
          gpx => dataptr(gp0_fc(n,1), i)
          gpy => dataptr(gp0_fc(n,2), i)
          spx => dataptr(s_fc(n,1), i)
          spy => dataptr(s_fc(n,2), i)
          cpx => dataptr(chi_fc(n,1), i)
          cpy => dataptr(chi_fc(n,2), i)
          vpx => dataptr(vel_bc_n(n,1), i)
          vpy => dataptr(vel_bc_n(n,2), i)
          lo = lwb(get_box(s_update(n), i))
          hi = upb(get_box(s_update(n), i))
          select case (dm)
          case (2)
             call mk_baro_fluxdiv_2d(sup(:,:,1,out_comp), ng_u, &
                                     gpx(:,:,1,1), gpy(:,:,1,1), ng_g, &
                                     spx(:,:,1,:), spy(:,:,1,:), ng_s, &
                                     cpx(:,:,1,1), cpy(:,:,1,1), ng_c, &
                                     vpx(:,:,1,1), vpy(:,:,1,1), ng_b, &
                                     lo, hi, dx(n,:), &
                                     the_bc_level(n)%phys_bc_level_array(i,:,:))
          case (3)
             gpz => dataptr(gp0_fc(n,3), i)
             spz => dataptr(s_fc(n,3), i)
             cpz => dataptr(chi_fc(n,3), i)
             vpz => dataptr(vel_bc_n(n,3), i)
             call mk_baro_fluxdiv_3d(sup(:,:,:,out_comp), ng_u, &
                                     gpx(:,:,:,1), gpy(:,:,:,1), gpz(:,:,:,1), ng_g, &
                                     spx(:,:,:,:), spy(:,:,:,:), spz(:,:,:,:), ng_s, &
                                     cpx(:,:,:,1), cpy(:,:,:,1), cpz(:,:,:,1), ng_c, &
                                     vpx(:,:,:,1), vpy(:,:,:,1), vpz(:,:,:,1), ng_b, &
                                     lo, hi, dx(n,:), &
                                     the_bc_level(n)%phys_bc_level_array(i,:,:))
          end select
       end do
    end do

  contains

    subroutine mk_baro_fluxdiv_2d(s_update,ng_u,gpx,gpy,ng_g,sx,sy,ng_s, &
                                            chix,chiy,ng_c,vel_bc_nx,vel_bc_ny,ng_b, &
                                            lo,hi,dx,bc)

      integer        , intent(in   ) :: lo(:), hi(:), ng_u, ng_g, ng_s, ng_c, ng_b
      real(kind=dp_t), intent(inout) ::  s_update(lo(1)-ng_u:,lo(2)-ng_u:)
      real(kind=dp_t), intent(in   ) ::       gpx(lo(1)-ng_g:,lo(2)-ng_g:)
      real(kind=dp_t), intent(in   ) ::       gpy(lo(1)-ng_g:,lo(2)-ng_g:)
      real(kind=dp_t), intent(in   ) ::        sx(lo(1)-ng_s:,lo(2)-ng_s:,:)
      real(kind=dp_t), intent(in   ) ::        sy(lo(1)-ng_s:,lo(2)-ng_s:,:)
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

      real(kind=dp_t) :: S_fac, c_fc, kp

      S_fac = 1.d0/rhobar(1)-1.d0/rhobar(2)

      ! x-faces
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
            ! fluxx = rho*chi*kp*gpx
            c_fc = sx(i,j,2)/sx(i,j,1)
            kp = S_fac*c_fc*(1.d0-c_fc)
            fluxx(i,j) = sx(i,j,1)*chix(i,j)*kp*gpx(i,j)
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
            ! fluxy = rho*chi*kp*gpy
            c_fc = sy(i,j,2)/sy(i,j,1)
            kp = S_fac*c_fc*(1.d0-c_fc)
            fluxy(i,j) = sy(i,j,1)*chiy(i,j)*kp*gpy(i,j)
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

    end subroutine mk_baro_fluxdiv_2d

    subroutine mk_baro_fluxdiv_3d(s_update,ng_u,gpx,gpy,gpz,ng_g,sx,sy,sz,ng_s, &
                                  chix,chiy,chiz,ng_c,vel_bc_nx,vel_bc_ny,vel_bc_nz,ng_b, &
                                  lo,hi,dx,bc)

      integer        , intent(in   ) :: lo(:), hi(:), ng_u, ng_g, ng_s, ng_c, ng_b
      real(kind=dp_t), intent(inout) ::  s_update(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
      real(kind=dp_t), intent(in   ) ::       gpx(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:)
      real(kind=dp_t), intent(in   ) ::       gpy(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:)
      real(kind=dp_t), intent(in   ) ::       gpz(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:)
      real(kind=dp_t), intent(in   ) ::        sx(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
      real(kind=dp_t), intent(in   ) ::        sy(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
      real(kind=dp_t), intent(in   ) ::        sz(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
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

      real(kind=dp_t) :: S_fac, c_fc, kp

      S_fac = 1.d0/rhobar(1)-1.d0/rhobar(2)

      ! x-faces
      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)+1
         ! fluxx = rho*chi*kp*gpx
         c_fc = sx(i,j,k,2)/sx(i,j,k,1)
         kp = S_fac*c_fc*(1.d0-c_fc)
         fluxx(i,j,k) = sx(i,j,k,1)*chix(i,j,k)*kp*gpx(i,j,k)
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
         ! fluxy = rho*chi*kp*gpy
         c_fc = sy(i,j,k,2)/sy(i,j,k,1)
         kp = S_fac*c_fc*(1.d0-c_fc)
         fluxy(i,j,k) = sy(i,j,k,1)*chiy(i,j,k)*kp*gpy(i,j,k)
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
         ! fluxz = rho*chi*kp*gpz
         c_fc = sz(i,j,k,2)/sz(i,j,k,1)
         kp = S_fac*c_fc*(1.d0-c_fc)
         fluxz(i,j,k) = sz(i,j,k,1)*chiz(i,j,k)*kp*gpz(i,j,k)
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

    end subroutine mk_baro_fluxdiv_3d

  end subroutine mk_baro_fluxdiv

end module mk_baro_fluxdiv_module
