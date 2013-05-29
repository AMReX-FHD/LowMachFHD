module set_inhomogeneous_vel_bcs_module

  use ml_layout_module
  use multifab_module
  use define_bc_module
  use bc_module
  use inhomogeneous_bc_val_module
  use probin_common_module, only: prob_lo, prob_hi

  implicit none

  private

  public :: set_inhomogeneous_vel_bcs, modify_traction_bcs

contains

  subroutine set_inhomogeneous_vel_bcs(mla,vel_bc_n,vel_bc_t,eta_ed,dx,the_bc_level)

    ! vel_bc_n(nlevs,dm) are the normal velocities
    ! in 2D, vel_bc_t(nlevs,2) respresents
    !   1. y-velocity bc on x-faces (nodal)
    !   2. x-velocity bc on y-faces (nodal)
    ! in 3D, vel_bc_t(nlevs,6) represents
    !   1. y-velocity bc on x-faces (nodal in y and x)
    !   2. z-velocity bc on x-faces (nodal in z and x)
    !   3. x-velocity bc on y-faces (nodal in x and y)
    !   4. z-velocity bc on y-faces (nodal in z and y)
    !   5. x-velocity bc on z-faces (nodal in x and z)
    !   6. y-velocity bc on z-faces (nodal in y and z)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: vel_bc_n(:,:)
    type(multifab) , intent(inout) :: vel_bc_t(:,:)
    type(multifab) , intent(in   ) :: eta_ed(:,:)   ! nodal in 2D, edge-based in 3D
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local
    integer :: n,nlevs,i,dm,ng_n,ng_t,ng_e
    integer :: lo(mla%dim),hi(mla%dim)
    real(kind=dp_t), pointer :: nxp(:,:,:,:)
    real(kind=dp_t), pointer :: nyp(:,:,:,:)
    real(kind=dp_t), pointer :: nzp(:,:,:,:)
    real(kind=dp_t), pointer :: t1p(:,:,:,:)
    real(kind=dp_t), pointer :: t2p(:,:,:,:)
    real(kind=dp_t), pointer :: t3p(:,:,:,:)
    real(kind=dp_t), pointer :: t4p(:,:,:,:)
    real(kind=dp_t), pointer :: t5p(:,:,:,:)
    real(kind=dp_t), pointer :: t6p(:,:,:,:)
    real(kind=dp_t), pointer :: enp(:,:,:,:)
    real(kind=dp_t), pointer :: ep1(:,:,:,:)
    real(kind=dp_t), pointer :: ep2(:,:,:,:)
    real(kind=dp_t), pointer :: ep3(:,:,:,:)

    nlevs = mla%nlevel
    dm = mla%dim

    ng_n = vel_bc_n(1,1)%ng
    ng_t = vel_bc_t(1,1)%ng
    ng_e = eta_ed(1,1)%ng

    do n=1,nlevs
       do i=1,nfabs(vel_bc_n(n,1))
          nxp => dataptr(vel_bc_n(n,1),i)
          nyp => dataptr(vel_bc_n(n,2),i)
          t1p => dataptr(vel_bc_t(n,1),i)
          t2p => dataptr(vel_bc_t(n,2),i)
          enp => dataptr(eta_ed(n,1),i)
          lo = lwb(get_box(vel_bc_n(n,1),i))
          hi = upb(get_box(vel_bc_n(n,1),i))
          select case (dm)
          case (2)
             call set_inhomogeneous_vel_bcs_2d(nxp(:,:,1,1),nyp(:,:,1,1),ng_n, &
                                               t1p(:,:,1,1),t2p(:,:,1,1),ng_t, &
                                               enp(:,:,1,1),ng_e, &
                                               lo,hi,dx(n,:), &
                                               the_bc_level(n)%adv_bc_level_array(i,:,:,:))
          case (3)
             nzp => dataptr(vel_bc_n(n,3),i)
             t3p => dataptr(vel_bc_t(n,3),i)
             t4p => dataptr(vel_bc_t(n,4),i)
             t5p => dataptr(vel_bc_t(n,5),i)
             t6p => dataptr(vel_bc_t(n,6),i)
             ep1 => dataptr(eta_ed(n,1),i)
             ep2 => dataptr(eta_ed(n,2),i)
             ep3 => dataptr(eta_ed(n,3),i)
             call set_inhomogeneous_vel_bcs_3d(nxp(:,:,:,1),nyp(:,:,:,1),nzp(:,:,:,1),ng_n, &
                                               t1p(:,:,:,1),t2p(:,:,:,1),t3p(:,:,:,1), &
                                               t4p(:,:,:,1),t5p(:,:,:,1),t6p(:,:,:,1),ng_t, &
                                               ep1(:,:,:,1),ep2(:,:,:,1),ep3(:,:,:,1),ng_e, &
                                               lo,hi,dx(n,:), &
                                               the_bc_level(n)%adv_bc_level_array(i,:,:,:))

          end select
       end do
    end do
    
  end subroutine set_inhomogeneous_vel_bcs

  subroutine set_inhomogeneous_vel_bcs_2d(vel_bc_nx,vel_bc_ny,ng_n, &
                                          vel_bc_tyx,vel_bc_txy,ng_t, &
                                          eta_nd,ng_e,lo,hi,dx,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng_n,ng_t,ng_e
    real(kind=dp_t), intent(inout) ::  vel_bc_nx(lo(1)-ng_n:,lo(2)-ng_n:)
    real(kind=dp_t), intent(inout) ::  vel_bc_ny(lo(1)-ng_n:,lo(2)-ng_n:)
    real(kind=dp_t), intent(inout) :: vel_bc_tyx(lo(1)-ng_t:,lo(2)-ng_t:)
    real(kind=dp_t), intent(inout) :: vel_bc_txy(lo(1)-ng_t:,lo(2)-ng_t:)
    real(kind=dp_t), intent(inout) ::     eta_nd(lo(1)-ng_e:,lo(2)-ng_e:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ) :: bc(:,:,:)

    ! local
    integer :: i,j
    real(kind=dp_t) :: x,y

    !!!!!!!!!!!!!!!!!!!!
    ! normal velocities
    !!!!!!!!!!!!!!!!!!!!

    ! xvel, lo x-faces
    if (bc(1,1,1) .eq. DIR_VEL) then
       x = prob_lo(1)
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          vel_bc_nx(lo(1),j) = inhomogeneous_bc_val_2d(1,x,y)
       end do
    end if

    ! xvel, hi x-faces
    if (bc(1,2,1) .eq. DIR_VEL) then
       x = prob_hi(1)
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          vel_bc_nx(hi(1)+1,j) = inhomogeneous_bc_val_2d(1,x,y)
       end do
    end if

    ! yvel, lo y-faces
    if (bc(2,1,2) .eq. DIR_VEL) then
       y = prob_lo(2)
       do i=lo(1),hi(1)
          x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
          vel_bc_ny(i,lo(2)) = inhomogeneous_bc_val_2d(2,x,y)
       end do
    end if

    ! yvel, hi y-faces
    if (bc(2,2,2) .eq. DIR_VEL) then
       y = prob_hi(2)
       do i=lo(1),hi(1)
          x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
          vel_bc_ny(i,hi(2)+1) = inhomogeneous_bc_val_2d(2,x,y)
       end do
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!
    ! transverse velocities
    !!!!!!!!!!!!!!!!!!!!!!!!

    ! yvel, lo x-faces
    if (bc(1,1,2) .eq. DIR_VEL .or. bc(1,1,2) .eq. DIR_TRACT) then
       x = prob_lo(1)
       do j=lo(2),hi(2)+1
          y = prob_lo(2) + dble(j)*dx(2)
          vel_bc_tyx(lo(1),j) = inhomogeneous_bc_val_2d(2,x,y) / eta_nd(lo(1),j)
       end do
    end if

    ! yvel, hi x-faces
    if (bc(1,2,2) .eq. DIR_VEL .or. bc(1,2,2) .eq. DIR_TRACT) then
       x = prob_hi(1)
       do j=lo(2),hi(2)+1
          y = prob_lo(2) + dble(j)*dx(2)
          vel_bc_tyx(hi(1)+1,j) = inhomogeneous_bc_val_2d(2,x,y) / eta_nd(hi(1)+1,j)
       end do
    end if

    ! xvel, lo y-faces
    if (bc(2,1,1) .eq. DIR_VEL .or. bc(2,1,1) .eq. DIR_TRACT) then
       y = prob_lo(2)
       do i=lo(1),hi(1)+1
          x = prob_lo(1) + dble(i)*dx(1)
          vel_bc_txy(i,lo(2)) = inhomogeneous_bc_val_2d(1,x,y) / eta_nd(i,lo(2))
       end do
    end if

    ! xvel, hi y-faces
    if (bc(2,2,1) .eq. DIR_VEL .or. bc(2,2,1) .eq. DIR_TRACT) then
       y = prob_hi(2)
       do i=lo(1),hi(1)+1
          x = prob_lo(1) + dble(i)*dx(1)
          vel_bc_txy(i,hi(2)+1) = inhomogeneous_bc_val_2d(1,x,y) / eta_nd(i,hi(2)+1)
       end do
    end if

  end subroutine set_inhomogeneous_vel_bcs_2d

  subroutine set_inhomogeneous_vel_bcs_3d(vel_bc_nx,vel_bc_ny,vel_bc_nz,ng_n, &
                                          vel_bc_tyx,vel_bc_tzx,vel_bc_txy,vel_bc_tzy, &
                                          vel_bc_txz,vel_bc_tyz,ng_t, &
                                          eta_xy,eta_xz,eta_yz,ng_e,lo,hi,dx,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng_n,ng_t,ng_e
    real(kind=dp_t), intent(inout) ::  vel_bc_nx(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:)
    real(kind=dp_t), intent(inout) ::  vel_bc_ny(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:)
    real(kind=dp_t), intent(inout) ::  vel_bc_nz(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:)
    real(kind=dp_t), intent(inout) :: vel_bc_tyx(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
    real(kind=dp_t), intent(inout) :: vel_bc_tzx(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
    real(kind=dp_t), intent(inout) :: vel_bc_txy(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
    real(kind=dp_t), intent(inout) :: vel_bc_tzy(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
    real(kind=dp_t), intent(inout) :: vel_bc_txz(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
    real(kind=dp_t), intent(inout) :: vel_bc_tyz(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
    real(kind=dp_t), intent(in   ) ::     eta_xy(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(in   ) ::     eta_xz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(in   ) ::     eta_yz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ) :: bc(:,:,:)

    ! local
    integer :: i,j,k
    real(kind=dp_t) :: x,y,z

    !!!!!!!!!!!!!!!!!!!!
    ! normal velocities
    !!!!!!!!!!!!!!!!!!!!

    ! xvel, lo x-faces
    if (bc(1,1,1) .eq. DIR_VEL) then
       x = prob_lo(1)
       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             vel_bc_nx(lo(1),j,k) = inhomogeneous_bc_val_3d(1,x,y,z)
          end do
       end do
    end if

    ! xvel, hi x-faces
    if (bc(1,2,1) .eq. DIR_VEL) then
       x = prob_hi(1)
       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             vel_bc_nx(hi(1)+1,j,k) = inhomogeneous_bc_val_3d(1,x,y,z)
          end do
       end do
    end if
    
    ! yvel, lo y-faces
    if (bc(2,1,2) .eq. DIR_VEL) then
       y = prob_lo(2)
       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             vel_bc_ny(i,lo(2),k) = inhomogeneous_bc_val_3d(2,x,y,z)
          end do
       end do
    end if

    ! yvel, hi y-faces
    if (bc(2,2,2) .eq. DIR_VEL) then
       y = prob_hi(2)
       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             vel_bc_ny(i,hi(2)+1,k) = inhomogeneous_bc_val_3d(2,x,y,z)
          end do
       end do
    end if
    
    ! zvel, lo z-faces
    if (bc(3,1,3) .eq. DIR_VEL) then
       z = prob_lo(3)
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             vel_bc_nz(i,j,lo(3)) = inhomogeneous_bc_val_3d(3,x,y,z)
          end do
       end do
    end if

    ! zvel, hi z-faces
    if (bc(3,2,3) .eq. DIR_VEL) then
       z = prob_hi(3)
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             vel_bc_nz(i,j,hi(3)+1) = inhomogeneous_bc_val_3d(3,x,y,z)
          end do
       end do
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!
    ! transverse velocities
    !!!!!!!!!!!!!!!!!!!!!!!!
    
    ! yvel, lo x-faces
    if (bc(1,1,2) .eq. DIR_VEL .or. bc(1,1,2) .eq. DIR_TRACT) then
       x = prob_lo(1)
       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do j=lo(2),hi(2)+1
             y = prob_lo(2) + dble(j)*dx(2)
             vel_bc_tyx(lo(1),j,k) = inhomogeneous_bc_val_3d(2,x,y,z) / eta_xy(lo(1),j,k)
          end do
       end do

    end if

    ! yvel, hi x-faces
    if (bc(1,2,2) .eq. DIR_VEL .or. bc(1,2,2) .eq. DIR_TRACT) then
       x = prob_hi(1)
       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do j=lo(2),hi(2)+1
             y = prob_lo(2) + dble(j)*dx(2)
             vel_bc_tyx(hi(1)+1,j,k) = inhomogeneous_bc_val_3d(2,x,y,z) / eta_xy(hi(1)+1,j,k)
          end do
       end do

    end if

    ! zvel, lo x-faces
    if (bc(1,1,3) .eq. DIR_VEL .or. bc(1,1,3) .eq. DIR_TRACT) then
       x = prob_lo(1)
       do k=lo(3),hi(3)+1
          z = prob_lo(3) + dble(k)*dx(3)
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             vel_bc_tzx(lo(1),j,k) = inhomogeneous_bc_val_3d(3,x,y,z) / eta_xz(lo(1),j,k)
          end do
       end do

    end if

    ! zvel, hi x-faces
    if (bc(1,2,3) .eq. DIR_VEL .or. bc(1,2,3) .eq. DIR_TRACT) then
       x = prob_hi(1)
       do k=lo(3),hi(3)+1
          z = prob_lo(3) + dble(k)*dx(3)
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             vel_bc_tzx(hi(1)+1,j,k) = inhomogeneous_bc_val_3d(3,x,y,z) / eta_xz(hi(1)+1,j,k)
          end do
       end do

    end if

    ! xvel, lo y-faces
    if (bc(2,1,1) .eq. DIR_VEL .or. bc(2,1,1) .eq. DIR_TRACT) then
       y = prob_lo(2)
       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do i=lo(1),hi(1)+1
             x = prob_lo(1) + dble(i)*dx(1)
             vel_bc_txy(i,lo(2),k) = inhomogeneous_bc_val_3d(1,x,y,z) / eta_xy(i,lo(2),k)
          end do
       end do

    end if

    ! xvel, hi y-faces
    if (bc(2,2,1) .eq. DIR_VEL .or. bc(2,2,1) .eq. DIR_TRACT) then
       y = prob_hi(2)
       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do i=lo(1),hi(1)+1
             x = prob_lo(1) + dble(i)*dx(1)
             vel_bc_txy(i,hi(2)+1,k) = inhomogeneous_bc_val_3d(1,x,y,z) / eta_xy(i,hi(2)+1,k)
          end do
       end do

    end if

    ! zvel, lo y-faces
    if (bc(2,1,3) .eq. DIR_VEL .or. bc(2,1,3) .eq. DIR_TRACT) then
       y = prob_lo(2)
       do k=lo(3),hi(3)+1
          z = prob_lo(3) + dble(k)*dx(3)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             vel_bc_tzy(i,lo(2),k) = inhomogeneous_bc_val_3d(3,x,y,z) / eta_yz(i,lo(2),k)
          end do
       end do

    end if

    ! zvel, hi y-faces
    if (bc(2,2,3) .eq. DIR_VEL .or. bc(2,2,3) .eq. DIR_TRACT) then
       y = prob_hi(2)
       do k=lo(3),hi(3)+1
          z = prob_lo(3) + dble(k)*dx(3)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             vel_bc_tzy(i,hi(2)+1,k) = inhomogeneous_bc_val_3d(3,x,y,z) / eta_yz(i,hi(2)+1,k)
          end do
       end do

    end if

    ! xvel, lo z-faces
    if (bc(3,1,1) .eq. DIR_VEL .or. bc(3,1,1) .eq. DIR_TRACT) then
       z = prob_lo(3)
       do i=lo(1),hi(1)+1
          x = prob_lo(1) + dble(i)*dx(1)
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             vel_bc_txz(i,j,lo(3)) = inhomogeneous_bc_val_3d(1,x,y,z) / eta_xz(i,j,lo(3))
          end do
       end do

    end if

    ! xvel, hi z-faces
    if (bc(3,2,1) .eq. DIR_VEL .or. bc(3,2,1) .eq. DIR_TRACT) then
       z = prob_hi(3)
       do i=lo(1),hi(1)+1
          x = prob_lo(1) + dble(i)*dx(1)
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             vel_bc_txz(i,j,hi(3)+1) = inhomogeneous_bc_val_3d(1,x,y,z) / eta_xz(i,j,hi(3)+1)
          end do
       end do

    end if

    ! yvel, lo z-faces
    if (bc(3,1,2) .eq. DIR_VEL .or. bc(3,1,2) .eq. DIR_TRACT) then
       z = prob_lo(3)
       do j=lo(2),hi(2)+1
          y = prob_lo(2) + dble(j)*dx(2)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             vel_bc_tyz(i,j,lo(3)) = inhomogeneous_bc_val_3d(2,x,y,z) / eta_yz(i,j,lo(3))
          end do
       end do

    end if

    ! yvel, hi z-faces
    if (bc(3,2,2) .eq. DIR_VEL .or. bc(3,2,2) .eq. DIR_TRACT) then
       z = prob_hi(3)
       do j=lo(2),hi(2)+1
          y = prob_lo(2) + dble(j)*dx(2)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             vel_bc_tyz(i,j,hi(3)+1) = inhomogeneous_bc_val_3d(2,x,y,z) / eta_yz(i,j,hi(3)+1)
          end do
       end do

    end if

  end subroutine set_inhomogeneous_vel_bcs_3d

  subroutine modify_traction_bcs(mla,vel_bc_n,vel_bc_t,dx,the_bc_level)

    ! This subroutine modifies the traction bc, f / eta, by subtracting
    ! off the dv_n/dt part for v_t

    ! vel_bc_n(nlevs,dm) are the normal velocities
    ! in 2D, vel_bc_t(nlevs,2) respresents
    !   1. y-velocity bc on x-faces (nodal)
    !   2. x-velocity bc on y-faces (nodal)
    ! in 3D, vel_bc_t(nlevs,6) represents
    !   1. y-velocity bc on x-faces (nodal in y and x)
    !   2. z-velocity bc on x-faces (nodal in z and x)
    !   3. x-velocity bc on y-faces (nodal in x and y)
    !   4. z-velocity bc on y-faces (nodal in z and y)
    !   5. x-velocity bc on z-faces (nodal in x and z)
    !   6. y-velocity bc on z-faces (nodal in y and z)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: vel_bc_n(:,:)
    type(multifab) , intent(inout) :: vel_bc_t(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local
    integer :: n,nlevs,i,dm,ng_n,ng_t
    integer :: lo(mla%dim),hi(mla%dim)
    real(kind=dp_t), pointer :: nxp(:,:,:,:)
    real(kind=dp_t), pointer :: nyp(:,:,:,:)
    real(kind=dp_t), pointer :: nzp(:,:,:,:)
    real(kind=dp_t), pointer :: t1p(:,:,:,:)
    real(kind=dp_t), pointer :: t2p(:,:,:,:)
    real(kind=dp_t), pointer :: t3p(:,:,:,:)
    real(kind=dp_t), pointer :: t4p(:,:,:,:)
    real(kind=dp_t), pointer :: t5p(:,:,:,:)
    real(kind=dp_t), pointer :: t6p(:,:,:,:)

    nlevs = mla%nlevel
    dm = mla%dim

    ng_n = vel_bc_n(1,1)%ng
    ng_t = vel_bc_t(1,1)%ng

    do n=1,nlevs
       do i=1,nfabs(vel_bc_n(n,1))
          nxp => dataptr(vel_bc_n(n,1),i)
          nyp => dataptr(vel_bc_n(n,2),i)
          t1p => dataptr(vel_bc_t(n,1),i)
          t2p => dataptr(vel_bc_t(n,2),i)
          lo = lwb(get_box(vel_bc_n(n,1),i))
          hi = upb(get_box(vel_bc_n(n,1),i))
          select case (dm)
          case (2)
             call modify_traction_bcs_2d(nxp(:,:,1,1),nyp(:,:,1,1),ng_n, &
                                               t1p(:,:,1,1),t2p(:,:,1,1),ng_t, &
                                               lo,hi,dx(n,:), &
                                               the_bc_level(n)%adv_bc_level_array(i,:,:,:))
          case (3)
             nzp => dataptr(vel_bc_n(n,3),i)
             t3p => dataptr(vel_bc_t(n,3),i)
             t4p => dataptr(vel_bc_t(n,4),i)
             t5p => dataptr(vel_bc_t(n,5),i)
             t6p => dataptr(vel_bc_t(n,6),i)
             call modify_traction_bcs_3d(nxp(:,:,:,1),nyp(:,:,:,1),nzp(:,:,:,1),ng_n, &
                                               t1p(:,:,:,1),t2p(:,:,:,1),t3p(:,:,:,1), &
                                               t4p(:,:,:,1),t5p(:,:,:,1),t6p(:,:,:,1),ng_t, &
                                               lo,hi,dx(n,:), &
                                               the_bc_level(n)%adv_bc_level_array(i,:,:,:))

          end select
       end do
    end do
    
  end subroutine modify_traction_bcs

  subroutine modify_traction_bcs_2d(vel_bc_nx,vel_bc_ny,ng_n, &
                                          vel_bc_tyx,vel_bc_txy,ng_t,lo,hi,dx,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng_n,ng_t
    real(kind=dp_t), intent(inout) ::  vel_bc_nx(lo(1)-ng_n:,lo(2)-ng_n:)
    real(kind=dp_t), intent(inout) ::  vel_bc_ny(lo(1)-ng_n:,lo(2)-ng_n:)
    real(kind=dp_t), intent(inout) :: vel_bc_tyx(lo(1)-ng_t:,lo(2)-ng_t:)
    real(kind=dp_t), intent(inout) :: vel_bc_txy(lo(1)-ng_t:,lo(2)-ng_t:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ) :: bc(:,:,:)

    ! local
    integer :: i,j
    real(kind=dp_t) :: x,y

    !!!!!!!!!!!!!!!!!!!!!!!!
    ! transverse velocities
    !!!!!!!!!!!!!!!!!!!!!!!!

    ! yvel, lo x-faces
    ! subtract dvx/dy
    if (bc(1,1,2) .eq. DIR_TRACT) then
       x = prob_lo(1)
       do j=lo(2),hi(2)+1
          y = prob_lo(2) + dble(j)*dx(2)
          vel_bc_tyx(lo(1),j) = vel_bc_tyx(lo(1),j) &
               - (vel_bc_nx(lo(1),j)-vel_bc_nx(lo(1),j-1)) / dx(2)
       end do
    end if

    ! yvel, hi x-faces
    ! subtract dvx/dy
    if (bc(1,2,2) .eq. DIR_TRACT) then
       x = prob_hi(1)
       do j=lo(2),hi(2)+1
          y = prob_lo(2) + dble(j)*dx(2)
          vel_bc_tyx(hi(1)+1,j) = vel_bc_tyx(hi(1)+1,j) &
               - (vel_bc_nx(hi(1)+1,j)-vel_bc_nx(hi(1)+1,j-1)) / dx(2)
       end do
    end if

    ! xvel, lo y-faces
    ! subtract dvy/dx
    if (bc(2,1,1) .eq. DIR_TRACT) then
       y = prob_lo(2)
       do i=lo(1),hi(1)+1
          x = prob_lo(1) + dble(i)*dx(1)
          vel_bc_txy(i,lo(2)) = vel_bc_txy(i,lo(2)) &
               - (vel_bc_ny(i,lo(2))-vel_bc_ny(i-1,lo(2))) / dx(1)

       end do
    end if

    ! xvel, hi y-faces
    ! subtract dvy/dx
    if (bc(2,2,1) .eq. DIR_TRACT) then
       y = prob_hi(2)
       do i=lo(1),hi(1)+1
          x = prob_lo(1) + dble(i)*dx(1)
          vel_bc_txy(i,hi(2)+1) = vel_bc_txy(i,hi(2)+1) &
               - (vel_bc_ny(i,hi(2)+1)-vel_bc_ny(i-1,hi(2)+1)) / dx(1)
       end do
    end if

  end subroutine modify_traction_bcs_2d

  subroutine modify_traction_bcs_3d(vel_bc_nx,vel_bc_ny,vel_bc_nz,ng_n, &
                                          vel_bc_tyx,vel_bc_tzx,vel_bc_txy,vel_bc_tzy, &
                                          vel_bc_txz,vel_bc_tyz,ng_t,lo,hi,dx,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng_n,ng_t
    real(kind=dp_t), intent(inout) ::  vel_bc_nx(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:)
    real(kind=dp_t), intent(inout) ::  vel_bc_ny(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:)
    real(kind=dp_t), intent(inout) ::  vel_bc_nz(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:)
    real(kind=dp_t), intent(inout) :: vel_bc_tyx(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
    real(kind=dp_t), intent(inout) :: vel_bc_tzx(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
    real(kind=dp_t), intent(inout) :: vel_bc_txy(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
    real(kind=dp_t), intent(inout) :: vel_bc_tzy(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
    real(kind=dp_t), intent(inout) :: vel_bc_txz(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
    real(kind=dp_t), intent(inout) :: vel_bc_tyz(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ) :: bc(:,:,:)

    ! local
    integer :: i,j,k
    real(kind=dp_t) :: x,y,z

    !!!!!!!!!!!!!!!!!!!!!!!!
    ! transverse velocities
    !!!!!!!!!!!!!!!!!!!!!!!!
    
    ! yvel, lo x-faces
    ! subtract dvx/dy
    if (bc(1,1,2) .eq. DIR_TRACT) then
       x = prob_lo(1)
       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do j=lo(2),hi(2)+1
             y = prob_lo(2) + dble(j)*dx(2)
             vel_bc_tyx(lo(1),j,k) = vel_bc_tyx(lo(1),j,k) &
                  - (vel_bc_nx(lo(1),j,k)-vel_bc_nx(lo(1),j-1,k)) / dx(2)
          end do
       end do

    end if

    ! yvel, hi x-faces
    ! subtract dvx/dy
    if (bc(1,2,2) .eq. DIR_TRACT) then
       x = prob_hi(1)
       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do j=lo(2),hi(2)+1
             y = prob_lo(2) + dble(j)*dx(2)
             vel_bc_tyx(hi(1)+1,j,k) = vel_bc_tyx(hi(1)+1,j,k) &
                  - (vel_bc_nx(hi(1)+1,j,k)-vel_bc_nx(hi(1)+1,j-1,k)) / dx(2)
          end do
       end do

    end if

    ! zvel, lo x-faces
    ! subtract dvx/dz
    if (bc(1,1,3) .eq. DIR_TRACT) then
       x = prob_lo(1)
       do k=lo(3),hi(3)+1
          z = prob_lo(3) + dble(k)*dx(3)
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             vel_bc_tzx(lo(1),j,k) = vel_bc_tzx(lo(1),j,k) &
                  - (vel_bc_nx(lo(1),j,k)-vel_bc_nx(lo(1),j,k-1)) / dx(3)
          end do
       end do

    end if

    ! zvel, hi x-faces
    ! subtract dvx/dz
    if (bc(1,2,3) .eq. DIR_TRACT) then
       x = prob_hi(1)
       do k=lo(3),hi(3)+1
          z = prob_lo(3) + dble(k)*dx(3)
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             vel_bc_tzx(hi(1)+1,j,k) = vel_bc_tzx(hi(1)+1,j,k) &
                  - (vel_bc_nx(hi(1)+1,j,k)-vel_bc_nx(hi(1)+1,j,k-1)) / dx(3)
          end do
       end do

    end if

    ! xvel, lo y-faces
    ! subtract dvy/dx
    if (bc(2,1,1) .eq. DIR_TRACT) then
       y = prob_lo(2)
       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do i=lo(1),hi(1)+1
             x = prob_lo(1) + dble(i)*dx(1)
             vel_bc_txy(i,lo(2),k) = vel_bc_txy(i,lo(2),k) &
                  - (vel_bc_ny(i,lo(2),k)-vel_bc_ny(i-1,lo(2),k)) / dx(1)
          end do
       end do

    end if

    ! xvel, hi y-faces
    ! subtract dvy/dx
    if (bc(2,2,1) .eq. DIR_TRACT) then
       y = prob_hi(2)
       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do i=lo(1),hi(1)+1
             x = prob_lo(1) + dble(i)*dx(1)
             vel_bc_txy(i,hi(2)+1,k) = vel_bc_txy(i,hi(2)+1,k) &
                  - (vel_bc_ny(i,hi(2)+1,k)-vel_bc_ny(i-1,hi(2)+1,k)) / dx(1)
          end do
       end do

    end if

    ! zvel, lo y-faces
    ! subtract dvy/dz
    if (bc(2,1,3) .eq. DIR_TRACT) then
       y = prob_lo(2)
       do k=lo(3),hi(3)+1
          z = prob_lo(3) + dble(k)*dx(3)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             vel_bc_tzy(i,lo(2),k) = vel_bc_tzy(i,lo(2),k) &
                  - (vel_bc_ny(i,lo(2),k)-vel_bc_ny(i,lo(2),k-1)) / dx(3)
          end do
       end do

    end if

    ! zvel, hi y-faces
    ! subtract dvy/dz
    if (bc(2,2,3) .eq. DIR_TRACT) then
       y = prob_hi(2)
       do k=lo(3),hi(3)+1
          z = prob_lo(3) + dble(k)*dx(3)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             vel_bc_tzy(i,hi(2)+1,k) = vel_bc_tzy(i,hi(2)+1,k) &
                  - (vel_bc_ny(i,hi(2)+1,k)-vel_bc_ny(i,hi(2)+1,k-1)) / dx(3)
          end do
       end do

    end if

    ! xvel, lo z-faces
    ! subtract dvz/dx
    if (bc(3,1,1) .eq. DIR_TRACT) then
       z = prob_lo(3)
       do i=lo(1),hi(1)+1
          x = prob_lo(1) + dble(i)*dx(1)
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             vel_bc_txz(i,j,lo(3)) = vel_bc_txz(i,j,lo(3)) &
                  - (vel_bc_nz(i,j,lo(3))-vel_bc_nz(i-1,j,lo(3))) / dx(1)
          end do
       end do

    end if

    ! xvel, hi z-faces
    ! subtract dvz/dx
    if (bc(3,2,1) .eq. DIR_TRACT) then
       z = prob_hi(3)
       do i=lo(1),hi(1)+1
          x = prob_lo(1) + dble(i)*dx(1)
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             vel_bc_txz(i,j,hi(3)+1) = vel_bc_txz(i,j,hi(3)+1) &
                  - (vel_bc_nz(i,j,hi(3)+1)-vel_bc_nz(i-1,j,hi(3)+1)) / dx(1)
          end do
       end do

    end if

    ! yvel, lo z-faces
    ! subtract dvz/dy
    if (bc(3,1,2) .eq. DIR_TRACT) then
       z = prob_lo(3)
       do j=lo(2),hi(2)+1
          y = prob_lo(2) + dble(j)*dx(2)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             vel_bc_tyz(i,j,lo(3)) = vel_bc_tyz(i,j,lo(3)) &
                  - (vel_bc_nz(i,j,lo(3))-vel_bc_nz(i,j-1,lo(3))) / dx(2)
          end do
       end do

    end if

    ! yvel, hi z-faces
    ! subtract dvz/dy
    if (bc(3,2,2) .eq. DIR_TRACT) then
       z = prob_hi(3)
       do j=lo(2),hi(2)+1
          y = prob_lo(2) + dble(j)*dx(2)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             vel_bc_tyz(i,j,hi(3)+1) = vel_bc_tyz(i,j,hi(3)+1) &
                  - (vel_bc_nz(i,j,hi(3)+1)-vel_bc_nz(i,j-1,hi(3)+1)) / dx(2)
          end do
       end do

    end if

  end subroutine modify_traction_bcs_3d

end module set_inhomogeneous_vel_bcs_module
