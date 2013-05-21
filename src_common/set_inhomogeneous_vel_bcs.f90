module set_inhomogeneous_vel_bcs_module

  use ml_layout_module
  use multifab_module
  use define_bc_module
  use bc_module
  use inhomogeneous_bc_val_module
  use probin_common_module, only: prob_lo, prob_hi

  implicit none

  private

  public :: set_inhomogeneous_vel_bcs

contains

  subroutine set_inhomogeneous_vel_bcs(mla,vel_bc,dx,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: vel_bc(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local
    integer :: n,nlevs,i,dm,ng_v
    integer :: lo(mla%dim),hi(mla%dim)
    real(kind=dp_t), pointer :: vxp(:,:,:,:)
    real(kind=dp_t), pointer :: vyp(:,:,:,:)
    real(kind=dp_t), pointer :: vzp(:,:,:,:)

    nlevs = mla%nlevel
    dm = mla%dim

    ng_v = vel_bc(1,1)%ng

    do n=1,nlevs
       do i=1,nfabs(vel_bc(n,1))
          vxp => dataptr(vel_bc(n,1),i)
          vyp => dataptr(vel_bc(n,2),i)
          lo = lwb(get_box(vel_bc(n,1),i))
          hi = upb(get_box(vel_bc(n,1),i))
          select case (dm)
          case (2)
             call set_inhomogeneous_vel_bcs_2d(vxp(:,:,1,:),vyp(:,:,1,:),ng_v, &
                                               lo,hi,dx(n,:), &
                                               the_bc_level(n)%adv_bc_level_array(i,:,:,:))
          case (3)
             vzp => dataptr(vel_bc(n,3),i)
             call set_inhomogeneous_vel_bcs_3d(vxp(:,:,:,:),vyp(:,:,:,:),vzp(:,:,:,:),ng_v, &
                                               lo,hi,dx(n,:), &
                                               the_bc_level(n)%adv_bc_level_array(i,:,:,:))

          end select
       end do
    end do
    
  end subroutine set_inhomogeneous_vel_bcs

  subroutine set_inhomogeneous_vel_bcs_2d(v_bcx,v_bcy,ng_v,lo,hi,dx,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng_v
    real(kind=dp_t), intent(inout) :: v_bcx(lo(1)-ng_v:,lo(2)-ng_v:,:)
    real(kind=dp_t), intent(inout) :: v_bcy(lo(1)-ng_v:,lo(2)-ng_v:,:)
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
          v_bcx(lo(1),j,1) = inhomogeneous_bc_val_2d(1,x,y)
       end do
    end if

    ! xvel, hi x-faces
    if (bc(1,2,1) .eq. DIR_VEL) then
       x = prob_hi(1)
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          v_bcx(hi(1)+1,j,1) = inhomogeneous_bc_val_2d(1,x,y)
       end do
    end if

    ! yvel, lo y-faces
    if (bc(2,1,2) .eq. DIR_VEL) then
       y = prob_lo(2)
       do i=lo(1),hi(1)
          x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
          v_bcy(i,lo(2),2) = inhomogeneous_bc_val_2d(2,x,y)
       end do
    end if

    ! yvel, hi y-faces
    if (bc(2,2,2) .eq. DIR_VEL) then
       y = prob_hi(2)
       do i=lo(1),hi(1)
          x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
          v_bcy(i,hi(2)+1,2) = inhomogeneous_bc_val_2d(2,x,y)
       end do
    end if

    !!!!!!!!!!!!!!!!!!!!
    ! transverse velocities
    !!!!!!!!!!!!!!!!!!!!

    ! xvel, lo y-faces
    if (bc(2,1,1) .eq. DIR_VEL .or. bc(2,1,1) .eq. DIR_TRACT) then
       y = prob_lo(2)
       do i=lo(1),hi(1)+1
          x = prob_lo(1) + dble(i)*dx(1)
          v_bcx(i,lo(2),2) = inhomogeneous_bc_val_2d(1,x,y)
       end do
    end if

    ! xvel, hi y-faces
    if (bc(2,2,1) .eq. DIR_VEL .or. bc(2,2,1) .eq. DIR_TRACT) then
       y = prob_hi(2)
       do i=lo(1),hi(1)+1
          x = prob_lo(1) + dble(i)*dx(1)
          v_bcx(i,hi(2),2) = inhomogeneous_bc_val_2d(1,x,y)
       end do
    end if

    ! yvel, lo x-faces
    if (bc(1,1,2) .eq. DIR_VEL .or. bc(1,1,2) .eq. DIR_TRACT) then
       x = prob_lo(1)
       do j=lo(2),hi(2)+1
          y = prob_lo(2) + dble(j)*dx(2)
          v_bcy(lo(1),j,1) = inhomogeneous_bc_val_2d(2,x,y)
       end do
    end if

    ! yvel, hi x-faces
    if (bc(1,2,2) .eq. DIR_VEL .or. bc(1,2,2) .eq. DIR_TRACT) then
       x = prob_hi(1)
       do j=lo(2),hi(2)+1
          y = prob_lo(2) + dble(j)*dx(2)
          v_bcy(hi(1),j,1) = inhomogeneous_bc_val_2d(2,x,y)
       end do
    end if

  end subroutine set_inhomogeneous_vel_bcs_2d

  subroutine set_inhomogeneous_vel_bcs_3d(v_bcx,v_bcy,v_bcz,ng_v,lo,hi,dx,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng_v
    real(kind=dp_t), intent(inout) :: v_bcx(lo(1)-ng_v:,lo(2)-ng_v:,lo(3)-ng_v:,:)
    real(kind=dp_t), intent(inout) :: v_bcy(lo(1)-ng_v:,lo(2)-ng_v:,lo(3)-ng_v:,:)
    real(kind=dp_t), intent(inout) :: v_bcz(lo(1)-ng_v:,lo(2)-ng_v:,lo(3)-ng_v:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ) :: bc(:,:,:)

    ! local
    integer :: i,j,k
    real(kind=dp_t) :: x,y,z

  end subroutine set_inhomogeneous_vel_bcs_3d

end module set_inhomogeneous_vel_bcs_module
