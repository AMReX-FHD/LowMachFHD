module multifab_physbc_bds_module

  use ml_layout_module
  use multifab_module
  use define_bc_module
  use bc_module
  use bl_error_module
  use inhomogeneous_bc_val_module
  use probin_common_module, only: prob_lo, prob_hi

  implicit none

  private

  public :: multifab_physbc_bds

contains

  subroutine multifab_physbc_bds(s,start_scomp,start_bccomp,num_comp,the_bc_level, &
                                 time_in,dx_in,prob_lo_in,prob_hi_in,increment_bccomp_in)

    ! this fills 2 ghost cells for conetration using cubic interpolation
    ! here, we extrapolate to the ghost cell-centers, not the physical boundary location

    type(multifab) , intent(inout) :: s
    integer        , intent(in   ) :: start_scomp,start_bccomp,num_comp
    type(bc_level) , intent(in   ) :: the_bc_level
    real(kind=dp_t), intent(in   ), optional :: time_in,dx_in(:),prob_lo_in(:),prob_hi_in(:)
    logical        ,  intent(in  ), optional :: increment_bccomp_in
   
    ! Local
    integer :: lo(get_dim(s)),hi(get_dim(s))
    integer :: i,ng,dm,scomp,bccomp
    logical :: increment_bccomp

    real(kind=dp_t), allocatable :: dx(:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)

    ng = nghost(s)
    dm = get_dim(s)

    allocate(dx(dm))

    if (present(dx_in)) then
       dx(1:dm) = dx_in(1:dm)
    else
       dx = 1.d0
    end if

    increment_bccomp = .true.
    if (present(increment_bccomp_in)) then
       increment_bccomp = increment_bccomp_in
    end if
    
    do i=1,nfabs(s)
       sp => dataptr(s,i)
       lo = lwb(get_box(s,i))
       hi = upb(get_box(s,i))
       do scomp=start_scomp,start_scomp+num_comp-1
          if (increment_bccomp) then
             bccomp = start_bccomp + (scomp-start_scomp)
          else
             bccomp = start_bccomp
          end if
          select case (dm)
          case (2)
             call physbc_bds_2d(sp(:,:,1,scomp), lo, hi, ng, &
                                the_bc_level%adv_bc_level_array(i,:,:,bccomp),bccomp, dx)
          case (3)
             call physbc_bds_3d(sp(:,:,:,scomp), lo, hi, ng, &
                                the_bc_level%adv_bc_level_array(i,:,:,bccomp),bccomp, dx)
          end select
       end do
    end do

    deallocate(dx)
 
  end subroutine multifab_physbc_bds

  subroutine physbc_bds_2d(s,lo,hi,ng,bc,bccomp,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) ::    s(lo(1)-ng:,lo(2)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer :: i,j
    real(kind=dp_t) :: x,y,fac

    fac = 1.d0/29.d0

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. FOEXTRAP .or. bc(1,1) .eq. HOEXTRAP) then
       ! cubic extrapolation using 3 interior cell averages and Neumann bc
       x = prob_lo(1)
       do j=lo(2)-ng,hi(2)+ng
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          s(lo(1)-1,j) = fac*( -39.d0*dx(1)*inhomogeneous_bc_val_2d(bccomp,x,y) &
                              + 15.d0*s(lo(1),j) +  18.d0*s(lo(1)+1,j) -  4.d0*s(lo(1)+2,j))
          s(lo(1)-2,j) = fac*(-264.d0*dx(1)*inhomogeneous_bc_val_2d(bccomp,x,y) &
                              -242.d0*s(lo(1),j) + 336.d0*s(lo(1)+1,j) - 65.d0*s(lo(1)+2,j))
       end do
    else if (bc(1,1) .eq. EXT_DIR) then
       ! dirichlet condition obtained from inhomogeneous_bc_val
       x = prob_lo(1)
       do j=lo(2)-ng,hi(2)+ng
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          s(lo(1)-1,j) =  -26.d0*inhomogeneous_bc_val_2d(bccomp,x,y) &
                         + 44.d0*s(lo(1),j) -  20.d0*s(lo(1)+1,j) +  3.d0*s(lo(1)+2,j)
          s(lo(1)-2,j) = -176.d0*inhomogeneous_bc_val_2d(bccomp,x,y) &
                        + 286.d0*s(lo(1),j) - 128.d0*s(lo(1)+1,j) + 19.d0*s(lo(1)+2,j)
       end do
    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_bds_2d: bc(1,1) =',bc(1,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. FOEXTRAP .or. bc(1,2) .eq. HOEXTRAP) then
       ! cubic extrapolation using 3 interior cell averages and Neumann bc
       x = prob_hi(1)
       do j=lo(2)-ng,hi(2)+ng
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          s(hi(1)+1,j) = fac*(  39.d0*dx(1)*inhomogeneous_bc_val_2d(bccomp,x,y) &
                              + 15.d0*s(hi(1),j) +  18.d0*s(hi(1)-1,j) -  4.d0*s(hi(1)-2,j))
          s(hi(1)+2,j) = fac*( 264.d0*dx(1)*inhomogeneous_bc_val_2d(bccomp,x,y) &
                              -242.d0*s(hi(1),j) + 336.d0*s(hi(1)-1,j) - 65.d0*s(hi(1)-2,j))
       end do
    else if (bc(1,2) .eq. EXT_DIR) then
       ! dirichlet condition obtained from inhomogeneous_bc_val
       x = prob_hi(1)
       do j=lo(2)-ng,hi(2)+ng
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          s(hi(1)+1,j) =  -26.d0*inhomogeneous_bc_val_2d(bccomp,x,y) &
                         + 44.d0*s(hi(1),j) -  20.d0*s(hi(1)-1,j) +  3.d0*s(hi(1)-2,j)
          s(hi(1)+2,j) = -176.d0*inhomogeneous_bc_val_2d(bccomp,x,y) &
                        + 286.d0*s(hi(1),j) - 128.d0*s(hi(1)-1,j) + 19.d0*s(hi(1)-2,j)
       end do
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_bds_2d: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. FOEXTRAP .or. bc(2,1) .eq. HOEXTRAP) then
       ! cubic extrapolation using 3 interior cell averages and Neumann bc
       y = prob_lo(2)
       do i=lo(1)-ng,hi(1)+ng
          x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
          s(i,lo(2)-1) = fac*( -39.d0*dx(2)*inhomogeneous_bc_val_2d(bccomp,x,y) &
                              + 15.d0*s(i,lo(2)) +  18.d0*s(i,lo(2)+1) -  4.d0*s(i,lo(2)+2))
          s(i,lo(2)-2) = fac*(-264.d0*dx(2)*inhomogeneous_bc_val_2d(bccomp,x,y) &
                              -242.d0*s(i,lo(2)) + 336.d0*s(i,lo(2)+1) - 65.d0*s(i,lo(2)+2))
       end do
    else if (bc(2,1) .eq. EXT_DIR) then
       ! cubic extrapolation using 3 interior cell averages and Dirichlet bc
       y = prob_lo(2)
       do i=lo(1)-ng,hi(1)+ng
          x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
          s(i,lo(2)-1) =  -26.d0*inhomogeneous_bc_val_2d(bccomp,x,y) &
                         + 44.d0*s(i,lo(2)) -  20.d0*s(i,lo(2)+1) +  3.d0*s(i,lo(2)+2)
          s(i,lo(2)-2) = -176.d0*inhomogeneous_bc_val_2d(bccomp,x,y) &
                        + 286.d0*s(i,lo(2)) - 128.d0*s(i,lo(2)+1) + 19.d0*s(i,lo(2)+2)
       end do
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_bds_2d: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. FOEXTRAP .or. bc(2,2) .eq. HOEXTRAP) then
       ! cubic extrapolation using 3 interior cell averages and Neumann bc
       y = prob_hi(2)
       do i=lo(1)-ng,hi(1)+ng
          x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
          s(i,hi(2)+1) = fac*( 39.d0*dx(2)*inhomogeneous_bc_val_2d(bccomp,x,y) &
                              + 15.d0*s(i,hi(2)) +  18.d0*s(i,hi(2)-1) -  4.d0*s(i,hi(2)-2))
          s(i,hi(2)+2) = fac*(264.d0*dx(2)*inhomogeneous_bc_val_2d(bccomp,x,y) &
                              -242.d0*s(i,hi(2)) + 336.d0*s(i,hi(2)-1) - 65.d0*s(i,hi(2)-2))
       end do
    else if (bc(2,2) .eq. EXT_DIR) then
       ! cubic extrapolation using 3 interior cell averages and Dirichlet bc
       y = prob_hi(2)
       do i=lo(1)-ng,hi(1)+ng
          x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
          s(i,hi(2)+1) =  -26.d0*inhomogeneous_bc_val_2d(bccomp,x,y) &
                         + 44.d0*s(i,hi(2)) -  20.d0*s(i,hi(2)-1) +  3.d0*s(i,hi(2)-2)
          s(i,hi(2)+2) = -176.d0*inhomogeneous_bc_val_2d(bccomp,x,y) &
                        + 286.d0*s(i,hi(2)) - 128.d0*s(i,hi(2)-1) + 19.d0*s(i,hi(2)-2)
       end do
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_bds_2d: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine physbc_bds_2d

  subroutine physbc_bds_3d(s,lo,hi,ng,bc,bccomp,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) ::    s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer :: i,j,k
    real(kind=dp_t) :: x,y,z,fac

    fac = 1.d0/29.d0

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. FOEXTRAP .or. bc(1,1) .eq. HOEXTRAP) then
       ! cubic extrapolation using 3 interior cell averages and Neumann bc
       x = prob_lo(1)
       do k=lo(3)-ng,hi(3)+ng
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do j=lo(2)-ng,hi(2)+ng
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             s(lo(1)-1,j,k) = fac*( -39.d0*dx(1)*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                                   + 15.d0*s(lo(1),j,k) +  18.d0*s(lo(1)+1,j,k) -  4.d0*s(lo(1)+2,j,k))
             s(lo(1)-2,j,k) = fac*(-264.d0*dx(1)*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                                   -242.d0*s(lo(1),j,k) + 336.d0*s(lo(1)+1,j,k) - 65.d0*s(lo(1)+2,j,k))
          end do
       end do
    else if (bc(1,1) .eq. EXT_DIR) then
       ! cubic extrapolation using 3 interior cell averages and Dirichlet bc
       x = prob_lo(1)
       do k=lo(3)-ng,hi(3)+ng
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do j=lo(2)-ng,hi(2)+ng
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             s(lo(1)-1,j,k) =  -26.d0*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                              + 44.d0*s(lo(1),j,k) -  20.d0*s(lo(1)+1,j,k) +  3.d0*s(lo(1)+2,j,k)
             s(lo(1)-2,j,k) = -176.d0*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                             + 286.d0*s(lo(1),j,k) - 128.d0*s(lo(1)+1,j,k) + 19.d0*s(lo(1)+2,j,k)
          end do
       end do
    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_bds_3d: bc(1,1) =',bc(1,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. FOEXTRAP .or. bc(1,2) .eq. HOEXTRAP) then
       ! cubic extrapolation using 3 interior cell averages and Neumann bc
       x = prob_hi(1)
       do k=lo(3)-ng,hi(3)+ng
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do j=lo(2)-ng,hi(2)+ng
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             s(hi(1)+1,j,k) = fac*(  39.d0*dx(1)*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                                   + 15.d0*s(hi(1),j,k) +  18.d0*s(hi(1)-1,j,k) -  4.d0*s(hi(1)-2,j,k))
             s(hi(1)+2,j,k) = fac*( 264.d0*dx(1)*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                                   -242.d0*s(hi(1),j,k) + 336.d0*s(hi(1)-1,j,k) - 65.d0*s(hi(1)-2,j,k))
          end do
       end do
    else if (bc(1,2) .eq. EXT_DIR) then
       ! cubic extrapolation using 3 interior cell averages and Dirichlet bc
       x = prob_hi(1)
       do k=lo(3)-ng,hi(3)+ng
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do j=lo(2)-ng,hi(2)+ng
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             s(hi(1)+1,j,k) =  -26.d0*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                              + 44.d0*s(hi(1),j,k) -  20.d0*s(hi(1)-1,j,k) +  3.d0*s(hi(1)-2,j,k)
             s(hi(1)+2,j,k) = -176.d0*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                             + 286.d0*s(hi(1),j,k) - 128.d0*s(hi(1)-1,j,k) + 19.d0*s(hi(1)-2,j,k)
          end do
       end do
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_bds_3d: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. FOEXTRAP .or. bc(2,1) .eq. HOEXTRAP) then
       ! cubic extrapolation using 3 interior cell averages and Neumann bc
       y = prob_lo(2)
       do k=lo(3)-ng,hi(3)+ng
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do i=lo(1)-ng,hi(1)+ng
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             s(i,lo(2)-1,k) = fac*( -39.d0*dx(2)*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                                   + 15.d0*s(i,lo(2),k) +  18.d0*s(i,lo(2)+1,k) -  4.d0*s(i,lo(2)+2,k))
             s(i,lo(2)-2,k) = fac*(-264.d0*dx(2)*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                                   -242.d0*s(i,lo(2),k) + 336.d0*s(i,lo(2)+1,k) - 65.d0*s(i,lo(2)+2,k))
          end do
       end do
    else if (bc(2,1) .eq. EXT_DIR) then
       ! cubic extrapolation using 3 interior cell averages and Dirichlet bc
       y = prob_lo(2)
       do k=lo(3)-ng,hi(3)+ng
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do i=lo(1)-ng,hi(1)+ng
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             s(i,lo(2)-1,k) =  -26.d0*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                              + 44.d0*s(i,lo(2),k) -  20.d0*s(i,lo(2)+1,k) +  3.d0*s(i,lo(2)+2,k)
             s(i,lo(2)-2,k) = -176.d0*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                             + 286.d0*s(i,lo(2),k) - 128.d0*s(i,lo(2)+1,k) + 19.d0*s(i,lo(2)+2,k)
          end do
       end do
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_bds_3d: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. FOEXTRAP .or. bc(2,2) .eq. HOEXTRAP) then
       ! cubic extrapolation using 3 interior cell averages and Neumann bc
       y = prob_hi(2)
       do k=lo(3)-ng,hi(3)+ng
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do i=lo(1)-ng,hi(1)+ng
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             s(i,hi(2)+1,k) = fac*(  39.d0*dx(2)*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                                   + 15.d0*s(i,hi(2),k) +  18.d0*s(i,hi(2)-1,k) -  4.d0*s(i,hi(2)-2,k))
             s(i,hi(2)+2,k) = fac*( 264.d0*dx(2)*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                                   -242.d0*s(i,hi(2),k) + 336.d0*s(i,hi(2)-1,k) - 65.d0*s(i,hi(2)-2,k))
          end do
       end do
    else if (bc(2,2) .eq. EXT_DIR) then
       ! cubic extrapolation using 3 interior cell averages and Dirichlet bc
       y = prob_hi(2)
       do k=lo(3)-ng,hi(3)+ng
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do i=lo(1)-ng,hi(1)+ng
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             s(i,hi(2)+1,k) =  -26.d0*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                              + 44.d0*s(i,hi(2),k) -  20.d0*s(i,hi(2)-1,k) +  3.d0*s(i,hi(2)-2,k)
             s(i,hi(2)+2,k) = -176.d0*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                             + 286.d0*s(i,hi(2),k) - 128.d0*s(i,hi(2)-1,k) + 19.d0*s(i,hi(2)-2,k)
          end do
       end do
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_bds_3d: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,1) .eq. FOEXTRAP .or. bc(3,1) .eq. HOEXTRAP) then
       ! cubic extrapolation using 3 interior cell averages and Neumann bc
       z = prob_lo(3)
       do j=lo(2)-ng,hi(2)+ng
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          do i=lo(1)-ng,hi(1)+ng
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             s(i,j,lo(3)-1) = fac*( -39.d0*dx(3)*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                                   + 15.d0*s(i,j,lo(3)) +  18.d0*s(i,j,lo(3)+1) -  4.d0*s(i,j,lo(3)+2))
             s(i,j,lo(3)-2) = fac*(-264.d0*dx(3)*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                                   -242.d0*s(i,j,lo(3)) + 336.d0*s(i,j,lo(3)+1) - 65.d0*s(i,j,lo(3)+2))
          end do
       end do
    else if (bc(3,1) .eq. EXT_DIR) then
       ! cubic extrapolation using 3 interior cell averages and Dirichlet bc
       z = prob_lo(3)
       do j=lo(2)-ng,hi(2)+ng
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          do i=lo(1)-ng,hi(1)+ng
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             s(i,j,lo(3)-1) =  -26.d0*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                              + 44.d0*s(i,j,lo(3)) -  20.d0*s(i,j,lo(3)+1) +  3.d0*s(i,j,lo(3)+2)
             s(i,j,lo(3)-2) = -176.d0*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                             + 286.d0*s(i,j,lo(3)) - 128.d0*s(i,j,lo(3)+1) + 19.d0*s(i,j,lo(3)+2)
          end do
       end do
    else if (bc(3,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_bds_3d: bc(3,1) =',bc(3,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,2) .eq. FOEXTRAP .or. bc(3,2) .eq. HOEXTRAP) then
       ! cubic extrapolation using 3 interior cell averages and Neumann bc
       z = prob_hi(3)
       do j=lo(2)-ng,hi(2)+ng
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          do i=lo(1)-ng,hi(1)+ng
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             s(i,j,hi(3)+1) = fac*(  39.d0*dx(3)*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                                   + 15.d0*s(i,j,hi(3)) +  18.d0*s(i,j,hi(3)-1) -  4.d0*s(i,j,hi(3)-2))
             s(i,j,hi(3)+2) = fac*( 264.d0*dx(3)*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                                   -242.d0*s(i,j,hi(3)) + 336.d0*s(i,j,hi(3)-1) - 65.d0*s(i,j,hi(3)-2))
          end do
       end do
    else if (bc(3,2) .eq. EXT_DIR) then
       ! cubic extrapolation using 3 interior cell averages and Dirichlet bc
       z = prob_hi(3)
       do j=lo(2)-ng,hi(2)+ng
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          do i=lo(1)-ng,hi(1)+ng
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             s(i,j,hi(3)+1) =  -26.d0*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                              + 44.d0*s(i,j,hi(3)) -  20.d0*s(i,j,hi(3)-1) +  3.d0*s(i,j,hi(3)-2)
             s(i,j,hi(3)+2) = -176.d0*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                             + 286.d0*s(i,j,hi(3)) - 128.d0*s(i,j,hi(3)-1) + 19.d0*s(i,j,hi(3)-2)
          end do
       end do
    else if (bc(3,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_bds_3d: bc(3,2) =',bc(3,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine physbc_bds_3d

end module multifab_physbc_bds_module
