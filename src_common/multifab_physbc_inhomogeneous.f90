module multifab_physbc_inhomogeneous_module

  use multifab_module
  use define_bc_module
  use bc_module
  use bl_error_module
  use probin_module

  implicit none

  private

  public ::  multifab_physbc_macvel_inhomogeneous, multifab_physbc_domainvel_inhomogeneous

contains

  subroutine multifab_physbc_macvel_inhomogeneous(s,scomp,bccomp,num_comp,the_bc_level,dx)

    type(multifab) , intent(inout) :: s
    integer        , intent(in   ) :: scomp,bccomp,num_comp
    type(bc_level) , intent(in   ) :: the_bc_level
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! Local
    integer                  :: lo(get_dim(s)),hi(get_dim(s))
    integer                  :: i,ng_s,dm
    real(kind=dp_t), pointer :: sp(:,:,:,:)

    if (scomp .ne. 1) then
       call bl_error('multifab_physbc_macvel_inhomogeneous expects scomp = 1')
    end if

    if (bccomp .gt. get_dim(s)) then
       call bl_error('multifab_physbc_macvel_inhomogeneous expects bccomp <= dm')
    end if

    if (num_comp .ne. 1) then
       call bl_error('multifab_physbc_macvel_inhomogeneous expects num_comp = 1')
    end if

    ng_s = nghost(s)
    dm = get_dim(s)

    do i=1,nfabs(s)
       sp => dataptr(s,i)
       lo = lwb(get_box(s,i))
       hi = upb(get_box(s,i))
       select case (dm)
       case (2)
          call physbc_macvel_inhomogeneous_2d(sp(:,:,1,scomp), ng_s, lo, hi, &
                                              the_bc_level%adv_bc_level_array(i,:,:,bccomp), &
                                              bccomp,dx)
       case (3)
          call physbc_macvel_inhomogeneous_3d(sp(:,:,:,scomp), ng_s, lo, hi, &
                                              the_bc_level%adv_bc_level_array(i,:,:,bccomp), &
                                              bccomp,dx)
       end select
    end do

  end subroutine multifab_physbc_macvel_inhomogeneous

  subroutine physbc_macvel_inhomogeneous_2d(s,ng_s,lo,hi,bc,bccomp,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    real(kind=dp_t), intent(inout) ::             s(lo(1)-ng_s:,lo(2)-ng_s:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer :: i,j
    real(kind=dp_t) :: x,y

    if (bccomp .ne. 1 .and. bccomp .ne. 2) then
       call bl_error('physbc_macvel_inhomogeneous_2d requires bccomp = 1 or 2')
    end if

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel_inhomogeneous
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          x = prob_lo(1)
          do j=lo(2),hi(2)+1
             y = prob_lo(2) + dble(j)*dx(2)
             s(lo(1)-ng_s:lo(1)-1,j) = 2.d0*inhomogeneous_bc_val_2d(bccomp,x,y) - s(lo(1),j)
          end do
       end if
    else if (bc(1,1) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          print *,'physbc_macvel_inhomogeneous_2d: bc(1,1) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          x = prob_lo(1)
          do j=lo(2),hi(2)+1
             y = prob_lo(2) + dble(j)*dx(2)
             s(lo(1)-ng_s:lo(1)-1,j) = s(lo(1),j) - dx(1)*inhomogeneous_bc_val_2d(bccomp,x,y)
          end do
       end if
    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_inhomogeneous_2d: bc(1,1) =',bc(1,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel_inhomogeneous
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          x = prob_hi(1)
          do j=lo(2),hi(2)+1
             y = prob_lo(2) + dble(j)*dx(2)
             s(hi(1)+1:hi(1)+ng_s,j) = 2.d0*inhomogeneous_bc_val_2d(bccomp,x,y) - s(hi(1),j)
          end do
       end if
    else if (bc(1,2) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          print *,'physbc_macvel_inhomogeneous_2d: bc(1,2) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          x = prob_hi(1)
          do j=lo(2),hi(2)+1
             y = prob_lo(2) + dble(j)*dx(2)
             s(hi(1)+1:hi(1)+ng_s,j) = s(hi(1),j) + dx(1)*inhomogeneous_bc_val_2d(bccomp,x,y)
          end do
       end if
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_inhomogeneous_2d: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          y = prob_lo(2)
          do i=lo(1),hi(1)+1
             x = prob_lo(1) + dble(i)*dx(1)
             s(i,lo(2)-ng_s:lo(2)-1) = 2.d0*inhomogeneous_bc_val_2d(bccomp,x,y) - s(i,lo(2))
          end do
       else if (bccomp .eq. 2) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel_inhomogeneous
       end if
    else if (bc(2,1) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          y = prob_lo(2)
          do i=lo(1),hi(1)+1
             x = prob_lo(1) + dble(i)*dx(1)
             s(i,lo(2)-ng_s:lo(2)-1) = s(i,lo(2)) - dx(1)*inhomogeneous_bc_val_2d(bccomp,x,y)
          end do
       else if (bccomp .eq. 2) then
          print *,'physbc_macvel_inhomogeneous_2d: bc(2,1) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_inhomogeneous_2d: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          y = prob_hi(2)
          do i=lo(1),hi(1)+1
             x = prob_lo(1) + dble(i)*dx(1)
             s(i,hi(2)+1:hi(2)+ng_s) = 2.d0*inhomogeneous_bc_val_2d(bccomp,x,y) - s(i,hi(2))
          end do
       else if (bccomp .eq. 2) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel_inhomogeneous
       end if
    else if (bc(2,2) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          y = prob_hi(2)
          do i=lo(1),hi(1)+1
             x = prob_lo(1) + dble(i)*dx(1)
             s(i,hi(2)+1:hi(2)+ng_s) = s(i,hi(2)) + dx(1)*inhomogeneous_bc_val_2d(bccomp,x,y)
          end do
       else if (bccomp .eq. 2) then
          print *,'physbc_macvel_inhomogeneous_2d: bc(2,2) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_inhomogeneous_2d: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine physbc_macvel_inhomogeneous_2d

  subroutine physbc_macvel_inhomogeneous_3d(s,ng_s,lo,hi,bc,bccomp,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer :: i,j,k
    real(kind=dp_t) :: x,y,z

    if (bccomp .ne. 1 .and. bccomp .ne. 2 .and. bccomp .ne. 3) then
       call bl_error('physbc_macvel_inhomogeneous_3d requires bccomp = 1, 2 or 3')
    end if

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel_inhomogeneous
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          x = prob_lo(1)
          do k=lo(3),hi(3)
             z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
             do j=lo(2),hi(2)+1
                y = prob_lo(3) + dble(j)*dx(2)
                s(lo(1)-ng_s:lo(1)-1,j,k) = 2.d0*inhomogeneous_bc_val_3d(bccomp,x,y,z) - s(lo(1),j,k)
             end do
          end do
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          x = prob_lo(1)
          do k=lo(3),hi(3)+1
             z = prob_lo(3) + dble(k)*dx(3)
             do j=lo(2),hi(2)
                y = prob_lo(3) + (dble(j)+0.5d0)*dx(2)
                s(lo(1)-ng_s:lo(1)-1,j,k) = 2.d0*inhomogeneous_bc_val_3d(bccomp,x,y,z) - s(lo(1),j,k)
             end do
          end do
       end if
    else if (bc(1,1) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          print *,'physbc_macvel_inhomogeneous_3d: bc(1,1) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          x = prob_lo(1)
          do k=lo(3),hi(3)
             z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
             do j=lo(2),hi(2)+1
                y = prob_lo(3) + dble(j)*dx(2)
                s(lo(1)-ng_s:lo(1)-1,j,k) = s(lo(1),j,k) - dx(1)*inhomogeneous_bc_val_3d(bccomp,x,y,z)
             end do
          end do
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          x = prob_lo(1)
          do k=lo(3),hi(3)+1
             z = prob_lo(3) + dble(k)*dx(3)
             do j=lo(2),hi(2)
                y = prob_lo(3) + (dble(j)+0.5d0)*dx(2)
                s(lo(1)-ng_s:lo(1)-1,j,k) = s(lo(1),j,k) - dx(1)*inhomogeneous_bc_val_3d(bccomp,x,y,z)
             end do
          end do
       end if
    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_inhomogeneous_3d: bc(1,1) =',bc(1,1),'for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel_inhomogeneous
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          x = prob_hi(1)
          do k=lo(3),hi(3)
             z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
             do j=lo(2),hi(2)+1
                y = prob_lo(3) + dble(j)*dx(2)
                s(hi(1)+1:hi(1)+ng_s,j,k) = 2.d0*inhomogeneous_bc_val_3d(bccomp,x,y,z) - s(hi(1),j,k)
             end do
          end do
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          x = prob_hi(1)
          do k=lo(3),hi(3)+1
             z = prob_lo(3) + dble(k)*dx(3)
             do j=lo(2),hi(2)
                y = prob_lo(3) + (dble(j)+0.5d0)*dx(2)
                s(hi(1)+1:hi(1)+ng_s,j,k) = 2.d0*inhomogeneous_bc_val_3d(bccomp,x,y,z) - s(hi(1),j,k)
             end do
          end do
       end if
    else if (bc(1,2) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          print *,'physbc_macvel_inhomogeneous_3d: bc(1,2) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          x = prob_hi(1)
          do k=lo(3),hi(3)
             z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
             do j=lo(2),hi(2)+1
                y = prob_lo(3) + dble(j)*dx(2)
                s(hi(1)+1:hi(1)+ng_s,j,k) = s(hi(1),j,k) + dx(1)*inhomogeneous_bc_val_3d(bccomp,x,y,z)
             end do
          end do
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          x = prob_hi(1)
          do k=lo(3),hi(3)+1
             z = prob_lo(3) + dble(k)*dx(3)
             do j=lo(2),hi(2)
                y = prob_lo(3) + (dble(j)+0.5d0)*dx(2)
                s(hi(1)+1:hi(1)+ng_s,j,k) = s(hi(1),j,k) + dx(1)*inhomogeneous_bc_val_3d(bccomp,x,y,z)
             end do
          end do
       end if
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_inhomogeneous_3d: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          y = prob_lo(2)
          do k=lo(3),hi(3)
             z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
             do i=lo(1),hi(1)+1
                x = prob_lo(1) + dble(i)*dx(1)
                s(i,lo(2)-ng_s:lo(2)-1,k) = 2.d0*inhomogeneous_bc_val_3d(bccomp,x,y,z) - s(i,lo(2),k)
             end do
          end do
       else if (bccomp .eq. 2) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel_inhomogeneous
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          y = prob_lo(2)
          do k=lo(3),hi(3)+1
             z = prob_lo(3) + dble(k)*dx(3)
             do i=lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                s(i,lo(2)-ng_s:lo(2)-1,k) = 2.d0*inhomogeneous_bc_val_3d(bccomp,x,y,z) - s(i,lo(2),k)
             end do
          end do
       end if
    else if (bc(2,1) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          y = prob_lo(2)
          do k=lo(3),hi(3)
             z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
             do i=lo(1),hi(1)+1
                x = prob_lo(1) + dble(i)*dx(1)
                s(i,lo(2)-ng_s:lo(2)-1,k) = s(i,lo(2),k) - dx(2)*inhomogeneous_bc_val_3d(bccomp,x,y,z)
             end do
          end do
       else if (bccomp .eq. 2) then
          print *,'physbc_macvel_inhomogeneous_3d: bc(2,1) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          y = prob_lo(2)
          do k=lo(3),hi(3)+1
             z = prob_lo(3) + dble(k)*dx(3)
             do i=lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                s(i,lo(2)-ng_s:lo(2)-1,k) = s(i,lo(2),k) - dx(2)*inhomogeneous_bc_val_3d(bccomp,x,y,z)
             end do
          end do
       end if
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_inhomogeneous_3d: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          y = prob_hi(2)
          do k=lo(3),hi(3)
             z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
             do i=lo(1),hi(1)+1
                x = prob_lo(1) + dble(i)*dx(1)
                s(i,hi(2)+1:hi(2)+ng_s,k) = 2.d0*inhomogeneous_bc_val_3d(bccomp,x,y,z) - s(i,hi(2),k)
             end do
          end do
       else if (bccomp .eq. 2) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel_inhomogeneous
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          y = prob_hi(2)
          do k=lo(3),hi(3)+1
             z = prob_lo(3) + dble(k)*dx(3)
             do i=lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                s(i,hi(2)+1:hi(2)+ng_s,k) = 2.d0*inhomogeneous_bc_val_3d(bccomp,x,y,z) - s(i,hi(2),k)
             end do
          end do
       end if
    else if (bc(2,2) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          y = prob_hi(2)
          do k=lo(3),hi(3)
             z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
             do i=lo(1),hi(1)+1
                x = prob_lo(1) + dble(i)*dx(1)
                s(i,hi(2)+1:hi(2)+ng_s,k) = s(i,hi(2),k) + dx(2)*inhomogeneous_bc_val_3d(bccomp,x,y,z)
             end do
          end do
       else if (bccomp .eq. 2) then
          print *,'physbc_macvel_inhomogeneous_3d: bc(2,2) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          y = prob_hi(2)
          do k=lo(3),hi(3)+1
             z = prob_lo(3) + dble(k)*dx(3)
             do i=lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                s(i,hi(2)+1:hi(2)+ng_s,k) = s(i,hi(2),k) + dx(2)*inhomogeneous_bc_val_3d(bccomp,x,y,z)
             end do
          end do
       end if
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_inhomogeneous_3d: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,1) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          z = prob_lo(3)
          do j=lo(2),hi(2)
             y = prob_lo(3) + (dble(j)+0.5d0)*dx(2)
             do i=lo(1),hi(1)+1
                x = prob_lo(1) + dble(i)*dx(1)
                s(i,j,lo(3)-ng_s:lo(3)-1) = 2.d0*inhomogeneous_bc_val_3d(bccomp,x,y,z) - s(i,j,lo(3))
             end do
          end do
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          z = prob_lo(3)
          do j=lo(2),hi(2)+1
             y = prob_lo(3) + dble(j)*dx(2)
             do i=lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                s(i,j,lo(3)-ng_s:lo(3)-1) = 2.d0*inhomogeneous_bc_val_3d(bccomp,x,y,z) - s(i,j,lo(3))
             end do
          end do
       else if (bccomp .eq. 3) then 
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel_inhomogeneous
       end if
    else if (bc(3,1) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          z = prob_lo(3)
          do j=lo(2),hi(2)
             y = prob_lo(3) + (dble(j)+0.5d0)*dx(2)
             do i=lo(1),hi(1)+1
                x = prob_lo(1) + dble(i)*dx(1)
                s(i,j,lo(3)-ng_s:lo(3)-1) = s(i,j,lo(3)) - dx(3)*inhomogeneous_bc_val_3d(bccomp,x,y,z)
             end do
          end do
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          z = prob_lo(3)
          do j=lo(2),hi(2)+1
             y = prob_lo(3) + dble(j)*dx(2)
             do i=lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                s(i,j,lo(3)-ng_s:lo(3)-1) = s(i,j,lo(3)) - dx(3)*inhomogeneous_bc_val_3d(bccomp,x,y,z)
             end do
          end do
       else if (bccomp .eq. 3) then 
          print *,'physbc_macvel_inhomogeneous_3d: bc(3,1) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(3,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_inhomogeneous_3d: bc(3,1) =',bc(3,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,2) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          z = prob_hi(3)
          do j=lo(2),hi(2)
             y = prob_lo(3) + (dble(j)+0.5d0)*dx(2)
             do i=lo(1),hi(1)+1
                x = prob_lo(1) + dble(i)*dx(1)
                s(i,j,hi(3)+1:hi(3)+ng_s) = 2.d0*inhomogeneous_bc_val_3d(bccomp,x,y,z) - s(i,j,hi(3))
             end do
          end do
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          z = prob_hi(3)
          do j=lo(2),hi(2)+1
             y = prob_lo(3) + dble(j)*dx(2)
             do i=lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                s(i,j,hi(3)+1:hi(3)+ng_s) = 2.d0*inhomogeneous_bc_val_3d(bccomp,x,y,z) - s(i,j,hi(3))
             end do
          end do
       else if (bccomp .eq. 3) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel_inhomogeneous
       end if
    else if (bc(3,2) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          z = prob_hi(3)
          do j=lo(2),hi(2)
             y = prob_lo(3) + (dble(j)+0.5d0)*dx(2)
             do i=lo(1),hi(1)+1
                x = prob_lo(1) + dble(i)*dx(1)
                s(i,j,hi(3)+1:hi(3)+ng_s) = s(i,j,hi(3)) + dx(3)*inhomogeneous_bc_val_3d(bccomp,x,y,z)
             end do
          end do
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          z = prob_hi(3)
          do j=lo(2),hi(2)+1
             y = prob_lo(3) + dble(j)*dx(2)
             do i=lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                s(i,j,hi(3)+1:hi(3)+ng_s) = s(i,j,hi(3)) + dx(3)*inhomogeneous_bc_val_3d(bccomp,x,y,z)
             end do
          end do
       else if (bccomp .eq. 3) then
          print *,'physbc_macvel_inhomogeneous_3d: bc(3,2) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(3,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_inhomogeneous_3d: bc(3,2) =',bc(3,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine physbc_macvel_inhomogeneous_3d

  subroutine multifab_physbc_domainvel_inhomogeneous(s,scomp,bccomp,num_comp,the_bc_level,dx)

    type(multifab) , intent(inout) :: s
    integer        , intent(in   ) :: scomp,bccomp,num_comp
    type(bc_level) , intent(in   ) :: the_bc_level
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! Local
    integer                  :: lo(get_dim(s)),hi(get_dim(s))
    integer                  :: i,ng_s,dm
    real(kind=dp_t), pointer :: sp(:,:,:,:)

    if (scomp .ne. 1) then
       call bl_error('multifab_physbc_domainvel_inhomogeneous expects scomp = 1')
    end if

    if (bccomp .gt. get_dim(s)) then
       call bl_error('multifab_physbc_domainvel_inhomogeneous expects bccomp <= dm')
    end if

    if (num_comp .ne. 1) then
       call bl_error('multifab_physbc_domainvel_inhomogeneous expects num_comp = 1')
    end if

    dm = get_dim(s)
    ng_s = nghost(s)

    do i=1,nfabs(s)
       sp => dataptr(s,i)
       lo = lwb(get_box(s,i))
       hi = upb(get_box(s,i))
       select case (dm)
       case (2)
          call physbc_domainvel_inhomogeneous_2d(sp(:,:,1,scomp), ng_s, &
                                                 lo, hi, &
                                                 the_bc_level%adv_bc_level_array(i,:,:,bccomp), &
                                                 bccomp,dx)
       case (3)
          call physbc_domainvel_inhomogeneous_3d(sp(:,:,:,scomp), ng_s, lo, hi, &
                                                 the_bc_level%adv_bc_level_array(i,:,:,bccomp), &
                                                 bccomp,dx)
       end select
    end do
 
  end subroutine multifab_physbc_domainvel_inhomogeneous

  subroutine physbc_domainvel_inhomogeneous_2d(s,ng_s,lo,hi,bc,bccomp,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    real(kind=dp_t), intent(inout) ::            s(lo(1)-ng_s:,lo(2)-ng_s:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! local
    integer :: i,j
    real(kind=dp_t) :: x,y

    if (bccomp .ne. 1 .and. bccomp .ne. 2) then
       call bl_error('physbc_domainvel_inhomogeneous_2d requires bccomp = 1 or 2')
    end if

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    ! staggered x-velocity
    if (bccomp .eq. 1) then
       if (bc(1,1) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          x = prob_lo(1)
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             s(lo(1),j) = inhomogeneous_bc_val_2d(bccomp,x,y)
          end do
       else if (bc(1,1) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_inhomogeneous_2d: bc(1,1) =',bc(1,1),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    ! staggered x-velocity
    if (bccomp .eq. 1) then
       if (bc(1,2) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          x = prob_hi(1)
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             s(hi(1)+1,j) = inhomogeneous_bc_val_2d(bccomp,x,y)
          end do
       else if (bc(1,2) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_inhomogeneous_2d: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    ! staggered y-velocity
    if (bccomp .eq. 2) then
       if (bc(2,1) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          y = prob_lo(2)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             s(i,lo(2)) = inhomogeneous_bc_val_2d(bccomp,x,y)
          end do
       else if (bc(2,1) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_inhomogeneous_2d: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    ! staggered y-velocity
    if (bccomp .eq. 2) then
       if (bc(2,2) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          y = prob_hi(2)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             s(i,hi(2)+1) = inhomogeneous_bc_val_2d(bccomp,x,y)
          end do
       else if (bc(2,2) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_inhomogeneous_2d: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    end if

  end subroutine physbc_domainvel_inhomogeneous_2d

  subroutine physbc_domainvel_inhomogeneous_3d(s,ng_s,lo,hi,bc,bccomp,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! local
    integer :: i,j,k
    real(kind=dp_t) :: x,y,z

    if (bccomp .ne. 1 .and. bccomp .ne. 2 .and. bccomp .ne. 3) then
       call bl_error('physbc_domainvel_inhomogeneous_3d requires bccomp = 1, 2, or 3')
    end if

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    ! staggered x-velocity
    if (bccomp .eq. 1) then
       if (bc(1,1) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          x = prob_lo(1)
          do k=lo(3),hi(3)
             z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
             do j=lo(2),hi(2)
                y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
                s(lo(1),j,k) = inhomogeneous_bc_val_3d(bccomp,x,y,z)
             end do
          end do
       else if (bc(1,1) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_inhomogeneous_3d: bc(1,1) =',bc(1,1),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    ! staggered x-velocity
    if (bccomp .eq. 1) then
       if (bc(1,2) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          x = prob_hi(1)
          do k=lo(3),hi(3)
             z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
             do j=lo(2),hi(2)
                y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
                s(hi(1)+1,j,k) = inhomogeneous_bc_val_3d(bccomp,x,y,z)
             end do
          end do
       else if (bc(1,2) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_inhomogeneous_3d: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    ! staggered y-velocity
    if (bccomp .eq. 2) then
       if (bc(2,1) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          y = prob_lo(2)
          do k=lo(3),hi(3)
             z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
             do i=lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                s(i,lo(2),k) = inhomogeneous_bc_val_3d(bccomp,x,y,z)
             end do
          end do
       else if (bc(2,1) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_inhomogeneous_3d: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    ! staggered y-velocity
    if (bccomp .eq. 2) then
       if (bc(2,2) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          y = prob_hi(2)
          do k=lo(3),hi(3)
             z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
             do i=lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                s(i,hi(2)+1,k) = inhomogeneous_bc_val_3d(bccomp,x,y,z)
             end do
          end do
       else if (bc(2,2) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_inhomogeneous_3d: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    end if

!!!!!!!!!!!!!!!!!!
! lo-z boundary
!!!!!!!!!!!!!!!!!!

    ! staggered z-velocity
    if (bccomp .eq. 3) then
       if (bc(3,1) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          z = prob_lo(3)
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             do i=lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                s(i,j,lo(3)) = inhomogeneous_bc_val_3d(bccomp,x,y,z)
             end do
          end do
       else if (bc(3,1) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_inhomogeneous_3d: bc(3,1) =',bc(3,1),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
      end if
   end if

!!!!!!!!!!!!!!!!!!
! hi-z boundary
!!!!!!!!!!!!!!!!!!

   ! staggered z-velocity
   if (bccomp .eq. 3) then
      if (bc(3,2) .eq. DIR_VEL) then
         ! set domain face value to Dirichlet value
          z = prob_hi(3)
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             do i=lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                s(i,j,hi(3)+1) = inhomogeneous_bc_val_3d(bccomp,x,y,z)
             end do
          end do
      else if (bc(3,2) .eq. INTERIOR) then
         ! either periodic or interior; do nothing
      else
         print *,'physbc_domainvel_inhomogeneous_3d: bc(3,2) =',bc(3,2),' for bccomp =',bccomp
         call bl_error('NOT SUPPORTED')
      end if
   end if

 end subroutine physbc_domainvel_inhomogeneous_3d
 
 function inhomogeneous_bc_val_2d(comp,x,y) result(val)

   use bl_constants_module, only: M_PI

   integer        , intent(in   ) :: comp
   real(kind=dp_t), intent(in   ) :: x,y
   real(kind=dp_t)                :: val

   ! local
   real(kind=dp_t) :: velx,vely,Lx,Ly,time, visc_coef

   visc_coef=coeff_mag(2)/coeff_mag(1)

   select case (prob_sol)
   case (1,2)

      ! Taylor vortex
      velx = 1.d0
      vely = 1.d0
      Lx = prob_hi(1)-prob_lo(1)
      Ly = prob_hi(2)-prob_lo(2)
      time = 0.d0

      if (comp .eq. 1) then
         val = velx - 2.d0*exp(-8.d0*M_PI**2*visc_coef*time/(Lx*Ly)) &
              *cos(2.d0*M_PI*(x-velx*time)/Lx)*sin(2.d0*M_PI*(y-vely*time)/Ly)
      else if (comp .eq. 2) then
         val = vely + 2.d0*exp(-8.d0*M_PI**2*visc_coef*time/(Lx*Ly)) &
              *sin(2.d0*M_PI*(x-velx*time)/Lx)*cos(2.d0*M_PI*(y-vely*time)/Ly)
      end if
   
   case (3)
      if (comp .eq. 1) then
         val = exp(-8.0d0*M_PI**2*time)*sin(2.d0*M_PI*x)*sin(2.d0*M_PI*y)
      else if (comp .eq. 2) then
         val = exp(-8.0d0*M_PI**2*time)*cos(2.d0*M_PI*x)*cos(2.d0*M_PI*y)
      end if

   case (20)

      ! lid driven cavity moving hi-y wall with +x velocity
      if (y .eq. prob_hi(2) .and. comp .eq. 1) then
         val = sin(M_PI*(x-prob_lo(1))/(prob_hi(1)-prob_lo(1)))
      else
         val = 0.d0
      end if
  
   case (31,32) ! div not free, for testing accuracy

      if (comp .eq. 1) then
         val = sin(2.d0*M_PI*x)*sin(2.d0*M_PI*y)
      elseif (comp .eq. 2) then
         val = cos(2.d0*M_PI*x)*sin(2.d0*M_PI*y)
      else 
         call bl_error('inhomogeneous_bc_val_3d comp should not be greater than  3')
      end if

   case default

      val = 0.d0

   end select

 end function inhomogeneous_bc_val_2d

 function inhomogeneous_bc_val_3d(comp,x,y,z) result(val)

   use bl_constants_module, only: M_PI

   integer        , intent(in   ) :: comp
   real(kind=dp_t), intent(in   ) :: x,y,z
   real(kind=dp_t)                :: val

   ! local
   real(kind=dp_t) :: velx,vely,Lx,Ly,time, visc_coef
   real(kind=dp_t) :: ufac,vfac,pfac,hx,hy,hz,freq

   visc_coef=coeff_mag(2)/coeff_mag(1)

   select case (prob_sol)
   case (1)

      ! Taylor vortex, quasi-2D
      velx = 1.d0
      vely = 1.d0
      Lx = prob_hi(1)-prob_lo(1)
      Ly = prob_hi(2)-prob_lo(2)
      time = 0.d0

      if (comp .eq. 1) then
         val = velx - 2.d0*exp(-8.d0*M_PI**2*visc_coef*time/(Lx*Ly)) &
              *cos(2.d0*M_PI*(x-velx*time)/Lx)*sin(2.d0*M_PI*(y-vely*time)/Ly)
      else if (comp .eq. 2) then
         val = vely + 2.d0*exp(-8.d0*M_PI**2*visc_coef*time/(Lx*Ly)) &
              *sin(2.d0*M_PI*(x-velx*time)/Lx)*cos(2.d0*M_PI*(y-vely*time)/Ly)
      else
         val = 0.d0
      end if

   case (2) 

      !ABC flow    
      freq  = 2.d0*M_PI

      ufac = dexp(-freq*freq*visc_coef*time)
      pfac = dexp(-2.0d0*freq*freq*visc_coef*time)

      hx = ABC_coefs(1)
      hy = ABC_coefs(2)
      hz = ABC_coefs(3)

      if (comp .eq. 1) then
         val = 1.0d0 + ufac*(hz*cos(freq*(y-time))+hx*sin(freq*(z-time)))
      else if (comp .eq. 2) then
         val = 1.0d0 + ufac*(hy*sin(freq*(x-time))+hx*cos(freq*(z-time)))
      else if (comp .eq. 3) then
         val = 1.0d0 + ufac*(hy*cos(freq*(x-time))+hz*sin(freq*(y-time)))
      else 
         call bl_error('inhomogeneous_bc_val_3d comp should not be greater than  3')
      end if

   case (20)

      ! lid driven cavity moving hi-z wall with +x and +y velocity
      if (z .eq. prob_hi(3) .and. (comp .eq. 1 .or. comp .eq. 2) ) then
         val = sin(M_PI*(x-prob_lo(1))/(prob_hi(1)-prob_lo(1))) &
              *sin(M_PI*(y-prob_lo(2))/(prob_hi(2)-prob_lo(2)))
      else
         val = 0.d0
      end if

   case (31,32) ! div not free, for testing accuracy

      if (comp .eq. 1) then
         val = sin(2.d0*M_PI*x)*sin(2.d0*M_PI*y)*sin(2.d0*M_PI*z)
      else if (comp .eq. 2) then
         val = cos(2.d0*M_PI*x)*sin(2.d0*M_PI*y)*sin(2.d0*M_PI*z)
      else if (comp .eq. 3) then
         val = sin(2.d0*M_PI*x)*sin(2.d0*M_PI*y)*cos(2.d0*M_PI*z)
      else 
         call bl_error('inhomogeneous_bc_val_3d comp should not be greater than  3')
      end if

   case default

      val = 0.d0

   end select

 end function inhomogeneous_bc_val_3d

end module multifab_physbc_inhomogeneous_module
