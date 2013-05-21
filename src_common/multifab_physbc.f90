module multifab_physbc_module

  use multifab_module
  use define_bc_module
  use bc_module
  use bl_error_module
  use inhomogeneous_bc_val_module
  use probin_common_module, only: prob_lo, prob_hi

  implicit none

  private

  public :: multifab_physbc, multifab_physbc_macvel, multifab_physbc_domainvel, &
       multifab_physbc_macvel_inhomogeneous, multifab_physbc_domainvel_inhomogeneous

contains

  subroutine multifab_physbc(s,scomp,bccomp,num_comp,the_bc_level)

    ! this fills ghost cells for rho and pressure/phi.
    ! as well as for transport coefficients (alpha/beta/gamma)

    type(multifab) , intent(inout) :: s
    integer        , intent(in   ) :: scomp,bccomp,num_comp
    type(bc_level) , intent(in   ) :: the_bc_level
   
    ! Local
    integer                  :: lo(get_dim(s)),hi(get_dim(s))
    integer                  :: i,ng,dm
    real(kind=dp_t), pointer :: sp(:,:,:,:)

    if (num_comp .ne. 1) then
       call bl_error('multifab_physbc expects num_comp = 1')
    end if

    ng = nghost(s)
    dm = get_dim(s)
    
    do i=1,nfabs(s)
       sp => dataptr(s,i)
       lo = lwb(get_box(s,i))
       hi = upb(get_box(s,i))
       select case (dm)
       case (2)
          call physbc_2d(sp(:,:,1,scomp), lo, hi, ng, &
                         the_bc_level%adv_bc_level_array(i,:,:,bccomp),bccomp)
       case (3)
          call physbc_3d(sp(:,:,:,scomp), lo, hi, ng, &
                         the_bc_level%adv_bc_level_array(i,:,:,bccomp),bccomp)
       end select
    end do
 
  end subroutine multifab_physbc

  subroutine physbc_2d(s,lo,hi,ng,bc,bccomp)

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) ::    s(lo(1)-ng:,lo(2)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp

    ! Local variables
    integer :: i,j

    if (bccomp .ne. 3 .and. bccomp .ne. 4) then
       call bl_error('physbc_2d requires bccomp = 3 (pressure)')
    end if

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. FOEXTRAP) then
       ! copy interior value
       do j=lo(2)-ng,hi(2)+ng
          s(lo(1)-ng:lo(1)-1,j) = s(lo(1),j)
       end do
    else if (bc(1,1) .eq. EXT_DIR) then
       call bl_error('physbc_2d: need to write Dirichlet primitive condition for bc(1,1)')
    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_2d: bc(1,1) =',bc(1,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. FOEXTRAP) then
       ! copy interior value
       do j=lo(2)-ng,hi(2)+ng
          s(hi(1)+1:hi(1)+ng,j) = s(hi(1),j)
       end do
    else if (bc(1,2) .eq. EXT_DIR) then
       call bl_error('physbc_2d: need to write Dirichlet primitive condition for bc(1,2)')
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_2d: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. FOEXTRAP) then
       ! copy interior value
       do i=lo(1)-ng,hi(1)+ng
          s(i,lo(2)-ng:lo(2)-1) = s(i,lo(2))
       end do
    else if (bc(2,1) .eq. EXT_DIR) then
       call bl_error('physbc_2d: need to write Dirichlet primitive condition for bc(2,1)')
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_2d: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. FOEXTRAP) then
       ! copy interior value
       do i=lo(1)-ng,hi(1)+ng
          s(i,hi(2)+1:hi(2)+ng) = s(i,hi(2))
       end do
    else if (bc(2,2) .eq. EXT_DIR) then
       call bl_error('physbc_2d: need to write Dirichlet primitive condition for bc(2,2)')
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_2d: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine physbc_2d

  subroutine physbc_3d(s,lo,hi,ng,bc,bccomp)

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) ::    s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp

    ! Local variables
    integer :: i,j,k

    if (bccomp .ne. 4 .and. bccomp .ne. 5) then
       call bl_error('physbc_3d requires bccomp = 4 (pressure)')
    end if

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. FOEXTRAP) then
       ! copy interior value
       do k=lo(3)-ng,hi(3)+ng
          do j=lo(2)-ng,hi(2)+ng
             s(lo(1)-ng:lo(1)-1,j,k) = s(lo(1),j,k)
          end do
       end do
    else if (bc(1,1) .eq. EXT_DIR) then
       call bl_error('physbc_3d: need to write Dirichlet primitive condition for bc(1,1)')
    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_3d: bc(1,1) =',bc(1,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. FOEXTRAP) then
       ! copy interior value
       do k=lo(3)-ng,hi(3)+ng
          do j=lo(2)-ng,hi(2)+ng
             s(hi(1)+1:hi(1)+ng,j,k) = s(hi(1),j,k)
          end do
       end do
    else if (bc(1,2) .eq. EXT_DIR) then
       call bl_error('physbc_3d: need to write Dirichlet primitive condition for bc(1,2)')
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_3d: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. FOEXTRAP) then
       ! copy interior value
       do k=lo(3)-ng,hi(3)+ng
          do i=lo(1)-ng,hi(1)+ng
             s(i,lo(2)-ng:lo(2)-1,k) = s(i,lo(2),k)
          end do
       end do
    else if (bc(2,1) .eq. EXT_DIR) then
       call bl_error('physbc_3d: need to write Dirichlet primitive condition for bc(2,1)')
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_3d: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. FOEXTRAP) then
       ! copy interior value
       do k=lo(3)-ng,hi(3)+ng
          do i=lo(1)-ng,hi(1)+ng
             s(i,hi(2)+1:hi(2)+ng,k) = s(i,hi(2),k)
          end do
       end do
    else if (bc(2,2) .eq. EXT_DIR) then
       call bl_error('physbc_3d: need to write Dirichlet primitive condition for bc(2,2)')
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_3d: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,1) .eq. FOEXTRAP) then
       ! copy interior value
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             s(i,j,lo(3)-ng:lo(3)-1) = s(i,j,lo(3))
          end do
       end do
    else if (bc(3,1) .eq. EXT_DIR) then
       call bl_error('physbc_3d: need to write Dirichlet primitive condition for bc(3,1)')
    else if (bc(3,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_3d: bc(3,1) =',bc(3,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,2) .eq. FOEXTRAP) then
       ! copy interior value
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             s(i,j,hi(3)+1:hi(3)+ng) = s(i,j,hi(3))
          end do
       end do
    else if (bc(3,2) .eq. EXT_DIR) then
       call bl_error('physbc_3d: need to write Dirichlet primitive condition for bc(3,2)')
    else if (bc(3,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_3d: bc(3,2) =',bc(3,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine physbc_3d

  subroutine multifab_physbc_macvel(s,scomp,bccomp,num_comp,the_bc_level,dx)

    type(multifab) , intent(inout) :: s
    integer        , intent(in   ) :: scomp,bccomp,num_comp
    type(bc_level) , intent(in   ) :: the_bc_level
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! Local
    integer                  :: lo(get_dim(s)),hi(get_dim(s))
    integer                  :: i,ng_s,dm
    real(kind=dp_t), pointer :: sp(:,:,:,:)

    if (scomp .ne. 1) then
       call bl_error('multifab_physbc_macvel expects scomp = 1')
    end if

    if (bccomp .gt. get_dim(s)) then
       call bl_error('multifab_physbc_macvel expects bccomp <= dm')
    end if

    if (num_comp .ne. 1) then
       call bl_error('multifab_physbc_macvel expects num_comp = 1')
    end if

    ng_s = nghost(s)
    dm = get_dim(s)
    
    do i=1,nfabs(s)
       sp => dataptr(s,i)
       lo = lwb(get_box(s,i))
       hi = upb(get_box(s,i))
       select case (dm)
       case (2)
          call physbc_macvel_2d(sp(:,:,1,scomp), ng_s, lo, hi, &
                                the_bc_level%adv_bc_level_array(i,:,:,bccomp), &
                                bccomp,dx)
       case (3)
          call physbc_macvel_3d(sp(:,:,:,scomp), ng_s, lo, hi, &
                                the_bc_level%adv_bc_level_array(i,:,:,bccomp), &
                                bccomp,dx)
       end select
    end do
 
  end subroutine multifab_physbc_macvel

  subroutine physbc_macvel_2d(s,ng_s,lo,hi,bc,bccomp,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng_s:,lo(2)-ng_s:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer :: i,j

    if (bccomp .ne. 1 .and. bccomp .ne. 2) then
       call bl_error('physbc_macvel_2d requires bccomp = 1 or 2')
    end if

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do j=lo(2),hi(2)+1
             s(lo(1)-ng_s:lo(1)-1,j) = -s(lo(1),j)
          end do
       end if
    else if (bc(1,1) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          print *,'physbc_macvel_2d: bc(1,1) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do j=lo(2),hi(2)+1
             s(lo(1)-ng_s:lo(1)-1,j) = s(lo(1),j)
          end do
       end if
    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_2d: bc(1,1) =',bc(1,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do j=lo(2),hi(2)+1
             s(hi(1)+1:hi(1)+ng_s,j) = -s(hi(1),j)
          end do
       end if
    else if (bc(1,2) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          print *,'physbc_macvel_2d: bc(1,2) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do j=lo(2),hi(2)+1
             s(hi(1)+1:hi(1)+ng_s,j) = s(hi(1),j)
          end do
       end if
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_2d: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do i=lo(1),hi(1)+1
             s(i,lo(2)-ng_s:lo(2)-1) = -s(i,lo(2))
          end do
       else if (bccomp .eq. 2) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel
       end if
    else if (bc(2,1) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do i=lo(1),hi(1)+1
             s(i,lo(2)-ng_s:lo(2)-1) = s(i,lo(2))
          end do
       else if (bccomp .eq. 2) then
          print *,'physbc_macvel_2d: bc(2,1) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_2d: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do i=lo(1),hi(1)+1
             s(i,hi(2)+1:hi(2)+ng_s) = -s(i,hi(2))
          end do
       else if (bccomp .eq. 2) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel
       end if
    else if (bc(2,2) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do i=lo(1),hi(1)+1
             s(i,hi(2)+1:hi(2)+ng_s) = s(i,hi(2))
          end do
       else if (bccomp .eq. 2) then
          print *,'physbc_macvel_2d: bc(2,2) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_2d: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine physbc_macvel_2d

  subroutine physbc_macvel_3d(s,ng_s,lo,hi,bc,bccomp,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer :: i,j,k

    if (bccomp .ne. 1 .and. bccomp .ne. 2 .and. bccomp .ne. 3) then
       call bl_error('physbc_macvel_3d requires bccomp = 1, 2 or 3')
    end if

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do k=lo(3)-ng_s,hi(3)+ng_s
             do j=lo(2),hi(2)+1
                s(lo(1)-ng_s:lo(1)-1,j,k) = -s(lo(1),j,k)
             end do
          end do
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do k=lo(3),hi(3)+1
             do j=lo(2)-ng_s,hi(2)+ng_s
                s(lo(1)-ng_s:lo(1)-1,j,k) = -s(lo(1),j,k)
             end do
          end do
       end if
    else if (bc(1,1) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          print *,'physbc_macvel_3d: bc(1,1) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do k=lo(3)-ng_s,hi(3)+ng_s
             do j=lo(2),hi(2)+1
                s(lo(1)-ng_s:lo(1)-1,j,k) = s(lo(1),j,k)
             end do
          end do
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do k=lo(3),hi(3)+1
             do j=lo(2)-ng_s,hi(2)+ng_s
                s(lo(1)-ng_s:lo(1)-1,j,k) = s(lo(1),j,k)
             end do
          end do
       end if
    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_3d: bc(1,1) =',bc(1,1),'for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do k=lo(3)-ng_s,hi(3)+ng_s
             do j=lo(2),hi(2)+1
                s(hi(1)+1:hi(1)+ng_s,j,k) = -s(hi(1),j,k)
             end do
          end do
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do k=lo(3),hi(3)+1
             do j=lo(2)-ng_s,hi(2)+ng_s
                s(hi(1)+1:hi(1)+ng_s,j,k) = -s(hi(1),j,k)
             end do
          end do
       end if
    else if (bc(1,2) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          print *,'physbc_macvel_3d: bc(1,2) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do k=lo(3)-ng_s,hi(3)+ng_s
             do j=lo(2),hi(2)+1
                s(hi(1)+1:hi(1)+ng_s,j,k) = s(hi(1),j,k)
             end do
          end do
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do k=lo(3),hi(3)+1
             do j=lo(2)-ng_s,hi(2)+ng_s
                s(hi(1)+1:hi(1)+ng_s,j,k) = s(hi(1),j,k)
             end do
          end do
       end if
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_3d: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do k=lo(3)-ng_s,hi(3)+ng_s
             do i=lo(1),hi(1)+1
                s(i,lo(2)-ng_s:lo(2)-1,k) = -s(i,lo(2),k)
             end do
          end do
       else if (bccomp .eq. 2) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do k=lo(3),hi(3)+1
             do i=lo(1)-ng_s,hi(1)+ng_s
                s(i,lo(2)-ng_s:lo(2)-1,k) = -s(i,lo(2),k)
             end do
          end do
       end if
    else if (bc(2,1) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do k=lo(3)-ng_s,hi(3)+ng_s
             do i=lo(1),hi(1)+1
                s(i,lo(2)-ng_s:lo(2)-1,k) = s(i,lo(2),k)
             end do
          end do
       else if (bccomp .eq. 2) then
          print *,'physbc_macvel_3d: bc(2,1) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do k=lo(3),hi(3)+1
             do i=lo(1)-ng_s,hi(1)+ng_s
                s(i,lo(2)-ng_s:lo(2)-1,k) = s(i,lo(2),k)
             end do
          end do
       end if
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_3d: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do k=lo(3)-ng_s,hi(3)+ng_s
             do i=lo(1),hi(1)+1
                s(i,hi(2)+1:hi(2)+ng_s,k) = -s(i,hi(2),k)
             end do
          end do
       else if (bccomp .eq. 2) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do k=lo(3),hi(3)+1
             do i=lo(1)-ng_s,hi(1)+ng_s
                s(i,hi(2)+1:hi(2)+ng_s,k) = -s(i,hi(2),k)
             end do
          end do
       end if
    else if (bc(2,2) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do k=lo(3)-ng_s,hi(3)+ng_s
             do i=lo(1),hi(1)+1
                s(i,hi(2)+1:hi(2)+ng_s,k) = s(i,hi(2),k)
             end do
          end do
       else if (bccomp .eq. 2) then
          print *,'physbc_macvel_3d: bc(2,2) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do k=lo(3),hi(3)+1
             do i=lo(1)-ng_s,hi(1)+ng_s
                s(i,hi(2)+1:hi(2)+ng_s,k) = s(i,hi(2),k)
             end do
          end do
       end if
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_3d: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,1) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do i=lo(1),hi(1)+1
             do j=lo(2)-ng_s,hi(2)+ng_s
                s(i,j,lo(3)-ng_s:lo(3)-1) = -s(i,j,lo(3))
             end do
          end do
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do i=lo(1)-ng_s,hi(1)+ng_s
             do j=lo(2),hi(2)+1
                s(i,j,lo(3)-ng_s:lo(3)-1) = -s(i,j,lo(3))
             end do
          end do
       else if (bccomp .eq. 3) then 
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel
       end if
    else if (bc(3,1) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do i=lo(1),hi(1)+1
             do j=lo(2)-ng_s,hi(2)+ng_s
                s(i,j,lo(3)-ng_s:lo(3)-1) = s(i,j,lo(3))
             end do
          end do
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do i=lo(1)-ng_s,hi(1)+ng_s
             do j=lo(2),hi(2)+1
                s(i,j,lo(3)-ng_s:lo(3)-1) = s(i,j,lo(3))
             end do
          end do
       else if (bccomp .eq. 3) then 
          print *,'physbc_macvel_3d: bc(3,1) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(3,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_3d: bc(3,1) =',bc(3,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,2) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do i=lo(1),hi(1)+1
             do j=lo(2)-ng_s,hi(2)+ng_s
                s(i,j,hi(3)+1:hi(3)+ng_s) = -s(i,j,hi(3))
             end do
          end do
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do i=lo(1)-ng_s,hi(1)+ng_s
             do j=lo(2),hi(2)+1
                s(i,j,hi(3)+1:hi(3)+ng_s) = -s(i,j,hi(3))
             end do
          end do
       else if (bccomp .eq. 3) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel
       end if
    else if (bc(3,2) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do i=lo(1),hi(1)+1
             do j=lo(2)-ng_s,hi(2)+ng_s
                s(i,j,hi(3)+1:hi(3)+ng_s) = s(i,j,hi(3))
             end do
          end do
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do i=lo(1)-ng_s,hi(1)+ng_s
             do j=lo(2),hi(2)+1
                s(i,j,hi(3)+1:hi(3)+ng_s) = s(i,j,hi(3))
             end do
          end do
       else if (bccomp .eq. 3) then
          print *,'physbc_macvel_3d: bc(3,2) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(3,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_3d: bc(3,2) =',bc(3,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine physbc_macvel_3d

  subroutine multifab_physbc_domainvel(s,scomp,bccomp,num_comp,the_bc_level,dx)

    type(multifab) , intent(inout) :: s
    integer        , intent(in   ) :: scomp,bccomp,num_comp
    type(bc_level) , intent(in   ) :: the_bc_level
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! Local
    integer                  :: lo(get_dim(s)),hi(get_dim(s))
    integer                  :: i,ng_s,dm
    real(kind=dp_t), pointer :: sp(:,:,:,:)

    if (scomp .ne. 1) then
       call bl_error('multifab_physbc_domainvel expects scomp = 1')
    end if

    if (bccomp .gt. get_dim(s)) then
       call bl_error('multifab_physbc_domainvel expects bccomp <= dm')
    end if

    if (num_comp .ne. 1) then
       call bl_error('multifab_physbc_domainvel expects num_comp = 1')
    end if

    ng_s = nghost(s)
    dm = get_dim(s)
    
    do i=1,nfabs(s)
       sp => dataptr(s,i)
       lo = lwb(get_box(s,i))
       hi = upb(get_box(s,i))
       select case (dm)
       case (2)
          call physbc_domainvel_2d(sp(:,:,1,scomp), ng_s, lo, hi, &
                                   the_bc_level%adv_bc_level_array(i,:,:,bccomp), &
                                   bccomp,dx)
       case (3)
          call physbc_domainvel_3d(sp(:,:,:,scomp), ng_s, lo, hi, &
                                   the_bc_level%adv_bc_level_array(i,:,:,bccomp), &
                                   bccomp,dx)
       end select
    end do
 
  end subroutine multifab_physbc_domainvel

  subroutine physbc_domainvel_2d(s,ng_s,lo,hi,bc,bccomp,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng_s:,lo(2)-ng_s:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp
    real(kind=dp_t), intent(in   ) :: dx(:)

    if (bccomp .ne. 1 .and. bccomp .ne. 2) then
       call bl_error('physbc_domainvel_2d requires bccomp = 1 or 2')
    end if

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    ! staggered x-velocity
    if (bccomp .eq. 1) then
       if (bc(1,1) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          s(lo(1),lo(2):hi(2)) = 0.d0
       else if (bc(1,1) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_2d: bc(1,1) =',bc(1,1),' for bccomp =',bccomp
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
          s(hi(1)+1,lo(2):hi(2)) = 0.d0
       else if (bc(1,2) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_2d: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
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
          s(lo(1):hi(1),lo(2)) = 0.d0
       else if (bc(2,1) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_2d: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
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
          s(lo(1):hi(1),hi(2)+1) = 0.d0
       else if (bc(2,2) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_2d: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    end if

  end subroutine physbc_domainvel_2d

  subroutine physbc_domainvel_3d(s,ng_s,lo,hi,bc,bccomp,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp
    real(kind=dp_t), intent(in   ) :: dx(:)

    if (bccomp .ne. 1 .and. bccomp .ne. 2 .and. bccomp .ne. 3) then
       call bl_error('physbc_domainvel_3d requires bccomp = 1, 2, or 3')
    end if

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    ! staggered x-velocity
    if (bccomp .eq. 1) then
       if (bc(1,1) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          s(lo(1),lo(2):hi(2),lo(3):hi(3)) = 0.d0
       else if (bc(1,1) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_3d: bc(1,1) =',bc(1,1),' for bccomp =',bccomp
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
          s(hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = 0.d0
       else if (bc(1,2) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_3d: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
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
          s(lo(1):hi(1),lo(2),lo(3):hi(3)) = 0.d0
       else if (bc(2,1) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_3d: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
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
          s(lo(1):hi(1),hi(2)+1,lo(3):hi(3)) = 0.d0
       else if (bc(2,2) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_3d: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
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
          s(lo(1):hi(1),lo(2):hi(2),lo(3)) = 0.d0
       else if (bc(3,1) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_3d: bc(3,1) =',bc(3,1),' for bccomp =',bccomp
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
         s(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = 0.d0
      else if (bc(3,2) .eq. INTERIOR) then
         ! either periodic or interior; do nothing
      else
         print *,'physbc_domainvel_3d: bc(3,2) =',bc(3,2),' for bccomp =',bccomp
         call bl_error('NOT SUPPORTED')
      end if
   end if

 end subroutine physbc_domainvel_3d

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

end module multifab_physbc_module
