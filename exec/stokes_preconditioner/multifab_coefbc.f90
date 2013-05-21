module multifab_coefbc_module

  use multifab_module
  use define_bc_module
  use bc_module
  use bl_error_module

  implicit none

  private

  public :: multifab_coefbc

contains

  subroutine multifab_coefbc(s,scomp,bccomp,num_comp,the_bc_level)

    ! this fills ghost cells for transport coefficitns (alpha/beta/gamma)

    type(multifab) , intent(inout) :: s
    integer        , intent(in   ) :: scomp,bccomp,num_comp
    type(bc_level) , intent(in   ) :: the_bc_level
   
    ! Local
    integer                  :: lo(get_dim(s)),hi(get_dim(s))
    integer                  :: i,ng,dm
    real(kind=dp_t), pointer :: sp(:,:,:,:)

    if (num_comp .ne. 1) then
       call bl_error('multifab_coefbc expects num_comp = 1')
    end if

    ng = nghost(s)
    dm = get_dim(s)
    
    do i=1,nfabs(s)
       sp => dataptr(s,i)
       lo = lwb(get_box(s,i))
       hi = upb(get_box(s,i))
       select case (dm)
       case (2)
          call coefbc_2d(sp(:,:,1,scomp), lo, hi, ng, &
                         the_bc_level%adv_bc_level_array(i,:,:,bccomp),bccomp)
       case (3)
          call coefbc_3d(sp(:,:,:,scomp), lo, hi, ng, &
                         the_bc_level%adv_bc_level_array(i,:,:,bccomp),bccomp)
       end select
    end do
 
  end subroutine multifab_coefbc

  subroutine coefbc_2d(s,lo,hi,ng,bc,bccomp)

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) ::    s(lo(1)-ng:,lo(2)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp

    ! Local variables
    integer :: i,j

    if (bccomp .ne. 3 .and. bccomp .ne. 4) then
       call bl_error('coefbc_2d requires bccomp = 3 (pressure)')
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
       ! average interior and ghost into ghost
       do j=lo(2)-ng,hi(2)+ng
          s(lo(1)-ng:lo(1)-1,j) = 0.5d0*(s(lo(1),j)+s(lo(1)-1,j))
       end do
    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'coefbc_2d: bc(1,1) =',bc(1,1),' for bccomp =',bccomp
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
       ! average interior and ghost into ghost
       do j=lo(2)-ng,hi(2)+ng
          s(hi(1)+1:hi(1)+ng,j) = 0.5d0*(s(hi(1),j)+s(hi(1)+1,j))
       end do
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'coefbc_2d: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
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
       ! average interior and ghost into ghost
       do i=lo(1)-ng,hi(1)+ng
          s(i,lo(2)-ng:lo(2)-1) = 0.5d0*(s(i,lo(2))+s(i,lo(2)-1))
       end do
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'coefbc_2d: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
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
       ! average interior and ghost into ghost
       do i=lo(1)-ng,hi(1)+ng
          s(i,hi(2)+1:hi(2)+ng) = 0.5d0*(s(i,hi(2))+s(i,hi(2)+1))
       end do
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'coefbc_2d: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine coefbc_2d

  subroutine coefbc_3d(s,lo,hi,ng,bc,bccomp)

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) ::    s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp

    ! Local variables
    integer :: i,j,k

    if (bccomp .ne. 4 .and. bccomp .ne. 5) then
       call bl_error('coefbc_3d requires bccomp = 4 (pressure)')
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
       ! average interior and ghost into ghost
       do k=lo(3)-ng,hi(3)+ng
          do j=lo(2)-ng,hi(2)+ng
             s(lo(1)-ng:lo(1)-1,j,k) = 0.5d0*(s(lo(1),j,k)+s(lo(1)-1,j,k))
          end do
       end do
    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'coefbc_3d: bc(1,1) =',bc(1,1),' for bccomp =',bccomp
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
       ! average interior and ghost into ghost
       do k=lo(3)-ng,hi(3)+ng
          do j=lo(2)-ng,hi(2)+ng
             s(hi(1)+1:hi(1)+ng,j,k) = 0.5d0*(s(hi(1),j,k)+s(hi(1)+1,j,k))
          end do
       end do
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'coefbc_3d: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
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
       ! average interior and ghost into ghost
       do k=lo(3)-ng,hi(3)+ng
          do i=lo(1)-ng,hi(1)+ng
             s(i,lo(2)-ng:lo(2)-1,k) = 0.5d0*(s(i,lo(2),k)+s(i,lo(2)-1,k))
          end do
       end do
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'coefbc_3d: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
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
       ! average interior and ghost into ghost
       do k=lo(3)-ng,hi(3)+ng
          do i=lo(1)-ng,hi(1)+ng
             s(i,hi(2)+1:hi(2)+ng,k) = 0.5d0*(s(i,hi(2),k)+s(i,hi(2)+1,k))
          end do
       end do
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'coefbc_3d: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
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
       ! average interior and ghost into ghost
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             s(i,j,lo(3)-ng:lo(3)-1) = 0.5d0*(s(i,j,lo(3))+s(i,j,lo(3)-1))
          end do
       end do
    else if (bc(3,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'coefbc_3d: bc(3,1) =',bc(3,1),' for bccomp =',bccomp
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
       ! average interior and ghost into ghost
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             s(i,j,hi(3)+1:hi(3)+ng) = 0.5d0*(s(i,j,hi(3))+s(i,j,hi(3)+1))
          end do
       end do
    else if (bc(3,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'coefbc_3d: bc(3,2) =',bc(3,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine coefbc_3d

end module multifab_coefbc_module
