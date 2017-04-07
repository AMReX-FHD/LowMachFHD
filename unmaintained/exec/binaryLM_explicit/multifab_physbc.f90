module multifab_physbc_module

  use multifab_module
  use define_bc_module
  use bl_types
  use bl_error_module
  use probin_module, only: wallspeed_lo, wallspeed_hi, boyce_bc, prob_lo, prob_hi, pmask, &
       rhobar, c_wall

  implicit none

  private

  public :: multifab_physbc, multifab_physbc_macvel, multifab_physbc_domainvel

contains

  subroutine multifab_physbc(s,start_scomp,start_bccomp,num_comp,the_bc_level)

    type(multifab) , intent(inout) :: s
    integer        , intent(in   ) :: start_scomp,start_bccomp,num_comp
    type(bc_level) , intent(in   ) :: the_bc_level
   
    ! Local
    integer                  :: lo(get_dim(s)),hi(get_dim(s))
    integer                  :: i,ng,dm,scomp,bccomp
    real(kind=dp_t), pointer :: sp(:,:,:,:)

    ng = nghost(s)
    dm = get_dim(s)
    
    do i=1,nfabs(s)
       sp => dataptr(s,i)
       lo = lwb(get_box(s,i))
       hi = upb(get_box(s,i))
       select case (dm)
       case (2)
          do scomp = start_scomp,start_scomp+num_comp-1
             bccomp = start_bccomp + scomp - start_scomp
             call physbc_2d(sp(:,:,1,scomp), lo, hi, ng, &
                            the_bc_level%adv_bc_level_array(i,:,:,bccomp),bccomp)
          end do
       case (3)
          do scomp = start_scomp,start_scomp+num_comp-1
             bccomp = start_bccomp + scomp - start_scomp
             call physbc_3d(sp(:,:,:,scomp), lo, hi, ng, &
                            the_bc_level%adv_bc_level_array(i,:,:,bccomp),bccomp)
          end do
       end select
    end do
 
  end subroutine multifab_physbc

  subroutine physbc_2d(s,lo,hi,ng,bc,icomp)

    use bl_constants_module
    use bc_module

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) ::    s(lo(1)-ng:,lo(2)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: icomp

    ! Local variables
    integer :: i,j

    if (icomp .ne. 3 .and. icomp .ne. 4) then
       call bl_error('physbc_2d requires icomp = 3 or 4')
    end if

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. FOEXTRAP) then
       ! cell-centered rho or c
       do j=lo(2)-ng,hi(2)+ng
          s(lo(1)-ng:lo(1)-1,j) = s(lo(1),j)
       end do
    else if (bc(1,1) .eq. EXT_DIR) then
       if (icomp .eq. 3) then
          ! cell-centered rho
          s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng) = &
               1.d0/(c_wall(1,1)/rhobar(1) + (1.d0-c_wall(1,1))/rhobar(2))
       else if (icomp .eq. 4) then
          ! cell-centered c
          s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng) = c_wall(1,1)
       end if
    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else 
       print *,'physbc_2d: bc(1,1) =',bc(1,1),' for icomp =',icomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. FOEXTRAP) then
       ! cell-centered rho or c
       do j=lo(2)-ng,hi(2)+ng
          s(hi(1)+1:hi(1)+ng,j) = s(hi(1),j)
       end do
    else if (bc(1,2) .eq. EXT_DIR) then
       if (icomp .eq. 3) then
          ! cell-centered rho
          s(hi(1)+1:hi(1)+ng,lo(2)-ng:hi(2)+ng) = &
               1.d0/(c_wall(1,2)/rhobar(1) + (1.d0-c_wall(1,2))/rhobar(2))
       else if (icomp .eq. 4) then
          ! cell-centered c
          s(hi(1)+1:hi(1)+ng,lo(2)-ng:hi(2)+ng) = c_wall(1,2)
       end if
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else 
       print *,'physbc_2d: bc(1,2) =',bc(1,2),' for icomp =',icomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. FOEXTRAP) then
       ! cell-centered rho or c
       do i=lo(1)-ng,hi(1)+ng
          s(i,lo(2)-ng:lo(2)-1) = s(i,lo(2))
       end do
    else if (bc(2,1) .eq. EXT_DIR) then
       if (icomp .eq. 3) then
          ! cell-centered rho
          s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = &
               1.d0/(c_wall(2,1)/rhobar(1) + (1.d0-c_wall(2,1))/rhobar(2))
       else if (icomp .eq. 4) then
          ! cell-centered c
          s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = c_wall(2,1)
       end if
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else 
       print *,'physbc_2d: bc(2,1) =',bc(2,1),' for icomp =',icomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. FOEXTRAP) then
       ! cell-centered rho or c
       do i=lo(1)-ng,hi(1)+ng
          s(i,hi(2)+1:hi(2)+ng) = s(i,hi(2))
       end do
    else if (bc(2,2) .eq. EXT_DIR) then
       if (icomp .eq. 3) then
          ! cell-centered rho
          s(lo(1)-ng:hi(1)+ng,hi(2)+1:hi(2)+ng) = &
               1.d0/(c_wall(2,2)/rhobar(1) + (1.d0-c_wall(2,2))/rhobar(2))
       else if (icomp .eq. 4) then
          ! cell-centered c
          s(lo(1)-ng:hi(1)+ng,hi(2)+1:hi(2)+ng) = c_wall(2,2)
       end if
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else 
       print *,'physbc_2d: bc(2,2) =',bc(2,2),' for icomp =',icomp
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine physbc_2d

  subroutine physbc_3d(s,lo,hi,ng,bc,icomp)

    use bl_constants_module
    use bc_module

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) ::    s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: icomp

    ! Local variables
    integer :: i,j,k

    if (icomp .ne. 4 .and. icomp .ne. 5) then
       call bl_error('physbc_3d requires icomp = 4 or 5')
    end if

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. FOEXTRAP) then
       ! cell-centered rho or c
       do k=lo(3)-ng,hi(3)+ng
          do j=lo(2)-ng,hi(2)+ng
             s(lo(1)-ng:lo(1)-1,j,k) = s(lo(1),j,k)
          end do
       end do
    else if (bc(1,1) .eq. EXT_DIR) then
       if (icomp .eq. 4) then
          ! cell-centered rho
          s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = &
               1.d0/(c_wall(1,1)/rhobar(1) + (1.d0-c_wall(1,1))/rhobar(2))
       else if (icomp .eq. 5) then
          ! cell-centered c
          s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = c_wall(1,1)
       end if
    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else 
       print *,'physbc_3d: bc(1,1) =',bc(1,1),' for icomp =',icomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. FOEXTRAP) then
       ! cell-centered rho or c
       do k=lo(3)-ng,hi(3)+ng
          do j=lo(2)-ng,hi(2)+ng
             s(hi(1)+1:hi(1)+ng,j,k) = s(hi(1),j,k)
          end do
       end do
    else if (bc(1,2) .eq. EXT_DIR) then
       if (icomp .eq. 4) then
          ! cell-centered rho
          s(hi(1)+1:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = &
               1.d0/(c_wall(1,2)/rhobar(1) + (1.d0-c_wall(1,2))/rhobar(2))
       else if (icomp .eq. 5) then
          ! cell-centered c
          s(hi(1)+1:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = c_wall(1,2)
       end if
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else 
       print *,'physbc_3d: bc(1,2) =',bc(1,2),' for icomp =',icomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. FOEXTRAP) then
       ! cell-centered rho or c
       do k=lo(3)-ng,hi(3)+ng
          do i=lo(1)-ng,hi(1)+ng
             s(i,lo(2)-ng:lo(2)-1,k) = s(i,lo(2),k)
          end do
       end do
    else if (bc(2,1) .eq. EXT_DIR) then
       if (icomp .eq. 4) then
          ! cell-centered rho
          s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1,lo(3)-ng:hi(3)+ng) = &
               1.d0/(c_wall(2,1)/rhobar(1) + (1.d0-c_wall(2,1))/rhobar(2))
       else if (icomp .eq. 5) then
          ! cell-centered c
          s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1,lo(3)-ng:hi(3)+ng) = c_wall(2,1)
       end if
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else 
       print *,'physbc_3d: bc(2,1) =',bc(2,1),' for icomp =',icomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. FOEXTRAP) then
       ! cell-centered rho or c
       do k=lo(3)-ng,hi(3)+ng
          do i=lo(1)-ng,hi(1)+ng
             s(i,hi(2)+1:hi(2)+ng,k) = s(i,hi(2),k)
          end do
       end do
    else if (bc(2,2) .eq. EXT_DIR) then
       if (icomp .eq. 4) then
          ! cell-centered rho
          s(lo(1)-ng:hi(1)+ng,hi(2)+1:hi(2)+ng,lo(3)-ng:hi(3)+ng) = &
               1.d0/(c_wall(2,2)/rhobar(1) + (1.d0-c_wall(2,2))/rhobar(2))
       else if (icomp .eq. 5) then
          ! cell-centered c
          s(lo(1)-ng:hi(1)+ng,hi(2)+1:hi(2)+ng,lo(3)-ng:hi(3)+ng) = c_wall(2,2)
       end if
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else 
       print *,'physbc_3d: bc(2,2) =',bc(2,2),' for icomp =',icomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,1) .eq. FOEXTRAP) then
       ! cell-centered rho or c
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             s(i,j,lo(3)-ng:lo(3)-1) = s(i,j,lo(3))
          end do
       end do
    else if (bc(3,1) .eq. EXT_DIR) then
       if (icomp .eq. 4) then
          ! cell-centered rho
          s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = &
               1.d0/(c_wall(3,1)/rhobar(1) + (1.d0-c_wall(3,1))/rhobar(2))
       else if (icomp .eq. 5) then
          ! cell-centered c
          s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = c_wall(3,1)
       end if
    else if (bc(3,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else 
       print *,'physbc_3d: bc(3,1) =',bc(3,1),' for icomp =',icomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,2) .eq. FOEXTRAP) then
       ! cell-centered rho or c
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             s(i,j,hi(3)+1:hi(3)+ng) = s(i,j,hi(3))
          end do
       end do
    else if (bc(3,2) .eq. EXT_DIR) then
       if (icomp .eq. 4) then
          ! cell-centered rho
          s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,hi(3)+1:hi(3)+ng) = &
               1.d0/(c_wall(3,2)/rhobar(1) + (1.d0-c_wall(3,2))/rhobar(2))
       else if (icomp .eq. 5) then
          ! cell-centered c
          s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,hi(3)+1:hi(3)+ng) = c_wall(3,2)
       end if
    else if (bc(3,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else 
       print *,'physbc_3d: bc(3,2) =',bc(3,2),' for icomp =',icomp
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine physbc_3d

  subroutine multifab_physbc_macvel(s,start_scomp,start_bccomp,num_comp,the_bc_level,dx)

    type(multifab) , intent(inout) :: s
    integer        , intent(in   ) :: start_scomp,start_bccomp,num_comp
    type(bc_level) , intent(in   ) :: the_bc_level
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! Local
    integer                  :: lo(get_dim(s)),hi(get_dim(s))
    integer                  :: i,ng_s,dm,scomp,bccomp
    real(kind=dp_t), pointer :: sp(:,:,:,:)

    ng_s = nghost(s)
    dm = get_dim(s)
    
    do i=1,nfabs(s)
       sp => dataptr(s,i)
       lo = lwb(get_box(s,i))
       hi = upb(get_box(s,i))
       select case (dm)
       case (2)
          do scomp = start_scomp,start_scomp+num_comp-1
             bccomp = start_bccomp + scomp - start_scomp
             call physbc_macvel_2d(sp(:,:,1,scomp), ng_s, lo, hi, &
                                   the_bc_level%adv_bc_level_array(i,:,:,bccomp),bccomp,dx)
          end do
       case (3)
          do scomp = start_scomp,start_scomp+num_comp-1
             bccomp = start_bccomp + scomp - start_scomp
             call physbc_macvel_3d(sp(:,:,:,scomp), ng_s, lo, hi, &
                                   the_bc_level%adv_bc_level_array(i,:,:,bccomp),bccomp,dx)
          end do
       end select
    end do
 
  end subroutine multifab_physbc_macvel

  subroutine physbc_macvel_2d(s,ng_s,lo,hi,bc,icomp,dx)

    use bl_constants_module
    use bc_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng_s:,lo(2)-ng_s:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: icomp
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer :: i,j
    real(kind=dp_t) :: x,speed

    if (icomp .ne. 1 .and. icomp .ne. 2) then
       call bl_error('physbc_macvel_2d requires icomp = 1 or 2')
    end if

    if (ng_s .lt. 1) then
       ! no need to do anything
       return
    end if

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. REFLECT_ODD) then
       if (icomp .eq. 2) then
          ! staggered y-velocity
          do i=1,ng_s
             s(lo(1)-i,lo(2)-ng_s:hi(2)+1+ng_s) = 2.d0*wallspeed_lo(1,1) &
                  - s(lo(1)+i-1,lo(2)-ng_s:hi(2)+1+ng_s)
          end do
       else
          print *,'physbc_macvel_2d: bc(1,1) = REFLECT_ODD for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(1,1) .eq. REFLECT_EVEN) then
       if (icomp .eq. 2) then
          ! staggered y-velocity
          do i=1,ng_s
             s(lo(1)-i,lo(2)-ng_s:hi(2)+1+ng_s) = s(lo(1)+i-1,lo(2)-ng_s:hi(2)+1+ng_s)
          end do
       else
          print *,'physbc_macvel_2d: bc(1,1) = REFLECT_EVEN for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(1,1) .eq. EXT_DIR) then
       if (icomp .eq. 1) then
          ! staggered x-velocity
          do i=1,ng_s
             s(lo(1)-i,lo(2)-ng_s:hi(2)+ng_s) = s(lo(1),lo(2)-ng_s:hi(2)+ng_s)
          end do
       else
          print *,'physbc_macvel_2d: bc(1,1) = EXT_DIR for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else 
       print *,'physbc_macvel_2d: bc(1,1) =',bc(1,1),' for icomp =',icomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. REFLECT_ODD) then
        if (icomp .eq. 2) then
          ! staggered y-velocity
          do i=1,ng_s
             s(hi(1)+i,lo(2)-ng_s:hi(2)+1+ng_s) = 2.d0*wallspeed_hi(1,1) &
                  - s(hi(1)+1-i,lo(2)-ng_s:hi(2)+1+ng_s)
          end do
       else
          print *,'physbc_macvel_2d: bc(1,2) = REFLECT_ODD for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(1,2) .eq. REFLECT_EVEN) then
       if (icomp .eq. 2) then
          ! staggered y-velocity
          do i=1,ng_s
             s(hi(1)+i,lo(2)-ng_s:hi(2)+1+ng_s) = s(hi(1)+1-i,lo(2)-ng_s:hi(2)+1+ng_s)
          end do
       else
          print *,'physbc_macvel_2d: bc(1,2) = REFLECT_EVEN for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(1,2) .eq. EXT_DIR) then
       if (icomp .eq. 1) then
          ! staggered x-velocity
          do i=1,ng_s
             s(hi(1)+1+i,lo(2)-ng_s:hi(2)+ng_s) = s(hi(1)+1,lo(2)-ng_s:hi(2)+ng_s)
          end do
       else
          print *,'physbc_macvel_2d: bc(1,2) = EXT_DIR for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else 
       print *,'physbc_macvel_2d: bc(1,2) =',bc(1,2),' for icomp =',icomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. REFLECT_ODD) then
       if (icomp .eq. 1) then
          ! staggered x-velocity
          do j=1,ng_s
             s(lo(1)-ng_s:hi(1)+1+ng_s,lo(2)-j) = 2.d0*wallspeed_lo(1,2) &
                  - s(lo(1)-ng_s:hi(1)+1+ng_s,lo(2)+j-1)
          end do
       else 
          print *,'physbc_macvel_2d: bc(2,1) = REFLECT_ODD for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(2,1) .eq. REFLECT_EVEN) then
       if (icomp .eq. 1) then
          ! staggered x-velocity
          do j=1,ng_s
             s(lo(1)-ng_s:hi(1)+1+ng_s,lo(2)-j) = s(lo(1)-ng_s:hi(1)+1+ng_s,lo(2)+j-1)
          end do
       else
          print *,'physbc_macvel_2d: bc(2,1) = REFLECT_EVEN for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(2,1) .eq. EXT_DIR) then
       if (icomp .eq. 2) then
          ! staggered y-velocity
          do j=1,ng_s
             s(lo(1)-ng_s:hi(1)+ng_s,lo(2)-j) = s(lo(1)-ng_s:hi(1)+ng_s,lo(2))
          end do
       else
          print *,'physbc_macvel_2d: bc(2,1) = EXT_DIR for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else 
       print *,'physbc_macvel_2d: bc(2,1) =',bc(2,1),' for icomp =',icomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. REFLECT_ODD) then
       if (icomp .eq. 1) then
          ! staggered x-velocity
          ! only apply boyce velocity tapering_s for non-periodic in x
          if (boyce_bc .and. (.not. pmask(1)) ) then
             do i=lo(1)-ng_s,hi(1)+1+ng_s
                x = prob_lo(1) + dble(i)*dx(1)
                ! This is from Boyce Griffith: Taper wall velocity toward the corner
                speed = wallspeed_hi(1,2)* &
                   0.5d0*(1.d0 + sin(2.d0*M_PI*(x-prob_lo(1))/prob_hi(1) - 0.5d0*M_PI))
                ! original stencil
                do j=1,ng_s
                   s(i,hi(2)+j) = 2.d0*speed - s(i,hi(2)+1-j)
                end do
                ! higher-order stencil
                ! s(i,hi(2)+1) = (8.d0/3.d0)*speed - 2.d0*s(i,hi(2)) + (1.d0/3.d0)*s(i,hi(2)-1)
             end do
          else
             do j=1,ng_s
                s(lo(1)-ng_s:hi(1)+1+ng_s,hi(2)+j) = 2.d0*wallspeed_hi(1,2) &
                     - s(lo(1)-ng_s:hi(1)+1+ng_s,hi(2)+1-j)
             end do
          end if
       else
          print *,'physbc_macvel_2d: bc(2,2) = REFLECT_ODD for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(2,2) .eq. REFLECT_EVEN) then
       if (icomp .eq. 1) then
          ! staggered x-velocity
          do j=1,ng_s
             s(lo(1)-ng_s:hi(1)+1+ng_s,hi(2)+j) = s(lo(1)-ng_s:hi(1)+1+ng_s,hi(2)+1-j)
          end do
       else
          print *,'physbc_macvel_2d: bc(2,2) = REFLECT_EVEN for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(2,2) .eq. EXT_DIR) then
       if (icomp .eq. 2) then
          ! staggered y-velocity
          do j=1,ng_s
             s(lo(1)-ng_s:hi(1)+ng_s,hi(2)+1+j) =  s(lo(1)-ng_s:hi(1)+ng_s,hi(2)+1)
          end do
       else
          print *,'physbc_macvel_2d: bc(2,2) = EXT_DIR for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else 
       print *,'physbc_macvel_2d: bc(2,2) =',bc(2,2),' for icomp =',icomp
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine physbc_macvel_2d

  subroutine physbc_macvel_3d(s,ng_s,lo,hi,bc,icomp,dx)

    use bl_constants_module
    use bc_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: icomp
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer :: i,j,k
    real(kind=dp_t) :: x,y,z,speed

    if (icomp .ne. 1 .and. icomp .ne. 2 .and. icomp .ne. 3) then
       call bl_error('physbc_macvel_3d requires icomp = 1, 2, or 3')
    end if

    if (ng_s .lt. 1) then
       ! no need to do anything
       return
    end if

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. REFLECT_ODD) then
       if (icomp .eq. 2) then
          ! staggered y-velocity
          do i=1,ng_s
             s(lo(1)-i,lo(2)-ng_s:hi(2)+1+ng_s,lo(3)-ng_s:hi(3)+ng_s) = 2.d0*wallspeed_lo(1,1) &
                  - s(lo(1)+i-1,lo(2)-ng_s:hi(2)+1+ng_s,lo(3)-ng_s:hi(3)+ng_s)
          end do
       else if (icomp .eq. 3) then
          ! staggered z-velocity
          do i=1,ng_s
             s(lo(1)-i,lo(2)-ng_s:hi(2)+ng_s,lo(3)-ng_s:hi(3)+1+ng_s) = 2.d0*wallspeed_lo(2,1) &
                  - s(lo(1)+i-1,lo(2)-ng_s:hi(2)+ng_s,lo(3)-ng_s:hi(3)+1+ng_s)
          end do
       else
          print *,'physbc_macvel_3d: bc(1,1) = REFLECT_ODD for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(1,1) .eq. REFLECT_EVEN) then
       if (icomp .eq. 2) then
          ! staggered y-velocity
          do i=1,ng_s
             s(lo(1)-i,lo(2)-ng_s:hi(2)+1+ng_s,lo(3)-ng_s:hi(3)+ng_s) = &
                  s(lo(1)+i-1,lo(2)-ng_s:hi(2)+1+ng_s,lo(3)-ng_s:hi(3)+ng_s)
          end do
       else if (icomp .eq. 3) then
          ! staggered z-velocity
          do i=1,ng_s
             s(lo(1)-i,lo(2)-ng_s:hi(2)+ng_s,lo(3)-ng_s:hi(3)+1+ng_s) = &
                  s(lo(1)+i-1,lo(2)-ng_s:hi(2)+ng_s,lo(3)-ng_s:hi(3)+1+ng_s)
          end do
       else
          print *,'physbc_macvel_3d: bc(1,1) = REFLECT_EVEN for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(1,1) .eq. EXT_DIR) then
       if (icomp .eq. 1) then
          ! staggered x-velocity
          do i=1,ng_s
             s(lo(1)-i,lo(2)-ng_s:hi(2)+ng_s,lo(3)-ng_s:hi(3)+ng_s) = &
                  s(lo(1),lo(2)-ng_s:hi(2)+ng_s,lo(3)-ng_s:hi(3)+ng_s)
          end do
       else
          print *,'physbc_macvel_3d: bc(1,1) = EXT_DIR for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if

    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else 
       print *,'physbc_macvel_3d: bc(1,1) =',bc(1,1),' for icomp =',icomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. REFLECT_ODD) then
       if (icomp .eq. 2) then
          ! staggered y-velocity
          do i=1,ng_s
             s(hi(1)+i,lo(2)-ng_s:hi(2)+1+ng_s,lo(3)-ng_s:hi(3)+ng_s) = 2.d0*wallspeed_hi(1,1) &
                  - s(hi(1)+1-i,lo(2)-ng_s:hi(2)+1+ng_s,lo(3)-ng_s:hi(3)+ng_s)
          end do
       else if (icomp .eq. 3) then
          ! staggered z-velocity
          do i=1,ng_s
             s(hi(1)+i,lo(2)-ng_s:hi(2)+ng_s,lo(3)-ng_s:hi(3)+1+ng_s) = 2.d0*wallspeed_hi(2,1) &
                  - s(hi(1)+1-i,lo(2)-ng_s:hi(2)+ng_s,lo(3)-ng_s:hi(3)+1+ng_s)
          end do
       else
          print *,'physbc_macvel_3d: bc(1,2) = REFLECT_ODD for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(1,2) .eq. REFLECT_EVEN) then
       if (icomp .eq. 2) then
          ! staggered y-velocity
          do i=1,ng_s
             s(hi(1)+i,lo(2)-ng_s:hi(2)+1+ng_s,lo(3)-ng_s:hi(3)+ng_s) = &
                  s(hi(1)+1-i,lo(2)-ng_s:hi(2)+1+ng_s,lo(3)-ng_s:hi(3)+ng_s)
          end do
       else if (icomp .eq. 3) then
          ! staggered z-velocity
          do i=1,ng_s
             s(hi(1)+i,lo(2)-ng_s:hi(2)+ng_s,lo(3)-ng_s:hi(3)+1+ng_s) = &
                  s(hi(1)+1-i,lo(2)-ng_s:hi(2)+ng_s,lo(3)-ng_s:hi(3)+1+ng_s)
          end do
       else
          print *,'physbc_macvel_3d: bc(1,2) = REFLECT_EVEN for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(1,2) .eq. EXT_DIR) then
       if (icomp .eq. 1) then
          ! staggered x-velocity
          do i=1,ng_s
             s(hi(1)+1+i,lo(2)-ng_s:hi(2)+ng_s,lo(3)-ng_s:hi(3)+ng_s) = &
                  s(hi(1)+1,lo(2)-ng_s:hi(2)+ng_s,lo(3)-ng_s:hi(3)+ng_s)  
          end do
       else
          print *,'physbc_macvel_3d: bc(1,2) = EXT_DIR for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else 
       print *,'physbc_macvel_3d: bc(1,2) =',bc(1,2),' for icomp =',icomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. REFLECT_ODD) then
       if (icomp .eq. 1) then
          ! staggered x-velocity
          do j=1,ng_s
             s(lo(1)-ng_s:hi(1)+1+ng_s,lo(2)-j,lo(3)-ng_s:hi(3)+ng_s) = 2.d0*wallspeed_lo(1,2) &
                  - s(lo(1)-ng_s:hi(1)+1+ng_s,lo(2)+j-1,lo(3)-ng_s:hi(3)+ng_s)
          end do
       else if (icomp .eq. 3) then
          ! staggered z-velocity
          do j=1,ng_s
             s(lo(1)-ng_s:hi(1)+ng_s,lo(2)-j,lo(3)-ng_s:hi(3)+1+ng_s) = 2.d0*wallspeed_lo(2,2) &
                  - s(lo(1)-ng_s:hi(1)+ng_s,lo(2)+j-1,lo(3)-ng_s:hi(3)+1+ng_s)
          end do
       else
          print *,'physbc_macvel_3d: bc(2,1) = REFLECT_ODD for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(2,1) .eq. REFLECT_EVEN) then
       if (icomp .eq. 1) then
          ! staggered x-velocity
          do j=1,ng_s
             s(lo(1)-ng_s:hi(1)+1+ng_s,lo(2)-j,lo(3)-ng_s:hi(3)+ng_s) = &
                  s(lo(1)-ng_s:hi(1)+1+ng_s,lo(2)+j-1,lo(3)-ng_s:hi(3)+ng_s)
          end do
       else if (icomp .eq. 3) then
          ! staggered z-velocity
          do j=1,ng_s
             s(lo(1)-ng_s:hi(1)+ng_s,lo(2)-j,lo(3)-ng_s:hi(3)+1+ng_s) = &
                  s(lo(1)-ng_s:hi(1)+ng_s,lo(2)+j-1,lo(3)-ng_s:hi(3)+1+ng_s)
          end do
       else
          print *,'physbc_macvel_3d: bc(2,1) = REFLECT_EVEN for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(2,1) .eq. EXT_DIR) then
       if (icomp .eq. 2) then
          ! staggered y-velocity
          do j=1,ng_s
             s(lo(1)-ng_s:hi(1)+ng_s,lo(2)-j,lo(3)-ng_s:hi(3)+ng_s) = &
                  s(lo(1)-ng_s:hi(1)+ng_s,lo(2),lo(3)-ng_s:hi(3)+ng_s)  
          end do
       else
          print *,'physbc_macvel_3d: bc(2,1) = EXT_DIR for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else 
       print *,'physbc_macvel_3d: bc(2,1) =',bc(2,1),' for icomp =',icomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. REFLECT_ODD) then
       if (icomp .eq. 1) then
          ! staggered x-velocity
          ! only apply boyce_bc if x and/or z is non-periodic
          if (boyce_bc .and. ( (.not. pmask(1)) .or. (.not. pmask(3)) ) ) then
             if ( (.not. pmask(1)) .and. (.not. pmask(3)) ) then
                ! if we are non-periodic in both x and z, taper velociy near x and z walls
                do i=lo(1)-ng_s,hi(1)+1+ng_s
                   x = prob_lo(1) + dble(i)*dx(1)
                   do k=lo(3)-ng_s,hi(3)+ng_s
                      z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
                      speed = wallspeed_hi(1,2)* &
                           0.5d0*(1.d0 + sin(2.d0*M_PI*(x-prob_lo(1))/prob_hi(1) - 0.5d0*M_PI))* &
                           0.5d0*(1.d0 + sin(2.d0*M_PI*(z-prob_lo(3))/prob_hi(3) - 0.5d0*M_PI))
                      ! original stencil
                      do j=1,ng_s
                         s(i,hi(2)+j,k) = 2.d0*speed - s(i,hi(2)+1-j,k)
                      end do
                      ! higher-order stencil
                      ! s(i,hi(2)+1,k) = (8.d0/3.d0)*speed - 2.d0*s(i,hi(2),k) + (1.d0/3.d0)*s(i,hi(2)-1,k)
                   end do
                end do
             else if (.not. pmask(1)) then
                ! if we are non-periodic in only x, taper velocities near x-walls
                do i=lo(1)-ng_s,hi(1)+1+ng_s
                   x = prob_lo(1) + dble(i)*dx(1)
                   speed = wallspeed_hi(1,2)* &
                        0.5d0*(1.d0 + sin(2.d0*M_PI*(x-prob_lo(1))/prob_hi(1) - 0.5d0*M_PI))
                   do k=lo(3)-ng_s,hi(3)+ng_s
                      do j=1,ng_s
                         s(i,hi(2)+j,k) = 2.d0*speed - s(i,hi(2)+1-j,k)
                      end do
                   end do
                end do
             else
                ! if we are non-periodic in only z, taper velocities near z-walls
                do k=lo(3)-ng_s,hi(3)+ng_s
                   z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
                   do i=lo(1)-ng_s,hi(1)+1+ng_s
                      speed = wallspeed_hi(1,2)* &
                           0.5d0*(1.d0 + sin(2.d0*M_PI*(z-prob_lo(3))/prob_hi(3) - 0.5d0*M_PI))
                      do j=1,ng_s
                         s(i,hi(2)+j,k) = 2.d0*speed - s(i,hi(2)+1-j,k)
                      end do
                   end do
                end do
             end if
          else
             do j=1,ng_s
                s(lo(1)-ng_s:hi(1)+1+ng_s,hi(2)+j,lo(3)-ng_s:hi(3)+ng_s) = 2.d0*wallspeed_hi(1,2) &
                     - s(lo(1)-ng_s:hi(1)+1+ng_s,hi(2)+1-j,lo(3)-ng_s:hi(3)+ng_s)
             end do
          end if
       else if (icomp .eq. 3) then
          ! staggered z-velocity
          ! only apply boyce_bc if x and/or z is non-periodic
          if (boyce_bc .and. ( (.not. pmask(1)) .or. (.not. pmask(3)) ) ) then
             if ( (.not. pmask(1)) .and. (.not. pmask(3)) ) then
                ! if we are non-periodic in both x and z, taper velociy near x and z walls
                do k=lo(3)-ng_s,hi(3)+1+ng_s
                   z = prob_lo(3) + dble(k)*dx(3)
                   do i=lo(1)-ng_s,hi(1)+ng_s
                      x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                      speed = wallspeed_hi(2,2)* &
                           0.5d0*(1.d0 + sin(2.d0*M_PI*(z-prob_lo(3))/prob_hi(3) - 0.5d0*M_PI))* &
                           0.5d0*(1.d0 + sin(2.d0*M_PI*(x-prob_lo(1))/prob_hi(1) - 0.5d0*M_PI))
                      ! original stencil
                      do j=1,ng_s
                         s(i,hi(2)+j,k) = 2.d0*speed - s(i,hi(2)+1-j,k)
                      end do
                      ! higher-order stencil
                      ! s(i,hi(2)+1,k) = (8.d0/3.d0)*speed - 2.d0*s(i,hi(2),k) + (1.d0/3.d0)*s(i,hi(2)-1,k)
                   end do
                end do
             else if (.not. pmask(1)) then
                ! if we are non-periodic in only x, taper velocities near x-walls
                do i=lo(1)-ng_s,hi(1)+ng_s
                   x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                   do k=lo(3)-ng_s,hi(3)+1+ng_s
                      speed = wallspeed_hi(2,2)* &
                           0.5d0*(1.d0 + sin(2.d0*M_PI*(x-prob_lo(1))/prob_hi(1) - 0.5d0*M_PI))
                      do j=1,ng_s
                         s(i,hi(2)+j,k) = 2.d0*speed - s(i,hi(2)+1-j,k)
                      end do
                   end do
                end do
             else
                ! if we are non-periodic in only z, taper velocities near z-walls
                do k=lo(3)-ng_s,hi(3)+1+ng_s
                   z = prob_lo(3) + dble(k)*dx(3)
                   speed = wallspeed_hi(2,2)* &
                        0.5d0*(1.d0 + sin(2.d0*M_PI*(z-prob_lo(3))/prob_hi(3) - 0.5d0*M_PI))
                   do i=lo(1)-ng_s,hi(1)+ng_s
                      do j=1,ng_s
                         s(i,hi(2)+j,k) = 2.d0*speed - s(i,hi(2)+1-j,k)
                      end do
                   end do
                end do
             end if
          else
             do j=1,ng_s
                s(lo(1)-ng_s:hi(1)+ng_s,hi(2)+j,lo(3)-ng_s:hi(3)+1+ng_s) = 2.d0*wallspeed_hi(2,2) &
                     - s(lo(1)-ng_s:hi(1)+ng_s,hi(2)+1-j,lo(3)-ng_s:hi(3)+1+ng_s)
             end do
          end if
       else
          print *,'physbc_macvel_3d: bc(2,2) = REFLECT_ODD for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(2,2) .eq. REFLECT_EVEN) then
       if (icomp .eq. 1) then
          ! staggered x-velocity
          do j=1,ng_s
             s(lo(1)-ng_s:hi(1)+1+ng_s,hi(2)+j,lo(3)-ng_s:hi(3)+ng_s) = &
                  s(lo(1)-ng_s:hi(1)+1+ng_s,hi(2)+1-j,lo(3)-ng_s:hi(3)+ng_s)
          end do
       else if (icomp .eq. 3) then
          ! staggered z-velocity
          do j=1,ng_s
             s(lo(1)-ng_s:hi(1)+ng_s,hi(2)+j,lo(3)-ng_s:hi(3)+1+ng_s) = &
                  s(lo(1)-ng_s:hi(1)+ng_s,hi(2)+1-j,lo(3)-ng_s:hi(3)+1+ng_s)
          end do
       else
          print *,'physbc_macvel_3d: bc(2,2) = REFLECT_EVEN for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(2,2) .eq. EXT_DIR) then
       if (icomp .eq. 2) then
          ! staggered y-velocity
          do j=1,ng_s
             s(lo(1)-ng_s:hi(1)+ng_s,hi(2)+1+j,lo(3)-ng_s:hi(3)+ng_s) = &
                  s(lo(1)-ng_s:hi(1)+ng_s,hi(2)+1,lo(3)-ng_s:hi(3)+ng_s)
          end do
       else
          print *,'physbc_macvel_3d: bc(2,2) = EXT_DIR for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else 
       print *,'physbc_macvel_3d: bc(2,2) =',bc(2,2),' for icomp =',icomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,1) .eq. REFLECT_ODD) then
       if (icomp .eq. 1) then
          ! staggered x-velocity
          do k=1,ng_s
             s(lo(1)-ng_s:hi(1)+1+ng_s,lo(2)-ng_s:hi(2)+ng_s,lo(3)-k) = 2.d0*wallspeed_lo(1,3) &
                  - s(lo(1)-ng_s:hi(1)+1+ng_s,lo(2)-ng_s:hi(2)+ng_s,lo(3)+k-1)
          end do
       else if (icomp .eq. 2) then
          ! staggered y-velocity
          do k=1,ng_s
             s(lo(1)-ng_s:hi(1)+ng_s,lo(2)-ng_s:hi(2)+1+ng_s,lo(3)-k) = 2.d0*wallspeed_lo(2,3) &
                  - s(lo(1)-ng_s:hi(1)+ng_s,lo(2)-ng_s:hi(2)+1+ng_s,lo(3)+k-1)
          end do
       else
          print *,'physbc_macvel_3d: bc(3,1) = REFLECT_ODD for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(3,1) .eq. REFLECT_EVEN) then
       if (icomp .eq. 1) then
          ! staggered x-velocity
          do k=1,ng_s
             s(lo(1)-ng_s:hi(1)+1+ng_s,lo(2)-ng_s:hi(2)+ng_s,lo(3)-k) = &
                  s(lo(1)-ng_s:hi(1)+1+ng_s,lo(2)-ng_s:hi(2)+ng_s,lo(3)+k-1)
          end do
       else if (icomp .eq. 2) then
          ! staggered y-velocity
          do k=1,ng_s
             s(lo(1)-ng_s:hi(1)+ng_s,lo(2)-ng_s:hi(2)+1+ng_s,lo(3)-k) = &
                  s(lo(1)-ng_s:hi(1)+ng_s,lo(2)-ng_s:hi(2)+1+ng_s,lo(3)+k-1)
          end do
       else
          print *,'physbc_macvel_3d: bc(3,1) = REFLECT_EVEN for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(3,1) .eq. EXT_DIR) then
       if (icomp .eq. 3) then
          ! staggered z-velocity
          do k=1,ng_s
             s(lo(1)-ng_s:hi(1)+ng_s,lo(2)-ng_s:hi(2)+ng_s,lo(3)-k) = &
                  s(lo(1)-ng_s:hi(1)+ng_s,lo(2)-ng_s:hi(2)+ng_s,lo(3)) 
          end do
       else
          print *,'physbc_macvel_3d: bc(3,1) = EXT_DIR for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(3,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else 
       print *,'physbc_macvel_3d: bc(3,1) =',bc(3,1),' for icomp =',icomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,2) .eq. REFLECT_ODD) then
       if (icomp .eq. 1) then
          ! staggered x-velocity
          ! only apply boyce_bc if x and/or y is non-periodic
          if (boyce_bc .and. ( (.not. pmask(1)) .or. (.not. pmask(2)) ) ) then
             if ( (.not. pmask(1)) .and. (.not. pmask(2)) ) then
                ! if we are non-periodic in both x and y, taper velociy near x and y walls
                do i=lo(1)-ng_s,hi(1)+1+ng_s
                   x = prob_lo(1) + dble(i)*dx(1)
                   do j=lo(2)-ng_s,hi(2)+ng_s
                      y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
                      speed = wallspeed_hi(1,3)* &
                           0.5d0*(1.d0 + sin(2.d0*M_PI*(x-prob_lo(1))/prob_hi(1) - 0.5d0*M_PI))* &
                           0.5d0*(1.d0 + sin(2.d0*M_PI*(y-prob_lo(2))/prob_hi(2) - 0.5d0*M_PI))
                      do k=1,ng_s
                         s(i,j,hi(3)+k) = 2.d0*speed - s(i,j,hi(3)+1-k)
                      end do
                   end do
                end do
             else if (.not. pmask(1)) then
                ! if we are non-periodic in only x, taper velocities near x-walls
                do i=lo(1)-ng_s,hi(1)+1+ng_s
                   x = prob_lo(1) + dble(i)*dx(1)
                   speed = wallspeed_hi(1,3)* &
                        0.5d0*(1.d0 + sin(2.d0*M_PI*(x-prob_lo(1))/prob_hi(1) - 0.5d0*M_PI))
                   do j=lo(2)-ng_s,hi(2)+ng_s
                      do k=1,ng_s
                         s(i,j,hi(3)+k) = 2.d0*speed - s(i,j,hi(3)+1-k)
                      end do
                   end do
                end do
             else
                ! if we are non-periodic in only y, taper velocities near y-walls
                do j=lo(2)-ng_s,hi(2)+ng_s
                   y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
                   do i=lo(1)-ng_s,hi(1)+1+ng_s
                      speed = wallspeed_hi(1,3)* &
                           0.5d0*(1.d0 + sin(2.d0*M_PI*(y-prob_lo(2))/prob_hi(2) - 0.5d0*M_PI))
                      do k=1,ng_s
                         s(i,j,hi(3)+k) = 2.d0*speed - s(i,j,hi(3)+1-k)
                      end do
                   end do
                end do
             end if
          else
             do k=1,ng_s
                s(lo(1)-ng_s:hi(1)+1+ng_s,lo(2)-ng_s:hi(2)+ng_s,hi(3)+k) = 2.d0*wallspeed_hi(1,3) &
                     - s(lo(1)-ng_s:hi(1)+1+ng_s,lo(2)-ng_s:hi(2)+ng_s,hi(3)+1-k)
             end do
          end if
       else if (icomp .eq. 2) then
          ! staggered y-velocity
          ! only apply boyce_bc if x and/or y is non-periodic
          if (boyce_bc .and. ( (.not. pmask(1)) .or. (.not. pmask(2)) ) ) then
             if ( (.not. pmask(1)) .and. (.not. pmask(2)) ) then
                ! if we are non-periodic in both x and y, taper velociy near x and y walls
                do j=lo(2)-ng_s,hi(2)+1+ng_s
                   y = prob_lo(2) + dble(j)*dx(2)
                   do i=lo(1)-ng_s,hi(1)+ng_s
                      x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                      speed = wallspeed_hi(2,3)* &
                           0.5d0*(1.d0 + sin(2.d0*M_PI*(y-prob_lo(2))/prob_hi(2) - 0.5d0*M_PI))* &
                           0.5d0*(1.d0 + sin(2.d0*M_PI*(x-prob_lo(1))/prob_hi(1) - 0.5d0*M_PI))
                      do k=1,ng_s
                         s(i,j,hi(3)+k) = 2.d0*speed - s(i,j,hi(3)+1-k)
                      end do
                   end do
                end do
             else if (.not. pmask(1)) then
                ! if we are non-periodic in only x, taper velocities near x-walls
                do i=lo(1)-ng_s,hi(1)+ng_s
                   x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                   do j=lo(2)-ng_s,hi(2)+1+ng_s
                      speed = wallspeed_hi(2,3)* &
                           0.5d0*(1.d0 + sin(2.d0*M_PI*(x-prob_lo(1))/prob_hi(1) - 0.5d0*M_PI))
                      do k=1,ng_s
                         s(i,j,hi(3)+k) = 2.d0*speed - s(i,j,hi(3)+1-k)
                      end do
                   end do
                end do
             else
                ! if we are non-periodic in only y, taper velocities near y-walls
                do j=lo(2)-ng_s,hi(2)+1+ng_s
                   y = prob_lo(2) + dble(j)*dx(2)
                   speed = wallspeed_hi(2,3)* &
                        0.5d0*(1.d0 + sin(2.d0*M_PI*(y-prob_lo(2))/prob_hi(2) - 0.5d0*M_PI))
                   do i=lo(1)-ng_s,hi(1)+ng_s
                      do k=1,ng_s
                         s(i,j,hi(3)+k) = 2.d0*speed - s(i,j,hi(3)+1-k)
                      end do
                   end do
                end do
             end if
          else
             do k=1,ng_s
                s(lo(1)-ng_s:hi(1)+ng_s,lo(2)-ng_s:hi(2)+1+ng_s,hi(3)+k) = 2.d0*wallspeed_hi(2,3) &
                     - s(lo(1)-ng_s:hi(1)+ng_s,lo(2)-ng_s:hi(2)+1+ng_s,hi(3)+1-k)
             end do
          end if
       else
          print *,'physbc_macvel_3d: bc(3,2) = REFLECT_ODD for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(3,2) .eq. REFLECT_EVEN) then
       if (icomp .eq. 1) then
          ! staggered x-velocity
          do k=1,ng_s
             s(lo(1)-ng_s:hi(1)+1+ng_s,lo(2)-ng_s:hi(2)+ng_s,hi(3)+k) = &
                  s(lo(1)-ng_s:hi(1)+1+ng_s,lo(2)-ng_s:hi(2)+ng_s,hi(3)+1-k)
          end do
       else if (icomp .eq. 2) then
          ! staggered y-velocity
          do k=1,ng_s
             s(lo(1)-ng_s:hi(1)+ng_s,lo(2)-ng_s:hi(2)+1+ng_s,hi(3)+k) = &
                  s(lo(1)-ng_s:hi(1)+ng_s,lo(2)-ng_s:hi(2)+1+ng_s,hi(3)+1-k)
          end do
       else
          print *,'physbc_macvel_3d: bc(3,2) = REFLECT_EVEN for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(3,2) .eq. EXT_DIR) then
       if (icomp .eq. 3) then
          ! staggered z-velocity
          do k=1,ng_s
             s(lo(1)-ng_s:hi(1)+ng_s,lo(2)-ng_s:hi(2)+ng_s,hi(3)+1+k) = &
                  s(lo(1)-ng_s:hi(1)+ng_s,lo(2)-ng_s:hi(2)+ng_s,hi(3)+1) 
          end do
       else
          print *,'physbc_macvel_3d: bc(3,2) = EXT_DIR for icomp =',icomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(3,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else 
       print *,'physbc_macvel_3d: bc(3,2) =',bc(3,2),' for icomp =',icomp
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine physbc_macvel_3d

  subroutine multifab_physbc_domainvel(s,start_scomp,start_bccomp,num_comp,the_bc_level)

    type(multifab) , intent(inout) :: s
    integer        , intent(in   ) :: start_scomp,start_bccomp,num_comp
    type(bc_level) , intent(in   ) :: the_bc_level

    ! Local
    integer                  :: lo(get_dim(s)),hi(get_dim(s))
    integer                  :: i,ng_s,dm,scomp,bccomp
    real(kind=dp_t), pointer :: sp(:,:,:,:)

    ng_s = nghost(s)

    dm = get_dim(s)
    
    do i=1,nfabs(s)
       sp => dataptr(s,i)
       lo = lwb(get_box(s,i))
       hi = upb(get_box(s,i))
       select case (dm)
       case (2)
          do scomp = start_scomp,start_scomp+num_comp-1
             bccomp = start_bccomp + scomp - start_scomp
             call physbc_domainvel_2d(sp(:,:,1,scomp), ng_s, lo, hi, &
                                      the_bc_level%adv_bc_level_array(i,:,:,bccomp),bccomp)
          end do
       case (3)
          do scomp = start_scomp,start_scomp+num_comp-1
             bccomp = start_bccomp + scomp - start_scomp
             call physbc_domainvel_3d(sp(:,:,:,scomp), ng_s, lo, hi, &
                                      the_bc_level%adv_bc_level_array(i,:,:,bccomp),bccomp)
          end do
       end select
    end do
 
  end subroutine multifab_physbc_domainvel

  subroutine physbc_domainvel_2d(s,ng_s,lo,hi,bc,icomp)

    use bl_constants_module
    use bc_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng_s:,lo(2)-ng_s:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: icomp

    if (icomp .ne. 1 .and. icomp .ne. 2) then
       call bl_error('physbc_domainvel_2d requires icomp = 1 or 2')
    end if

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. EXT_DIR) then
       if (icomp .eq. 1) then
          ! staggered x-velocity
          s(lo(1),lo(2)-ng_s:hi(2)+ng_s) = 0.d0
       end if
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. EXT_DIR) then
       if (icomp .eq. 1) then
          ! staggered x-velocity
          s(hi(1)+1,lo(2)-ng_s:hi(2)+ng_s) = 0.d0
       end if
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. EXT_DIR) then
       if (icomp .eq. 2) then
          ! staggered y-velocity
          s(lo(1)-ng_s:hi(1)+ng_s,lo(2)) = 0.d0
       end if
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. EXT_DIR) then
       if (icomp .eq. 2) then
          ! staggered y-velocity
          s(lo(1)-ng_s:hi(1)+ng_s,hi(2)+1) = 0.d0
       end if
    end if

  end subroutine physbc_domainvel_2d

  subroutine physbc_domainvel_3d(s,ng_s,lo,hi,bc,icomp)

    use bl_constants_module
    use bc_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: icomp

    if (icomp .ne. 1 .and. icomp .ne. 2 .and. icomp .ne. 3) then
       call bl_error('physbc_domainvel_3d requires icomp = 1, 2, or 3')
    end if

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. EXT_DIR) then
       if (icomp .eq. 1) then
          ! staggered x-velocity
          s(lo(1),lo(2)-ng_s:hi(2)+ng_s,lo(3)-ng_s:hi(3)+ng_s) = 0.d0
       end if
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. EXT_DIR) then
       if (icomp .eq. 1) then
          ! staggered x-velocity
          s(hi(1)+1,lo(2)-ng_s:hi(2)+ng_s,lo(3)-ng_s:hi(3)+ng_s) = 0.d0
       end if
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. EXT_DIR) then
       if (icomp .eq. 2) then
          ! staggered y-velocity
          s(lo(1)-ng_s:hi(1)+ng_s,lo(2),lo(3)-ng_s:hi(3)+ng_s) = 0.d0
       end if
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. EXT_DIR) then
       if (icomp .eq. 2) then
          ! staggered y-velocity
          s(lo(1)-ng_s:hi(1)+ng_s,hi(2)+1,lo(3)-ng_s:hi(3)+ng_s) = 0.d0
       end if
    end if

!!!!!!!!!!!!!!!!!!
! lo-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,1) .eq. EXT_DIR) then
       if (icomp .eq. 3) then
          ! staggered z-velocity
         s(lo(1)-ng_s:hi(1)+ng_s,lo(2)-ng_s:hi(2)+ng_s,lo(3)) = 0.d0
      end if
   end if

!!!!!!!!!!!!!!!!!!
! hi-z boundary
!!!!!!!!!!!!!!!!!!

   if (bc(3,2) .eq. EXT_DIR) then
       if (icomp .eq. 3) then
          ! staggered z-velocity
         s(lo(1)-ng_s:hi(1)+ng_s,lo(2)-ng_s:hi(2)+ng_s,hi(3)+1) = 0.d0
      end if
   end if

 end subroutine physbc_domainvel_3d

end module multifab_physbc_module
