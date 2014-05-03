module multifab_zero_edgeval_module

  use multifab_module
  use define_bc_module
  use bc_module
  use bl_error_module

  implicit none

  private

  public :: multifab_zero_edgeval

contains

  subroutine multifab_zero_edgeval(edge,start_comp,start_bccomp,num_comp,the_bc_level)

    ! vel_bc_n(nlevs,dm) are the normal velocities

    type(multifab) , intent(inout) :: edge(:)
    integer        , intent(in   ) :: start_comp,start_bccomp,num_comp
    type(bc_level) , intent(in   ) :: the_bc_level

    ! Local
    integer                  :: lo(get_dim(edge(1))),hi(get_dim(edge(1)))
    integer                  :: ng_e,i,dm,comp,bccomp
    real(kind=dp_t), pointer :: epx(:,:,:,:), epy(:,:,:,:), epz(:,:,:,:)

    ng_e = nghost(edge(1))
    dm = get_dim(edge(1))
    
    do i=1,nfabs(edge(1))
       epx => dataptr(edge(1),i)
       epy => dataptr(edge(2),i)
       lo = lwb(get_box(edge(1),i))
       hi = upb(get_box(edge(2),i))
       do comp=start_comp,start_comp+num_comp-1
          bccomp = start_bccomp + comp - start_comp
          select case (dm)
          case (2)
             call zero_edgeval_2d(epx(:,:,1,comp), epy(:,:,1,comp), ng_e, lo, hi, &
                                  the_bc_level%adv_bc_level_array(i,:,:,bccomp))
          case (3)
             epz => dataptr(edge(3),i)
             call zero_edgeval_3d(epx(:,:,:,comp), epy(:,:,:,comp), epz(:,:,:,comp), ng_e, lo, hi, &
                                  the_bc_level%adv_bc_level_array(i,:,:,bccomp))
          end select
       end do
    end do
 
  end subroutine multifab_zero_edgeval

  subroutine zero_edgeval_2d(edgex,edgey,ng_e,lo,hi,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng_e
    real(kind=dp_t), intent(inout) :: edgex(lo(1)-ng_e:,lo(2)-ng_e:)
    real(kind=dp_t), intent(inout) :: edgey(lo(1)-ng_e:,lo(2)-ng_e:)
    integer        , intent(in   ) :: bc(:,:)

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. FOEXTRAP .or. bc(1,1) .eq. HOEXTRAP .or. bc(1,1) .eq. EXT_DIR) then
       edgex(lo(1),lo(2):hi(2)) = 0.d0
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. FOEXTRAP .or. bc(1,2) .eq. HOEXTRAP .or. bc(1,2) .eq. EXT_DIR) then
       edgex(hi(1)+1,lo(2):hi(2)) = 0.d0
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. FOEXTRAP .or. bc(2,1) .eq. HOEXTRAP .or. bc(2,1) .eq. EXT_DIR) then
       edgey(lo(1):hi(1),lo(2)) = 0.d0
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. FOEXTRAP .or. bc(2,2) .eq. HOEXTRAP .or. bc(2,2) .eq. EXT_DIR) then
       edgey(lo(1):hi(1),hi(2)+1) = 0.d0
    end if

  end subroutine zero_edgeval_2d

  subroutine zero_edgeval_3d(edgex,edgey,edgez,ng_e,lo,hi,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng_e
    real(kind=dp_t), intent(inout) :: edgex(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(inout) :: edgey(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(inout) :: edgez(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    integer        , intent(in   ) :: bc(:,:)

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. FOEXTRAP .or. bc(1,1) .eq. HOEXTRAP .or. bc(1,1) .eq. EXT_DIR) then
       edgex(lo(1),lo(2):hi(2),lo(3):hi(3)) = 0.d0
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. FOEXTRAP .or. bc(1,2) .eq. HOEXTRAP .or. bc(1,2) .eq. EXT_DIR) then
       edgex(hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = 0.d0
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. FOEXTRAP .or. bc(2,1) .eq. HOEXTRAP .or. bc(2,1) .eq. EXT_DIR) then
       edgey(lo(1):hi(1),lo(2),lo(3):hi(3)) = 0.d0
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. FOEXTRAP .or. bc(2,2) .eq. HOEXTRAP .or. bc(2,2) .eq. EXT_DIR) then
       edgey(lo(1):hi(1),hi(2)+1,lo(3):hi(3)) = 0.d0
    end if

!!!!!!!!!!!!!!!!!!
! lo-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,1) .eq. FOEXTRAP .or. bc(3,1) .eq. HOEXTRAP .or. bc(3,1) .eq. EXT_DIR) then
       edgez(lo(1):hi(1),lo(2):hi(2),lo(3)) = 0.d0
    end if

!!!!!!!!!!!!!!!!!!
! hi-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,2) .eq. FOEXTRAP .or. bc(3,2) .eq. HOEXTRAP .or. bc(3,2) .eq. EXT_DIR) then
       edgez(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = 0.d0
    end if

  end subroutine zero_edgeval_3d

end module multifab_zero_edgeval_module
