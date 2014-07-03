module fill_rho_ghost_cells_module

  use multifab_module
  use define_bc_module
  use bc_module
  use probin_binarylm_module, only: c_bc
  use probin_common_module, only: rhobar

  implicit none

  private

  public :: fill_rho_ghost_cells

contains

  subroutine fill_rho_ghost_cells(prim,the_bc_level)

    type(multifab) , intent(inout) :: prim
    type(bc_level) , intent(in   ) :: the_bc_level

    ! local
    integer :: i,dm,ng
    integer :: lo(get_dim(prim)),hi(get_dim(prim))
    real(kind=dp_t), pointer :: pp(:,:,:,:)

    dm = get_dim(prim)
    ng = prim%ng

    do i=1,nfabs(prim)
       pp => dataptr(prim,i)
       lo = lwb(get_box(prim,i))
       hi = upb(get_box(prim,i))
       select case (dm)
       case (2)
          call fill_rho_ghost_cells_2d(pp(:,:,1,:),ng,lo,hi, &
                                       the_bc_level%adv_bc_level_array(i,:,:,scal_bc_comp))
       case (3)
          call fill_rho_ghost_cells_3d(pp(:,:,:,:),ng,lo,hi, &
                                       the_bc_level%adv_bc_level_array(i,:,:,scal_bc_comp))
       end select
    end do

  end subroutine fill_rho_ghost_cells

  subroutine fill_rho_ghost_cells_2d(prim,ng,lo,hi,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) :: prim(lo(1)-ng:,lo(2)-ng:,:)
    integer        , intent(in   ) :: bc(:,:)

    ! local
    integer :: i,j
    real(kind=dp_t) :: c_bc

    if (bc(1,1) .eq. FOEXTRAP .or. bc(1,1) .eq. HOEXTRAP .or. bc(1,1) .eq. EXT_DIR) then
       do j=lo(2)-ng,hi(2)+ng
          c_bc = prim(lo(1)-1,j,2)
          prim(lo(1)-ng:lo(1)-1,j,1) = 1.d0/(c_bc/rhobar(1) + (1.d0-c_bc)/rhobar(2))
       end do
    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'fill_rho_ghost_cells_2d: bc(1,1) =',bc(1,1)
       call bl_error('NOT SUPPORTED')
    end if

    if (bc(1,2) .eq. FOEXTRAP .or. bc(1,2) .eq. HOEXTRAP .or. bc(1,2) .eq. EXT_DIR) then
       do j=lo(2)-ng,hi(2)+ng
          c_bc = prim(hi(1)+1,j,2)
          prim(hi(1)+1:hi(1)+ng,j,1) = 1.d0/(c_bc/rhobar(1) + (1.d0-c_bc)/rhobar(2))
       end do
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'fill_rho_ghost_cells_2d: bc(1,1) =',bc(1,1)
       call bl_error('NOT SUPPORTED')
    end if

    if (bc(2,1) .eq. FOEXTRAP .or. bc(2,1) .eq. HOEXTRAP .or. bc(2,1) .eq. EXT_DIR) then
       do i=lo(1)-ng,hi(1)+ng
          c_bc = prim(i,lo(2)-1,2)
          prim(i,lo(2)-ng:lo(2)-1,1) = 1.d0/(c_bc/rhobar(1) + (1.d0-c_bc)/rhobar(2))
       end do
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'fill_rho_ghost_cells_2d: bc(2,1) =',bc(1,1)
       call bl_error('NOT SUPPORTED')
    end if

    if (bc(2,2) .eq. FOEXTRAP .or. bc(2,2) .eq. HOEXTRAP .or. bc(2,2) .eq. EXT_DIR) then
       do i=lo(1)-ng,hi(1)+ng
          c_bc = prim(i,hi(2)+1,2)
          prim(i,hi(2)+1:hi(2)+ng,1) = 1.d0/(c_bc/rhobar(1) + (1.d0-c_bc)/rhobar(2))
       end do
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'fill_rho_ghost_cells_2d: bc(2,2) =',bc(1,1)
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine fill_rho_ghost_cells_2d

  subroutine fill_rho_ghost_cells_3d(prim,ng,lo,hi,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) :: prim(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    integer        , intent(in   ) :: bc(:,:)

    ! local
    integer :: i,j,k
    real(kind=dp_t) :: c_bc

    if (bc(1,1) .eq. FOEXTRAP .or. bc(1,1) .eq. HOEXTRAP .or. bc(1,1) .eq. EXT_DIR) then
       do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          c_bc = prim(lo(1)-1,j,k,2)
          prim(lo(1)-ng:lo(1)-1,j,k,1) = 1.d0/(c_bc/rhobar(1) + (1.d0-c_bc)/rhobar(2))
       end do
       end do
    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'fill_rho_ghost_cells_3d: bc(1,1) =',bc(1,1)
       call bl_error('NOT SUPPORTED')
    end if

    if (bc(1,2) .eq. FOEXTRAP .or. bc(1,2) .eq. HOEXTRAP .or. bc(1,2) .eq. EXT_DIR) then
       do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          c_bc = prim(hi(1)+1,j,k,2)
          prim(hi(1)+1:hi(1)+ng,j,k,1) = 1.d0/(c_bc/rhobar(1) + (1.d0-c_bc)/rhobar(2))
       end do
       end do
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'fill_rho_ghost_cells_3d: bc(1,1) =',bc(1,1)
       call bl_error('NOT SUPPORTED')
    end if

    if (bc(2,1) .eq. FOEXTRAP .or. bc(2,1) .eq. HOEXTRAP .or. bc(2,1) .eq. EXT_DIR) then
       do k=lo(3)-ng,hi(3)+ng
       do i=lo(1)-ng,hi(1)+ng
          c_bc = prim(i,lo(2)-1,k,2)
          prim(i,lo(2)-ng:lo(2)-1,k,1) = 1.d0/(c_bc/rhobar(1) + (1.d0-c_bc)/rhobar(2))
       end do
       end do
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'fill_rho_ghost_cells_3d: bc(2,1) =',bc(1,1)
       call bl_error('NOT SUPPORTED')
    end if

    if (bc(2,2) .eq. FOEXTRAP .or. bc(2,2) .eq. HOEXTRAP .or. bc(2,2) .eq. EXT_DIR) then
       do k=lo(3)-ng,hi(3)+ng
       do i=lo(1)-ng,hi(1)+ng
          c_bc = prim(i,hi(2)+1,k,2)
          prim(i,hi(2)+1:hi(2)+ng,k,1) = 1.d0/(c_bc/rhobar(1) + (1.d0-c_bc)/rhobar(2))
       end do
       end do
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'fill_rho_ghost_cells_3d: bc(2,2) =',bc(1,1)
       call bl_error('NOT SUPPORTED')
    end if

    if (bc(3,1) .eq. FOEXTRAP .or. bc(3,1) .eq. HOEXTRAP .or. bc(3,1) .eq. EXT_DIR) then
       do j=lo(2)-ng,hi(2)+ng
       do i=lo(1)-ng,hi(1)+ng
          c_bc = prim(i,j,lo(3)-1,2)
          prim(i,j,lo(3)-ng:lo(3)-1,1) = 1.d0/(c_bc/rhobar(1) + (1.d0-c_bc)/rhobar(2))
       end do
       end do
    else if (bc(3,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'fill_rho_ghost_cells_3d: bc(3,1) =',bc(1,1)
       call bl_error('NOT SUPPORTED')
    end if

    if (bc(3,2) .eq. FOEXTRAP .or. bc(3,2) .eq. HOEXTRAP .or. bc(3,2) .eq. EXT_DIR) then
       do j=lo(2)-ng,hi(2)+ng
       do i=lo(1)-ng,hi(1)+ng
          c_bc = prim(i,j,hi(3)+1,2)
          prim(i,j,hi(3)+1:hi(3)+ng,1) = 1.d0/(c_bc/rhobar(1) + (1.d0-c_bc)/rhobar(2))
       end do
       end do
    else if (bc(3,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'fill_rho_ghost_cells_3d: bc(3,2) =',bc(1,1)
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine fill_rho_ghost_cells_3d

end module fill_rho_ghost_cells_module
