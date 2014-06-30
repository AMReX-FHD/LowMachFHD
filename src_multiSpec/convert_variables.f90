module convert_variables_module

  use multifab_module
  use ml_layout_module
  use probin_common_module, only: rhobar
  use probin_multispecies_module, only: nspecies

  implicit none

  private

  public :: convert_rho_to_c
  
contains

  subroutine convert_rho_to_c(mla,rho,c,rho_to_c)
    
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rho(:)
    type(multifab) , intent(inout) ::   c(:)
    logical        , intent(in   ) :: rho_to_c

    ! local
    integer :: nlevs,n,i,dm,ng_r,ng_c
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: rp(:,:,:,:)
    real(kind=dp_t), pointer :: cp(:,:,:,:)

    nlevs = mla%nlevel
    dm = mla%dim

    ng_r = rho(1)%ng
    ng_c =   c(1)%ng

    if (.not.(rho_to_c)) then
       if (ng_c .lt. ng_c) then
          call bl_error('convert_variables, not enough ghost cells in prim')
       end if
    end if
    
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          rp => dataptr(rho(n), i)
          cp => dataptr(  c(n), i)
          lo =  lwb(get_box(rho(n), i))
          hi =  upb(get_box(rho(n), i))
          select case (dm)
          case (2)
             if (rho_to_c) then
                ! rho to c - NO GHOST CELLS
                call rho_to_c_2d(rp(:,:,1,:),ng_r,cp(:,:,1,:),ng_c,lo,hi)
             else
                ! c to rho - INCLUDING GHOST CELLS
                call c_to_rho_2d(rp(:,:,1,:),ng_r,cp(:,:,1,:),ng_c,lo,hi)
             end if
          case (3)
             if (rho_to_c) then
                ! rho to c - NO GHOST CELLS
                call rho_to_c_3d(rp(:,:,:,:),ng_r,cp(:,:,:,:),ng_c,lo,hi)
             else
                ! c to rho - INCLUDING GHOST CELLS
                call c_to_rho_3d(rp(:,:,:,:),ng_r,cp(:,:,:,:),ng_c,lo,hi)
             end if
          end select
       end do
    end do
    
  end subroutine convert_rho_to_c

  subroutine rho_to_c_2d(rho,ng_r,c,ng_c,lo,hi)

    ! rho to c - NO GHOST CELLS

    integer        , intent(in   ) :: lo(:), hi(:), ng_r, ng_c
    real(kind=dp_t), intent(in   ) :: rho(lo(1)-ng_r:,lo(2)-ng_r:,:)
    real(kind=dp_t), intent(inout) ::   c(lo(1)-ng_c:,lo(2)-ng_c:,:)

    ! local
    integer :: i,j,n

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          do n=1,nspecies
             c(i,j,n) = rho(i,j,n)/rhobar(n)
          end do
       end do
    end do

  end subroutine rho_to_c_2d

  subroutine rho_to_c_3d(rho,ng_r,c,ng_c,lo,hi)

    ! rho to c - NO GHOST CELLS

    integer        , intent(in   ) :: lo(:), hi(:), ng_r, ng_c
    real(kind=dp_t), intent(in   ) :: rho(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:,:)
    real(kind=dp_t), intent(inout) ::   c(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:,:)

    ! local
    integer :: i,j,k,n

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             do n=1,nspecies
                c(i,j,k,n) = rho(i,j,k,n)/rhobar(n)
             end do
          end do
       end do
    end do

  end subroutine rho_to_c_3d

  subroutine c_to_rho_2d(rho,ng_r,c,ng_c,lo,hi)

    ! c to rho - INCLUDING GHOST CELLS
    
    integer        , intent(in   ) :: lo(:), hi(:), ng_r, ng_c
    real(kind=dp_t), intent(inout) :: rho(lo(1)-ng_r:,lo(2)-ng_r:,:)
    real(kind=dp_t), intent(in   ) ::   c(lo(1)-ng_c:,lo(2)-ng_c:,:)

    ! local
    integer :: i,j,n

    do j=lo(2)-ng_r,hi(2)+ng_r
       do i=lo(1)-ng_r,hi(1)+ng_r
          do n=1,nspecies
             rho(i,j,n) = c(i,j,n)*rhobar(n)
          end do
       end do
    end do

  end subroutine c_to_rho_2d

  subroutine c_to_rho_3d(rho,ng_r,c,ng_c,lo,hi)

    ! c to rho - INCLUDING GHOST CELLS

    integer        , intent(in   ) :: lo(:), hi(:), ng_r, ng_c
    real(kind=dp_t), intent(inout) :: rho(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:,:)
    real(kind=dp_t), intent(in   ) ::   c(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:,:)

    ! local
    integer :: i,j,k,n

    do k=lo(3)-ng_r,hi(3)+ng_r
       do j=lo(2)-ng_r,hi(2)+ng_r
          do i=lo(1)-ng_r,hi(1)+ng_r
             do n=1,nspecies
                rho(i,j,k,n) = c(i,j,k,n)*rhobar(n)
             end do
          end do
       end do
    end do

  end subroutine c_to_rho_3d

end module convert_variables_module
