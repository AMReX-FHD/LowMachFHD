module convert_variables_module

  use multifab_module
  use ml_layout_module
  use probin_common_module, only: rhobar
  use probin_multispecies_module, only: nspecies

  implicit none

  private

  public :: convert_cons_to_prim
  
contains

  subroutine convert_cons_to_prim(mla,rho,c,cons_to_prim)
    
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rho(:)
    type(multifab) , intent(inout) ::   c(:)
    logical        , intent(in   ) :: cons_to_prim

    ! local
    integer :: nlevs,n,i,dm,ng_r,ng_c
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: rp(:,:,:,:)
    real(kind=dp_t), pointer :: cp(:,:,:,:)

    nlevs = mla%nlevel

    ng_r = rho(1)%ng
    ng_c =   c(1)%ng

    if (.not.(cons_to_prim)) then
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
             if (cons_to_prim) then
                ! cons to prim - NO GHOST CELLS

             else
                ! prim to cons - INCLUDING GHOST CELLS

             end if
          case (3)
             if (cons_to_prim) then
                ! cons to prim - NO GHOST CELLS

             else
                ! prim to cons - INCLUDING GHOST CELLS

             end if
          end select
       end do
    end do
    
  end subroutine convert_cons_to_prim

  subroutine cons_to_prim_2d(rho,ng_r,c,ng_c,lo,hi)

    ! cons to prim - NO GHOST CELLS

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

  end subroutine cons_to_prim_2d

  subroutine cons_to_prim_3d(rho,ng_r,c,ng_c,lo,hi)

    ! cons to prim - NO GHOST CELLS

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

  end subroutine cons_to_prim_3d

  subroutine prim_to_cons_2d(rho,ng_r,c,ng_c,lo,hi)

    ! prim to cons - INCLUDING GHOST CELLS
    
    integer        , intent(in   ) :: lo(:), hi(:), ng_r, ng_c
    real(kind=dp_t), intent(inout) :: rho(lo(1)-ng_r:,lo(2)-ng_r:,:)
    real(kind=dp_t), intent(in   ) ::   c(lo(1)-ng_c:,lo(2)-ng_c:,:)

    ! local
    integer :: i,j,n

    do j=lo(2)-ng_c,hi(2)+ng_c
       do i=lo(1)-ng_c,hi(1)+ng_c
          do n=1,nspecies
             rho(i,j,n) = c(i,j,n)*rhobar(n)
          end do
       end do
    end do

  end subroutine prim_to_cons_2d

  subroutine prim_to_cons_3d(rho,ng_r,c,ng_c,lo,hi)

    ! prim to cons - INCLUDING GHOST CELLS

    integer        , intent(in   ) :: lo(:), hi(:), ng_r, ng_c
    real(kind=dp_t), intent(inout) :: rho(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:,:)
    real(kind=dp_t), intent(in   ) ::   c(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:,:)

    ! local
    integer :: i,j,k,n

    do k=lo(3)-ng_c,hi(3)+ng_c
       do j=lo(2)-ng_c,hi(2)+ng_c
          do i=lo(1)-ng_c,hi(1)+ng_c
             do n=1,nspecies
                rho(i,j,k,n) = c(i,j,k,n)*rhobar(n)
             end do
          end do
       end do
    end do

  end subroutine prim_to_cons_3d

end module convert_variables_module
