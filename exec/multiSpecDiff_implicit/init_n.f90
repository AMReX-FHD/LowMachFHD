module init_n_module

  use bl_types
  use ml_layout_module
  use probin_common_module, only: prob_lo, prob_hi, prob_type
  use probin_multispecies_module, only: nspecies
 
  implicit none

  private

  public :: init_n


contains

  subroutine init_n(mla,n_init,dx)

    ! initialize rho_i and umac in the valid region
    ! we first initialize c_i in the valid region
    ! then enforce that sum(c_i)=1 by overwriting the final concentration,
    ! and then use the EOS to compute rho_i

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: n_init(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
 
    ! local variables
    integer                        :: lo(mla%dim), hi(mla%dim)
    integer                        :: i, dm, n, nlevs, ng_n
    real(kind=dp_t), pointer       :: np(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "init_n")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_n = n_init(1)%ng

    ! looping over boxes 
    do n=1,nlevs
       do i=1,nfabs(n_init(n))
          np => dataptr(n_init(n),i)
          lo = lwb(get_box(n_init(n),i))
          hi = upb(get_box(n_init(n),i))
          select case (dm)
          case (2)
             call init_n_2d(np(:,:,1,:),ng_n,lo,hi,dx(n,:))
          case (3)
             call init_n_3d(np(:,:,:,:),ng_n,lo,hi,dx(n,:))
          end select
       end do
    end do

    call destroy(bpt)

  end subroutine init_n

  subroutine init_n_2d(n_init,ng_n,lo,hi,dx)

    integer          :: lo(:), hi(:), ng_n
    real(kind=dp_t)  :: n_init(lo(1)-ng_n:,lo(2)-ng_n:,:)
    real(kind=dp_t)  :: dx(:)
 
    ! local varables
    integer         :: i,j,comp
    real(kind=dp_t) :: x,y,r,cen(2),sum

    cen(1:2) = 0.5d0*(prob_lo(1:2)+prob_hi(1:2))

    do j=lo(2),hi(2)
       y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
       do i=lo(1),hi(1)
          x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)

          r = sqrt((x-cen(1))**2 + (y-cen(2))**2)

          n_init(i,j,1) = 0.5d0*exp(-r**2)
          n_init(i,j,2) = 0.25d0

       end do
    end do

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          sum = 0.d0
          do comp=1,nspecies-1
             sum = sum + n_init(i,j,comp)
          end do
          n_init(i,j,nspecies) = 1.d0-sum
       end do
    end do

  end subroutine init_n_2d

  subroutine init_n_3d(n_init,ng_n,lo,hi,dx)

    integer          :: lo(:), hi(:), ng_n
    real(kind=dp_t)  :: n_init(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:,:)
    real(kind=dp_t)  :: dx(:)
 
    ! local varables
    integer         :: i,j,k,comp
    real(kind=dp_t) :: x,y,z,r,cen(3),sum

    cen(1:3) = 0.5d0*(prob_lo(1:3)+prob_hi(1:3))

    do k=lo(3),hi(3)
       z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)

             r = sqrt((x-cen(1))**2 + (y-cen(2))**2 + (z-cen(3))**2)
             
             n_init(i,j,k,1) = 0.5d0*exp(-r**2)
             n_init(i,j,k,2) = 0.25d0

          end do
       end do
    end do

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             sum = 0.d0
             do comp=1,nspecies-1
                sum = sum + n_init(i,j,k,comp)
             end do
             n_init(i,j,k,nspecies) = 1.d0-sum
          end do
       end do
    end do

  end subroutine init_n_3d

end module init_n_module
