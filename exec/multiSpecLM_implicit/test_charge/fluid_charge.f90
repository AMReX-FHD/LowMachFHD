module fluid_charge_module
 
  use ml_layout_module
  use probin_multispecies_module, only: charge_per_mass, nspecies, charge_per_mass

  implicit none

  private

  public :: compute_total_charge
  
contains

  subroutine compute_total_charge(mla,rho,charge)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: rho(:)
    type(multifab) , intent(inout) :: charge(:)

    ! local variables
    integer :: i,n,dm,nlevs
    integer :: ng_1,ng_2
    integer :: lo(mla%dim),hi(mla%dim)

    ! pointers into multifabs
    real(kind=dp_t), pointer :: dp1(:,:,:,:)
    real(kind=dp_t), pointer :: dp2(:,:,:,:)

    dm = mla%dim
    nlevs = mla%nlevel

    ng_1 = rho(1)%ng
    ng_2 = charge(1)%ng

    do n=1,nlevs
       do i=1,nfabs(charge(n))
          dp1 => dataptr(rho(n),i)
          dp2 => dataptr(charge(n),i)
          lo = lwb(get_box(charge(n),i))
          hi = upb(get_box(charge(n),i))
          select case (dm)
          case (2)
             call compute_total_charge_2d(dp1(:,:,1,:),ng_1,dp2(:,:,1,1),ng_2,lo,hi)
          case (3)
             call compute_total_charge_3d(dp1(:,:,:,:),ng_1,dp2(:,:,:,1),ng_2,lo,hi)
          end select
       end do
    end do

  contains

    subroutine compute_total_charge_2d(rho,ng_1,charge,ng_2,lo,hi)
      
      integer          :: lo(:),hi(:),ng_1,ng_2,ng_3
      real(kind=dp_t)  ::    rho(lo(1)-ng_1:,lo(2)-ng_1:,:)
      real(kind=dp_t)  :: charge(lo(1)-ng_2:,lo(2)-ng_2:)

      ! local variables
      integer :: i,j,n

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)

            charge(i,j) = 0.d0
            do n=1,nspecies
               charge(i,j) = charge(i,j) + rho(i,j,n)*charge_per_mass(n)
            end do

         end do
      end do

    end subroutine compute_total_charge_2d

    subroutine compute_total_charge_3d(rho,ng_1,charge,ng_2,lo,hi)
      
      integer          :: lo(:),hi(:),ng_1,ng_2
      real(kind=dp_t)  ::    rho(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:,:)
      real(kind=dp_t)  :: charge(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_2:)

      ! local variables
      integer :: i,j,k,n

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

               charge(i,j,k) = 0.d0
               do n=1,nspecies
                  charge(i,j,k) = charge(i,j,k) + rho(i,j,k,n)*charge_per_mass(n)
               end do

            end do
         end do
      end do

    end subroutine compute_total_charge_3d

  end subroutine compute_total_charge


  subroutine mass_charge_force(mla)

    type(ml_layout), intent(in   ) :: mla

  end subroutine mass_charge_force


  subroutine momentum_charge_force(mla)

    type(ml_layout), intent(in   ) :: mla

  end subroutine momentum_charge_force


end module fluid_charge_module
