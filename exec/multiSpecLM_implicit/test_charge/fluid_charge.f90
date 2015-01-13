module fluid_charge_module
 
  use ml_layout_module
  use convert_stag_module
  use define_bc_module
  use probin_multispecies_module, only: charge_per_mass, nspecies, charge_per_mass, rho_part_bc_comp
  use probin_common_module, only: molmass, k_B

  implicit none

  private

  public :: compute_total_charge, compute_charge_coef, momentum_charge_force
  
contains

  ! compute total charge = rho y^T dot z = rho_i dot z
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
      
      integer          :: lo(:),hi(:),ng_1,ng_2
      real(kind=dp_t)  ::    rho(lo(1)-ng_1:,lo(2)-ng_1:,:)
      real(kind=dp_t)  :: charge(lo(1)-ng_2:,lo(2)-ng_2:)

      ! local variables
      integer :: i,j,comp

      do j=lo(2)-ng_2,hi(2)+ng_2
         do i=lo(1)-ng_2,hi(1)+ng_2

            charge(i,j) = 0.d0
            do comp=1,nspecies
               charge(i,j) = charge(i,j) + rho(i,j,comp)*charge_per_mass(comp)
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

      do k=lo(3)-ng_2,hi(3)+ng_2
         do j=lo(2)-ng_2,hi(2)+ng_2
            do i=lo(1)-ng_2,hi(1)+ng_2

               charge(i,j,k) = 0.d0
               do n=1,nspecies
                  charge(i,j,k) = charge(i,j,k) + rho(i,j,k,n)*charge_per_mass(n)
               end do

            end do
         end do
      end do

    end subroutine compute_total_charge_3d

  end subroutine compute_total_charge

  ! compute cell-centered mass diffusion coefficients due to charge fluid
  ! charge_coef = (rho/(n k_B T)) (z - charge*vector_of_ones)
  subroutine compute_charge_coef(mla,rho,rhotot,Temp,charge,charge_coef)

    type(ml_layout), intent(in   ) :: mla
    type(multifab ), intent(in   ) :: rho(:)
    type(multifab ), intent(in   ) :: rhotot(:)
    type(multifab ), intent(in   ) :: Temp(:)
    type(multifab ), intent(in   ) :: charge(:)
    type(multifab ), intent(in   ) :: charge_coef(:)

    ! local variables
    integer :: i,n,dm,nlevs
    integer :: ng_1,ng_2,ng_3,ng_4,ng_5
    integer :: lo(mla%dim),hi(mla%dim)

    ! pointers into multifabs
    real(kind=dp_t), pointer :: dp1(:,:,:,:)
    real(kind=dp_t), pointer :: dp2(:,:,:,:)
    real(kind=dp_t), pointer :: dp3(:,:,:,:)
    real(kind=dp_t), pointer :: dp4(:,:,:,:)
    real(kind=dp_t), pointer :: dp5(:,:,:,:)

    dm = mla%dim
    nlevs = mla%nlevel

    ng_1 = rho(1)%ng
    ng_2 = rhotot(1)%ng
    ng_3 = Temp(1)%ng
    ng_4 = charge(1)%ng
    ng_5 = charge_coef(1)%ng

    do n=1,nlevs
       do i=1,nfabs(charge_coef(n))
          dp1 => dataptr(rho(n),i)
          dp2 => dataptr(rhotot(n),i)
          dp3 => dataptr(Temp(n),i)
          dp4 => dataptr(charge(n),i)
          dp5 => dataptr(charge_coef(n),i)
          lo = lwb(get_box(charge_coef(n),i))
          hi = upb(get_box(charge_coef(n),i))
          select case (dm)
          case (2)
             call compute_charge_coef_2d(dp1(:,:,1,:),ng_1, &
                                         dp2(:,:,1,1),ng_2, &
                                         dp3(:,:,1,1),ng_3, &
                                         dp4(:,:,1,1),ng_4, &
                                         dp5(:,:,1,:),ng_5, lo,hi)
          case (3)
             call compute_charge_coef_3d(dp1(:,:,:,:),ng_1, &
                                         dp2(:,:,:,1),ng_2, &
                                         dp3(:,:,:,1),ng_3, &
                                         dp4(:,:,:,1),ng_4, &
                                         dp5(:,:,:,:),ng_5, lo,hi)
          end select
       end do
    end do

  contains

    subroutine compute_charge_coef_2d(rho,ng_1,rhotot,ng_2,Temp,ng_3, &
                                      charge,ng_4,charge_coef,ng_5,lo,hi)
      
      integer         :: lo(:),hi(:),ng_1,ng_2,ng_3,ng_4,ng_5
      real(kind=dp_t) ::         rho(lo(1)-ng_1:,lo(2)-ng_1:,:)
      real(kind=dp_t) ::      rhotot(lo(1)-ng_2:,lo(2)-ng_2:)
      real(kind=dp_t) ::        Temp(lo(1)-ng_3:,lo(2)-ng_3:)
      real(kind=dp_t) ::      charge(lo(1)-ng_4:,lo(2)-ng_4:)
      real(kind=dp_t) :: charge_coef(lo(1)-ng_5:,lo(2)-ng_5:,:)

      ! local variables
      integer :: i,j,comp
      real(kind=dp_t) :: n

      do j=lo(2)-1,hi(2)+1
      do i=lo(1)-1,hi(1)+1

         n = 0.d0
         do comp=1,nspecies
            n = n + rho(i,j,comp)/molmass(comp)
         end do
            
         do comp=1,nspecies
            charge_coef(i,j,comp) = (rhotot(i,j)/(n*k_B*Temp(i,j))) &
                 * (charge_per_mass(comp) - charge(i,j))
         end do

      end do
      end do

    end subroutine compute_charge_coef_2d

    subroutine compute_charge_coef_3d(rho,ng_1,rhotot,ng_2,Temp,ng_3, &
                                      charge,ng_4,charge_coef,ng_5,lo,hi)
      
      integer         :: lo(:),hi(:),ng_1,ng_2,ng_3,ng_4,ng_5
      real(kind=dp_t) ::         rho(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:,:)
      real(kind=dp_t) ::      rhotot(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_2:)
      real(kind=dp_t) ::        Temp(lo(1)-ng_3:,lo(2)-ng_3:,lo(3)-ng_3:)
      real(kind=dp_t) ::      charge(lo(1)-ng_4:,lo(2)-ng_4:,lo(3)-ng_4:)
      real(kind=dp_t) :: charge_coef(lo(1)-ng_5:,lo(2)-ng_5:,lo(3)-ng_5:,:)

      ! local variables
      integer :: i,j,k,comp
      real(kind=dp_t) :: n

      do k=lo(3)-1,hi(3)+1
      do j=lo(2)-1,hi(2)+1
      do i=lo(1)-1,hi(1)+1

         n = 0.d0
         do comp=1,nspecies
            n = n + rho(i,j,k,comp)/molmass(comp)
         end do
            
         do comp=1,nspecies
            charge_coef(i,j,k,comp) = (rhotot(i,j,k)/(n*k_B*Temp(i,j,k))) &
                 * (charge_per_mass(comp) - charge(i,j,k))
         end do

      end do
      end do
      end do

    end subroutine compute_charge_coef_3d

  end subroutine compute_charge_coef

  ! increment the momentum charge force by -(1/2)(charge * grad_Epot)^old - (1/2)(charge * grad_Epot)^new
  subroutine momentum_charge_force(mla,mom_charge_force,charge_old,charge_new,grad_Epot_old,grad_Epot_new, &
                                   the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab ), intent(inout) :: mom_charge_force(:,:)
    type(multifab ), intent(in   ) :: charge_old(:)
    type(multifab ), intent(in   ) :: charge_new(:)
    type(multifab ), intent(in   ) :: grad_Epot_old(:,:)
    type(multifab ), intent(in   ) :: grad_Epot_new(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: i,n,dm,nlevs

    type(multifab) :: charge_edge(mla%nlevel,mla%dim)

    nlevs = mla%nlevel
    dm = mla%dim

    do n=1,nlevs
       do i=1,dm
          call multifab_build_edge(charge_edge(n,i),mla%la(n),1,0,i)
       end do
    end do

    ! put charge^old on faces
    call average_cc_to_face(nlevs,charge_old,charge_edge,1,rho_part_bc_comp,1, &
                            the_bc_tower%bc_tower_array)

    ! increment mom_charge_force by -(charge * grad Epot)^old
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_c(charge_edge(n,i),1,grad_Epot_old(n,i),1,0)
          call multifab_mult_mult_s_c(charge_edge(n,i),1,0.5d0,1,0)
          call multifab_sub_sub_c(mom_charge_force(n,i),1,charge_edge(n,i),1,1,0)
       end do
    end do

    ! put charge^new on faces
    call average_cc_to_face(nlevs,charge_new,charge_edge,1,rho_part_bc_comp,1, &
                            the_bc_tower%bc_tower_array)

    ! increment mom_charge_force by -(charge * grad Epot)^new
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_c(charge_edge(n,i),1,grad_Epot_new(n,i),1,0)
          call multifab_mult_mult_s_c(charge_edge(n,i),1,0.5d0,1,0)
          call multifab_sub_sub_c(mom_charge_force(n,i),1,charge_edge(n,i),1,1,0)
       end do
    end do
    
    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(charge_edge(n,i))
       end do
    end do

  end subroutine momentum_charge_force


end module fluid_charge_module
