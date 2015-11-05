module fluid_charge_module
 
  use ml_layout_module
  use convert_stag_module
  use define_bc_module
  use bc_module
  use probin_common_module, only: molmass, k_B, total_volume
  use probin_multispecies_module, only: nspecies
  use probin_charged_module, only: charge_per_mass

  implicit none

  private

  public :: compute_total_charge, compute_charge_coef, momentum_charge_force, &
            enforce_charge_neutrality
  
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
            charge_coef(i,j,comp) = (rho(i,j,comp)/(n*k_B*Temp(i,j))) &
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
            charge_coef(i,j,k,comp) = (rho(i,j,k,comp)/(n*k_B*Temp(i,j,k))) &
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
    call average_cc_to_face(nlevs,charge_old,charge_edge,1,c_bc_comp,1, &
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
    call average_cc_to_face(nlevs,charge_new,charge_edge,1,c_bc_comp,1, &
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

  subroutine enforce_charge_neutrality(mla,rho)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rho(:)

    ! local
    type(multifab) :: charge(mla%nlevel)

    integer :: i,dm,nlevs,n,ng_r
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t) :: charge_temp

    ! pointers into multifab
    real(kind=dp_t), pointer :: rp(:,:,:,:)

    nlevs = mla%nlevel
    dm = mla%dim

    ng_r = rho(1)%ng

    if (nlevs .ne. 1) then
       call bl_error("enforce_charge_neutrality only works with nlevs=1")
    end if

    do n=1,nlevs
       call multifab_build(charge(n),mla%la(n),1,0)
    end do

    ! compute total charge in each cell
    call compute_total_charge(mla,rho,charge)

    ! integrate charge over domain
    charge_temp = multifab_sum_c(charge(1),1,1)

    ! divide total charge by # of zones
!    charge_temp = charge_temp / total_volume

    print*,'charge before',charge_temp
    
    ! for positively charged zones, pick the positive species with the largest rho and
    ! subtract density.  Pick the negative species with the largest rho_i and add density
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          rp => dataptr(rho(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          select case (dm)
          case (2)
             call enforce_charge_neutrality_2d(rp(:,:,1,:),ng_r,lo,hi,charge_temp)
          case (3)
          end select
       end do
    end do

    ! compute total charge in each cell
    call compute_total_charge(mla,rho,charge)

    ! integrate charge over domain
    charge_temp = multifab_sum_c(charge(1),1,1)

    ! divide total charge by # of zones
!    charge_temp = charge_temp / total_volume

    print*,'charge after',charge_temp

    do n=1,nlevs
       call multifab_destroy(charge(n))
    end do

  contains

    subroutine enforce_charge_neutrality_2d(rho,ng_r,lo,hi,net_charge)
      
      integer         :: lo(:),hi(:),ng_r
      real(kind=dp_t) :: rho(lo(1)-ng_r:,lo(2)-ng_r:,:)
      real(kind=dp_t) :: net_charge

      ! local variables
      integer :: i,j,comp,negative_comp,positive_comp
      logical :: is_positive(nspecies),is_negative(nspecies)
      real(kind=dp_t) :: rho_temp

      is_positive = .false.
      is_negative = .false.
      do comp=1,nspecies
         if (charge_per_mass(comp) .gt. 0.d0) then
            is_positive(comp) = .true.
         end if
         if (charge_per_mass(comp) .lt. 0.d0) then
            is_negative(comp) = .true.
         end if
      end do

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)

            if (i .eq. 0 .and. j .eq. 0) then

            ! find the negatively charged species with largest rho_i
            rho_temp = -1.d0
            negative_comp = -1
            do comp=1,nspecies
               if (is_negative(comp)) then
                  if (rho(i,j,comp) .gt. rho_temp) then
                     rho_temp = rho(i,j,comp)
                     negative_comp = comp
                  end if
               end if
            end do

            ! find the positively charged species with largest rho_i
            rho_temp = -1.d0
            positive_comp = -1
            do comp=1,nspecies
               if (is_positive(comp)) then
                  if (rho(i,j,comp) .gt. rho_temp) then
                     rho_temp = rho(i,j,comp)
                     positive_comp = comp
                  end if
               end if
            end do

            if (net_charge .lt. 0.d0) then
               rho(i,j,negative_comp) = rho(i,j,negative_comp) &
                    - abs(0.5d0*net_charge/charge_per_mass(negative_comp))
               rho(i,j,positive_comp) = rho(i,j,positive_comp) &
                    + abs(0.5d0*net_charge/charge_per_mass(positive_comp))
            else
               rho(i,j,positive_comp) = rho(i,j,positive_comp) &
                    - abs(0.5d0*net_charge/charge_per_mass(positive_comp))
               rho(i,j,negative_comp) = rho(i,j,negative_comp) &
                    + abs(0.5d0*net_charge/charge_per_mass(negative_comp))
            end if

            end if

         end do
      end do

    end subroutine enforce_charge_neutrality_2d

  end subroutine enforce_charge_neutrality





end module fluid_charge_module
