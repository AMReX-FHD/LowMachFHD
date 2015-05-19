module init_energy_module

  use ml_layout_module
  use probin_energy_module
  use probin_common_module, only: prob_type
  use probin_multispecies_module, only: nspecies, c_init, T_init
 
  implicit none

  private

  public :: init_energy


contains

  subroutine init_energy(mla,umac,rhotot,rho,rhoh,Temp,p0)

    ! initialize umac, rho, rho_i, rhoh, Temp, and p0
    ! what you supply is problem dependent, except that concentrations are required
    ! so we can enforce that sum(c_i)=1 by overwriting the final concentration,

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(inout) :: rhotot(:)
    type(multifab) , intent(inout) :: rho(:)
    type(multifab) , intent(inout) :: rhoh(:)
    type(multifab) , intent(inout) :: Temp(:)
    real(kind=dp_t), intent(inout) :: p0
 
    ! local variables
    integer                        :: lo(mla%dim), hi(mla%dim)
    integer                        :: n,nlevs,i,dm,ng_u,ng_1,ng_2,ng_3,ng_4
    real(kind=dp_t), pointer       ::  up(:,:,:,:)
    real(kind=dp_t), pointer       ::  vp(:,:,:,:)
    real(kind=dp_t), pointer       ::  wp(:,:,:,:)
    real(kind=dp_t), pointer       :: dp1(:,:,:,:)
    real(kind=dp_t), pointer       :: dp2(:,:,:,:)
    real(kind=dp_t), pointer       :: dp3(:,:,:,:)
    real(kind=dp_t), pointer       :: dp4(:,:,:,:)

    dm = mla%dim
    nlevs = mla%nlevel

    ng_u = umac(1,1)%ng
    ng_1 = rhotot(1)%ng
    ng_2 = rho(1)%ng
    ng_3 = rhoh(1)%ng
    ng_4 = Temp(1)%ng

    ! looping over boxes 
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          up  => dataptr(umac(n,1),i)
          vp  => dataptr(umac(n,2),i)
          dp1 => dataptr(rhotot(n),i)
          dp2 => dataptr(   rho(n),i)
          dp3 => dataptr(  rhoh(n),i)
          dp4 => dataptr(  Temp(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          select case (dm)
          case (2)
             call init_energy_2d(up(:,:,1,1),vp(:,:,1,1),ng_u, &
                                 dp1(:,:,1,1),ng_1,dp2(:,:,1,:),ng_2, &
                                 dp3(:,:,1,1),ng_3,dp4(:,:,1,1),ng_4, &
                                 p0,lo,hi)
          case (3)
             wp => dataptr(umac(n,3),i)
             call init_energy_3d(up(:,:,:,1),vp(:,:,:,1),wp(:,:,:,1),ng_u, &
                                 dp1(:,:,:,1),ng_1,dp2(:,:,:,:),ng_2, &
                                 dp3(:,:,:,1),ng_3,dp4(:,:,:,1),ng_4, &
                                 p0,lo,hi)
          end select
       end do
    end do

  end subroutine init_energy

  subroutine init_energy_2d(umac,vmac,ng_u,rhotot,ng_1,rho,ng_2,rhoh,ng_3,Temp,ng_4,p0,lo,hi)

    integer         :: lo(2), hi(2), ng_u, ng_1, ng_2, ng_3, ng_4
    real(kind=dp_t) ::   umac(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t) ::   vmac(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t) :: rhotot(lo(1)-ng_1:,lo(2)-ng_1:)
    real(kind=dp_t) ::    rho(lo(1)-ng_2:,lo(2)-ng_2:,:)
    real(kind=dp_t) ::   rhoh(lo(1)-ng_3:,lo(2)-ng_3:)
    real(kind=dp_t) ::   Temp(lo(1)-ng_4:,lo(2)-ng_4:)
    real(kind=dp_t) :: p0

 
    ! local varables
    integer :: i,j,n
    real(kind=dp_t) :: sum,hk(nspecies)

    integer :: iwrk
    real(kind=dp_t) :: rwrk

    select case (abs(prob_type))
    
    case (1)

       !=============================================================
       ! Provide (p0,T,w)
       ! w and T are constant and supplied by c_init(1,:) and T_init(1)
       ! zero velocity
       !=============================================================
 
       p0 = p0_in

       umac = 0.d0
       vmac = 0.d0

       ! specify constant concentration and temperature
       ! store concentrations in rho for now
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rho(i,j,1:nspecies) = c_init(1,1:nspecies)
             Temp(i,j) = T_init(1)
          end do
       end do

       ! overwrite final concentration so they sum to 1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             sum = 0
             do n=1,nspecies-1
                sum = sum + rho(i,j,n)
             end do
             rho(i,j,nspecies) = 1.d0 - sum
          end do
       end do

       ! compute density from the EOS
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             call CKRHOY(p0,Temp(i,j),rho(i,j,:),iwrk,rwrk,rhotot(i,j))
          end do
       end do

       ! compute enthalpy from the EOS and compute rhoh
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             call CKHMS(Temp(i,j),iwrk,rwrk,hk)
             rhoh(i,j) = 0.d0
             do n=1,nspecies
                ! note: rho still holds concentrations
                rhoh(i,j) = rhoh(i,j) + rho(i,j,n)*hk(n)
             end do
             rhoh(i,j) = rhoh(i,j)*rhotot(i,j)
          end do
       end do

       ! convert concentrations to rho_i
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             do n=1,nspecies
                rho(i,j,n) = rho(i,j,n)*rhotot(i,j)
             end do
          end do
       end do       

    case default

       call bl_error("Desired prob_type not supported in 2D")

    end select

  end subroutine init_energy_2d

  subroutine init_energy_3d(umac,vmac,wmac,ng_u,rhotot,ng_1,rho,ng_2,rhoh,ng_3,Temp,ng_4, &
                            p0,lo,hi)

    integer         :: lo(3), hi(3), ng_u, ng_1, ng_2, ng_3, ng_4
    real(kind=dp_t) ::   umac(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t) ::   vmac(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t) ::   wmac(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t) :: rhotot(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
    real(kind=dp_t) ::    rho(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_2:,:)
    real(kind=dp_t) ::   rhoh(lo(1)-ng_3:,lo(2)-ng_3:,lo(3)-ng_3:)
    real(kind=dp_t) ::   Temp(lo(1)-ng_4:,lo(2)-ng_4:,lo(3)-ng_4:)
    real(kind=dp_t) :: p0

 
    ! local varables
    integer :: i,j,k,n
    real(kind=dp_t) :: sum,hk(nspecies)

    integer :: iwrk
    real(kind=dp_t) :: rwrk

    select case (abs(prob_type))
    
    case (1)

       !=============================================================
       ! Provide (p0,T,w)
       ! w and T are constant and supplied by c_init(1,:) and T_init(1)
       ! zero velocity
       !=============================================================
 
       p0 = p0_in
       
       umac = 0.d0
       vmac = 0.d0
       wmac = 0.d0

       ! specify constant concentration and temperature
       ! store concentrations in rho for now
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rho(i,j,k,1:nspecies) = c_init(1,1:nspecies)
                Temp(i,j,k) = T_init(1)
             end do
          end do
       end do

       ! overwrite final concentration so they sum to 1
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                sum = 0
                do n=1,nspecies-1
                   sum = sum + rho(i,j,k,n)
                end do
                rho(i,j,k,nspecies) = 1.d0 - sum
             end do
          end do
       end do

       ! compute density from the EOS
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                call CKRHOY(p0,Temp(i,j,k),rho(i,j,k,:),iwrk,rwrk,rhotot(i,j,k))
             end do
          end do
       end do

       ! compute enthalpy from the EOS and compute rhoh
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                call CKHMS(Temp(i,j,k),iwrk,rwrk,hk)
                rhoh(i,j,k) = 0.d0
                do n=1,nspecies
                   ! note: rho still holds concentrations
                   rhoh(i,j,k) = rhoh(i,j,k) + rho(i,j,k,n)*hk(n)
                end do
                rhoh(i,j,k) = rhoh(i,j,k)*rhotot(i,j,k)
             end do
          end do
       end do

       ! convert concentrations to rho_i
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                do n=1,nspecies
                   rho(i,j,k,n) = rho(i,j,k,n)*rhotot(i,j,k)
                end do
             end do
          end do
       end do

    case default

       call bl_error("Desired prob_type not supported in 2D")

    end select

  end subroutine init_energy_3d

end module init_energy_module
