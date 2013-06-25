module diffusive_flux_module

  use multifab_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use div_and_grad_module

  implicit none

  private

  public :: diffusive_flux

contains

  ! Donev: This needs to be changed substantially to include multiple components in rho and flux
  ! We can keep this routine single-level and call it level-by-level in advance.
  
  subroutine diffusive_flux(rho,flux,dx,the_bc_tower)

    type(multifab) , intent(in   ) :: rho
    type(multifab) , intent(inout) :: flux(:)
    real(kind=dp_t), intent(in   ) :: dx
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: lo(rho%dim), hi(rho%dim)
    integer :: dm, ng_p, ng_f, i

    real(kind=dp_t), pointer ::  pp(:,:,:,:)
    real(kind=dp_t), pointer :: fxp(:,:,:,:)
    real(kind=dp_t), pointer :: fyp(:,:,:,:)
    real(kind=dp_t), pointer :: fzp(:,:,:,:)

    dm   = rho%dim     ! dimensionality
    ng_p = rho%ng      ! number of ghost cells for rho
    ng_f = flux(1)%ng  ! number of ghost cells for flux 

    ! Donev:
    ! You should use compute_grad from src_common/div_and_grad.f90
    ! to compute the gradient of rho (later, of mole fractions x)
    ! and then from that compute the fluxes

    do i=1,nfabs(rho)  ! loop over owned boxes 
       pp  => dataptr(rho,i)
       fxp => dataptr(flux(1),i)
       fyp => dataptr(flux(2),i)
       lo = lwb(get_box(rho,i))
       hi = upb(get_box(rho,i))
       select case(dm)
       case (2)
          ! Donev: Here you need something like pp(:,:,1,1:nspecies)
          ! and similarly for fluxes
          ! You should pass the gradients to the flux routine
          ! You will also need to pass diffusion coefficients eventually
          
          call compute_flux_2d(pp(:,:,1,1), ng_p, &
                               fxp(:,:,1,1),  fyp(:,:,1,1), ng_f, &
                               lo, hi, dx, &
                               the_bc_tower%bc_tower_array(1)%adv_bc_level_array(i,:,:,1))
       case (3)
          fzp => dataptr(flux(3),i)
          call compute_flux_3d(pp(:,:,:,1), ng_p, &
                               fxp(:,:,:,1),  fyp(:,:,:,1), fzp(:,:,:,1), ng_f, &
                               lo, hi, dx, &
                               the_bc_tower%bc_tower_array(1)%adv_bc_level_array(i,:,:,1))
       end select
    end do

  end subroutine compute_flux

  subroutine compute_flux_2d(rho, ng_p, fluxx, fluxy, ng_f, lo, hi, dx, adv_bc)

    integer          :: lo(2), hi(2), ng_p, ng_f
    ! Donev: These arrays need to be made one rank higher for multispecies
    double precision ::   rho(lo(1)-ng_p:,lo(2)-ng_p:)
    double precision :: fluxx(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision :: fluxy(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision :: dx
    integer          :: adv_bc(:,:)

    ! local variables
    integer i,j

    
  end subroutine compute_flux_2d

  subroutine compute_flux_3d(rho, ng_p, fluxx, fluxy, fluxz, ng_f, &
                             lo, hi, dx, adv_bc)

    integer          :: lo(3), hi(3), ng_p, ng_f
    double precision ::   rho(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    double precision :: fluxx(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: fluxy(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: fluxz(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: dx
    integer          :: adv_bc(:,:)

    ! local variables
    integer i,j,k

  end subroutine compute_flux_3d

end module diffusive_flux_module

