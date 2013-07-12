module diffusive_flux_module

  use multifab_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use div_and_grad_module
  use probin_multispecies_module
  use ml_layout_module

  implicit none

  private

  public :: diffusive_flux

contains
 
  subroutine diffusive_flux(mla,rho,flux,dx,the_bc_level) 

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: rho(:)
    type(multifab) , intent(inout) :: flux(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! Donev: or, do
    ! integer :: nspec
    ! nspec = rho%nc ! Number of components

    call compute_grad(mla,rho,flux,dx,1,scal_bc_comp,1,nspecies,the_bc_level)
    
    ! Donev: If there is a diff_coeff:
    !call multifab_mult_mult_s_c(flux,1,diff_coeff,nspecies,0)
    
    ! If there are nspecies diff_coeffs:
    ! do specie=1,nspecies
    !    call multifab_mult_mult_s_c(flux,specie,diff_coeffs(specie),1,0)
    ! end do
    
    ! Donev: If grad(temperature)
    !call compute_grad(mla,temperature,flux,dx,1,scal_bc_comp+n_species,n_species+1,1,the_bc_level)

  end subroutine diffusive_flux

end module diffusive_flux_module

