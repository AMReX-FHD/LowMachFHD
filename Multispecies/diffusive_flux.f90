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
 
  subroutine diffusive_flux(mla,rho,diff_coeffs,flux,dx,the_bc_level) 

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: rho(:)
    type(multifab) , intent(in   ) :: diff_coeffs(:)
    type(multifab) , intent(inout) :: flux(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local variables
    integer n,i,dm,nlevs

    dm    = mla%dim     ! dimensionality
    nlevs = mla%nlevel  ! number of levels 
 
    call compute_grad(mla,rho,flux,dx,1,scal_bc_comp,1,nspecies,the_bc_level)

    ! Multiply flux with diffusion constants 
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_c(flux(n,i),1,diff_coeffs(n),1,nspecies,0)
       end do
    end do
    
    ! Donev: If grad(temperature)
    !call compute_grad(mla,temperature,flux,dx,1,scal_bc_comp+n_species,n_species+1,1,the_bc_level)

  end subroutine diffusive_flux

end module diffusive_flux_module

