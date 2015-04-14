module initial_projection_module

  use ml_layout_module
  use bndry_reg_module
  use ml_solve_module
  use multifab_physbc_module
  use bc_module
  use define_bc_module
  use energy_eos_module
  use energy_eos_wrapper_module
  use convert_rhoc_to_c_module
  use convert_rhoh_to_h_module
  use mass_fluxdiv_energy_module
  use rhoh_fluxdiv_energy_module
  use macproject_module
  use convert_stag_module
  use div_and_grad_module
  use mk_advective_s_fluxdiv_module
  use mass_flux_utilities_module
  use probin_multispecies_module, only: nspecies
  use probin_common_module, only: n_cells

  use fabio_module

  implicit none

  private

  public :: initial_projection

contains

  subroutine initial_projection(mla,umac,rho,rhotot,rhoh,p0, &
                                gradp_baro,Temp,dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(inout) :: rho(:)
    type(multifab) , intent(inout) :: rhotot(:)
    type(multifab) , intent(inout) :: rhoh(:)
    real(kind=dp_t), intent(inout) :: p0
    type(multifab) , intent(inout) :: gradp_baro(:,:)
    type(multifab) , intent(inout) :: Temp(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables






  end subroutine initial_projection

end module initial_projection_module
