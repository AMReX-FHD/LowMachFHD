module advance_timestep_module

  use ml_layout_module
  use define_bc_module
  use advance_diffusion_module

  implicit none

  private

  public :: advance_timestep

contains

  subroutine advance_timestep(mla,n_old,n_new,dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: n_old(:)
    type(multifab) , intent(inout) :: n_new(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    call advance_diffusion(mla,n_old,n_new,dx,dt,the_bc_tower)

  end subroutine advance_timestep

end module advance_timestep_module
