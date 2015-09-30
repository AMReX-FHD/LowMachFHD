module advance_timestep_module

  use ml_layout_module
  use define_bc_module
  use advance_diffusion_module
  use advance_reaction_module
  use probin_reactdiff_module, only: nspecies, splitting_type

  implicit none

  private

  public :: advance_timestep

contains

  subroutine advance_timestep(mla,n_old,n_new,dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: n_old(:)
    type(multifab) , intent(inout) :: n_new(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: n,nlevs

    type(bl_prof_timer),save :: bpt

    call build(bpt,"advance_timestep")

    nlevs = mla%nlevel

    if (splitting_type .eq. 0) then
       ! D + R

       call advance_diffusion(mla,n_old,n_new,dx,dt,the_bc_tower)
       call advance_reaction (mla,n_new,n_old,dx,dt,the_bc_tower)  ! swap n_new/n_old to avoid calling copy()
       do n=1,nlevs
          call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,n_new(n)%ng) ! make sure n_new contains the new state
       end do

    else if (splitting_type .eq. 1) then
       ! (1/2)R + D + (1/2)R

       call advance_reaction (mla,n_old,n_new,dx,0.5d0*dt,the_bc_tower)
       call advance_diffusion(mla,n_new,n_old,dx,      dt,the_bc_tower) ! swap n_new/n_old to avoid calling copy()
       call advance_reaction (mla,n_old,n_new,dx,0.5d0*dt,the_bc_tower) ! swap n_new/n_old to avoid calling copy()

    else if (splitting_type .eq. 2) then
       ! (1/2)D + R + (1/2)D

       call advance_diffusion(mla,n_old,n_new,dx,0.5d0*dt,the_bc_tower)
       call advance_reaction (mla,n_new,n_old,dx,      dt,the_bc_tower) ! swap n_new/n_old to avoid calling copy()
       call advance_diffusion(mla,n_old,n_new,dx,0.5d0*dt,the_bc_tower) ! swap n_new/n_old to avoid calling copy()

    else
       call bl_error("advance_timestep: invalid splitting_type")
    end if

    call destroy(bpt)

  end subroutine advance_timestep

end module advance_timestep_module
