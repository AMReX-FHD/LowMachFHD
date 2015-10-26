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

    type(multifab) :: ext_src_d(mla%nlevel)
    type(multifab) :: ext_src_r(mla%nlevel)

    type(bl_prof_timer),save :: bpt

    call build(bpt,"advance_timestep")

    nlevs = mla%nlevel

    ! external source term for diffusion/reaction solvers for inhomogeneous bc algorithm
    do n=1,nlevs
       call multifab_build(ext_src_d(n),mla%la(n),nspecies,0)
       call multifab_build(ext_src_r(n),mla%la(n),nspecies,0)
       call setval(ext_src_d(n),0.d0,all=.true.)
       call setval(ext_src_r(n),0.d0,all=.true.)
    end do

    if (splitting_type .eq. 0) then
       ! D + R

       call advance_diffusion(mla,n_old,n_new,ext_src_d,dx,dt,the_bc_tower)
       call advance_reaction (mla,n_new,n_old,ext_src_r,dx,dt,the_bc_tower)  ! swap n_new/n_old to avoid calling copy()
       do n=1,nlevs
          call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,n_new(n)%ng) ! make sure n_new contains the new state
       end do

    else if (splitting_type .eq. 1) then
       ! (1/2)R + D + (1/2)R

       call advance_reaction (mla,n_old,n_new,ext_src_r,dx,0.5d0*dt,the_bc_tower)
       call advance_diffusion(mla,n_new,n_old,ext_src_d,dx,dt      ,the_bc_tower) ! swap n_new/n_old to avoid calling copy()
       call advance_reaction (mla,n_old,n_new,ext_src_r,dx,0.5d0*dt,the_bc_tower) ! swap n_new/n_old to avoid calling copy()

    else if (splitting_type .eq. 2) then
       ! (1/2)D + R + (1/2)D

       call advance_diffusion(mla,n_old,n_new,ext_src_d,dx,0.5d0*dt,the_bc_tower)
       call advance_reaction (mla,n_new,n_old,ext_src_r,dx,dt      ,the_bc_tower) ! swap n_new/n_old to avoid calling copy()
       call advance_diffusion(mla,n_old,n_new,ext_src_d,dx,0.5d0*dt,the_bc_tower) ! swap n_new/n_old to avoid calling copy()

    else if (splitting_type .eq. 3) then
       ! (1/2)D + R + (1/2)D with inhomogeneous boundary conditions
       ! under development and can eventually be merged in with splitting_type=2

       ! (1/2)D + R + (1/2)D

       call advance_diffusion(mla,n_old,n_new,ext_src_d,dx,0.5d0*dt,the_bc_tower)
       call advance_reaction (mla,n_new,n_old,ext_src_r,dx,dt      ,the_bc_tower) ! swap n_new/n_old to avoid calling copy()
       call advance_diffusion(mla,n_old,n_new,ext_src_d,dx,0.5d0*dt,the_bc_tower) ! swap n_new/n_old to avoid calling copy()

    else
       call bl_error("advance_timestep: invalid splitting_type")
    end if

    do n=1,nlevs
       call multifab_destroy(ext_src_d(n))
       call multifab_destroy(ext_src_r(n))
    end do


    call destroy(bpt)

  end subroutine advance_timestep

end module advance_timestep_module
