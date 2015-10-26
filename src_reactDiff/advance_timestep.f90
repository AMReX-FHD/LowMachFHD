module advance_timestep_module

  use ml_layout_module
  use define_bc_module
  use advance_diffusion_module
  use advance_reaction_module
  use multifab_physbc_module
  use bc_module
  use probin_reactdiff_module, only: nspecies, splitting_type, n_bc, reaction_type, use_Poisson_rng

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
    integer :: n,nlevs,dm

    type(multifab) :: ext_src_d(mla%nlevel)
    type(multifab) :: ext_src_r(mla%nlevel)
    type(multifab) :: zerofab(mla%nlevel)
    type(multifab) :: z(mla%nlevel)

    ! n_i boundary conditions (dir,lohi,species)
    real(kind=dp_t) :: n_bc_temp(mla%dim,2,nspecies)

    integer :: reaction_type_temp, use_Poisson_rng_temp

    type(bl_prof_timer),save :: bpt

    call build(bpt,"advance_timestep")

    nlevs = mla%nlevel
    dm = mla%dim

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

       do n=1,nlevs
          call multifab_build(z(n),mla%la(n),nspecies,0)
          call multifab_build(zerofab(n),mla%la(n),nspecies,0)
          call setval(zerofab(n),0.d0,all=.true.)
       end do

       ! compute z with new-time boundary conditions
       !
       !

       do n=1,nlevs
          ! temporary hack to test n_bc=1 case
          ! in general we will solve div D_k grad z_k = 0
          call setval(z(n),1.d0,all=.true.)
       end do

       ! save boundary conditions
       n_bc_temp(1:dm,1:2,1:nspecies) = n_bc(1:dm,1:2,1:nspecies)

       ! make boundary conditions homogeneous
       n_bc(1:dm,1:2,1:nspecies) = 0.d0
       
       ! compute \tilde{n}
       do n=1,nlevs
          call multifab_sub_sub_c(n_old(n),1,z(n),1,nspecies,0)
          call multifab_fill_boundary(n_old(n))
          call multifab_physbc(n_old(n),1,scal_bc_comp,nspecies, &
                               the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
       end do

       ! setup reaction_type so the reactions are done deterministically
       reaction_type_temp = reaction_type
       use_Poisson_rng_temp = use_Poisson_rng

       reaction_type = 1
       use_Poisson_rng = -1

       ! store reactions rates for z in ext_src_d
       ! store (negative) reaction rates for z in ext_src_r
       call advance_reaction(mla,n_old,ext_src_d,zerofab,dx,dt,the_bc_tower,return_rates_in=.true.)
       do n=1,nlevs
          call multifab_copy_c(ext_src_r(n),1,ext_src_d(n),1,nspecies,0)
          call multifab_mult_mult_s_c(ext_src_r(n),1,-1.d0,nspecies,0)
       end do

       ! restore reaction_type
       reaction_type = reaction_type
       use_Poisson_rng = use_Poisson_rng

       ! advance diffusion
       call advance_diffusion(mla,n_old,n_new,ext_src_d,dx,0.5d0*dt,the_bc_tower)

       ! advance reaction
       call advance_reaction (mla,n_new,n_old,ext_src_r,dx,dt      ,the_bc_tower) ! swap n_new/n_old to avoid calling copy()

       ! advance diffusion
       call advance_diffusion(mla,n_old,n_new,ext_src_d,dx,0.5d0*dt,the_bc_tower) ! swap n_new/n_old to avoid calling copy()

       ! restore boundary conditions
       n_bc(1:dm,1:2,1:nspecies) = n_bc_temp(1:dm,1:2,1:nspecies)

       ! (later) recompute z with new-time boundary conditions
       !
       !

       ! n_new = n_new + z
       do n=1,nlevs
          call multifab_plus_plus_c(n_new(n),1,z(n),1,nspecies,0)
          call multifab_fill_boundary(n_new(n))
          call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                               the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
       end do

       do n=1,nlevs
          call multifab_destroy(z(n))
          call multifab_destroy(zerofab(n))
       end do

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
