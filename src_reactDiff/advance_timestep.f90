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

  subroutine advance_timestep(mla,n_old,n_new,n_steady,dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: n_old(:)
    type(multifab) , intent(inout) :: n_new(:)
    type(multifab) , intent(in   ) :: n_steady(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: n,nlevs,dm

    type(multifab) :: fz(mla%nlevel)

    type(bl_prof_timer),save :: bpt

    call build(bpt,"advance_timestep")

    nlevs = mla%nlevel
    dm = mla%dim

    select case(splitting_type)
    case(0)
       ! D + R

       call advance_diffusion(mla,n_old,n_new,dx,dt,the_bc_tower)
       call advance_reaction (mla,n_new,n_old,dx,dt,the_bc_tower)  ! swap n_new/n_old to avoid calling copy()
       do n=1,nlevs
          call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,n_new(n)%ng) ! make sure n_new contains the new state
       end do

    case(1)
       ! (1/2)R + D + (1/2)R

       call advance_reaction (mla,n_old,n_new,dx,0.5d0*dt,the_bc_tower)
       call advance_diffusion(mla,n_new,n_old,dx,dt      ,the_bc_tower) ! swap n_new/n_old to avoid calling copy()
       call advance_reaction (mla,n_old,n_new,dx,0.5d0*dt,the_bc_tower) ! swap n_new/n_old to avoid calling copy()

    case(2)
       ! (1/2)D + R + (1/2)D

       call advance_diffusion(mla,n_old,n_new,dx,0.5d0*dt,the_bc_tower)
       call advance_reaction (mla,n_new,n_old,dx,dt      ,the_bc_tower) ! swap n_new/n_old to avoid calling copy()
       call advance_diffusion(mla,n_old,n_new,dx,0.5d0*dt,the_bc_tower) ! swap n_new/n_old to avoid calling copy()

    case(3)
       ! (1/2)D + R + (1/2)D with inhomogeneous, time-independent boundary conditions
       ! under development

       ! external source term for diffusion/reaction solvers for inhomogeneous bc algorithm
       do n=1,nlevs
          call multifab_build(fz(n),mla%la(n),nspecies,0)
       end do
       
       ! store reactions rates for n_steady in fz
       ! the input time step does not matter as the reaction_type/use_Poisson_rng settings are
       ! returning an explicit rate in units of number_density/time
       call advance_reaction(mla,n_steady,fz,dx,dt,the_bc_tower,return_rates_in=.true.)

       ! This code should be identical to case=2 just passing in nonzero external sources:
       call advance_diffusion(mla,n_old,n_new,dx,0.5d0*dt,the_bc_tower,ext_src_in=fz)
       call advance_reaction (mla,n_new,n_old,dx,dt      ,the_bc_tower,ext_src_in=fz) ! swap n_new/n_old to avoid calling copy()
       call advance_diffusion(mla,n_old,n_new,dx,0.5d0*dt,the_bc_tower,ext_src_in=fz) ! swap n_new/n_old to avoid calling copy()
              
       do n=1,nlevs
          call multifab_destroy(fz(n))
       end do

    case default
       call bl_error("advance_timestep: invalid splitting_type")
    end select

    call destroy(bpt)

  end subroutine advance_timestep

end module advance_timestep_module
