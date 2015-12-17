module advance_timestep_module

  use ml_layout_module
  use define_bc_module
  use advance_diffusion_module
  use advance_reaction_module
  use advance_reaction_diffusion_module
  use chemical_rates_module
  use multifab_physbc_module
  use bc_module
  use probin_reactdiff_module, only: nspecies, temporal_integrator, n_bc, reaction_type, &
                                     use_Poisson_rng, inhomogeneous_bc_fix

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

    type(multifab) :: Rn_steady(mla%nlevel)

    type(bl_prof_timer),save :: bpt

    call build(bpt,"advance_timestep")

    nlevs = mla%nlevel
    dm = mla%dim

    ! external source term for diffusion/reaction solvers for inhomogeneous bc algorithm
    do n=1,nlevs
       call multifab_build(Rn_steady(n),mla%la(n),nspecies,0)
       call multifab_setval(Rn_steady(n),0.d0,all=.true.)
    end do


    if (inhomogeneous_bc_fix) then

       if (temporal_integrator .lt. 0) then
          call bl_error("inhomogeneous_bc_fix not supported for temporal_integrator < 0")
       end if

       ! store reactions rates for n_steady in Rn_steady
       ! the input time step does not matter as the reaction_type/use_Poisson_rng settings are
       ! returning an explicit rate in units of number_density/time
       call deterministic_chemical_rates(mla,n_steady,Rn_steady,dx,dt)

    end if


    if (temporal_integrator .lt. 0) then  ! unsplitting schmes

       call advance_reaction_diffusion(mla,n_old,n_new,dx,dt,the_bc_tower)
   
    else if (temporal_integrator .eq. 0) then  ! D + R

       call advance_diffusion(mla,n_old,n_new,dx,dt,the_bc_tower,ext_src_in=Rn_steady)
       call advance_reaction (mla,n_new,n_old,dx,dt,the_bc_tower,ext_src_in=Rn_steady)  ! swap n_new/n_old to avoid calling copy()
       do n=1,nlevs
          call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,n_new(n)%ng)  ! make sure n_new contains the new state
       end do

    else if (temporal_integrator .eq. 1) then  ! (1/2)R + D + (1/2)R

       call advance_reaction (mla,n_old,n_new,dx,0.5d0*dt,the_bc_tower,ext_src_in=Rn_steady)
       call advance_diffusion(mla,n_new,n_old,dx,dt      ,the_bc_tower,ext_src_in=Rn_steady)  ! swap n_new/n_old to avoid calling copy()
       call advance_reaction (mla,n_old,n_new,dx,0.5d0*dt,the_bc_tower,ext_src_in=Rn_steady)  ! swap n_new/n_old to avoid calling copy()

    else if (temporal_integrator .eq. 2) then  ! (1/2)D + R + (1/2)D

       call advance_diffusion(mla,n_old,n_new,dx,0.5d0*dt,the_bc_tower,ext_src_in=Rn_steady)
       call advance_reaction (mla,n_new,n_old,dx,dt      ,the_bc_tower,ext_src_in=Rn_steady)  ! swap n_new/n_old to avoid calling copy()
       call advance_diffusion(mla,n_old,n_new,dx,0.5d0*dt,the_bc_tower,ext_src_in=Rn_steady)  ! swap n_new/n_old to avoid calling copy()

    else
       call bl_error("advance_timestep: invalid temporal_integrator")
    end if

              
    do n=1,nlevs
       call multifab_destroy(Rn_steady(n))
    end do

    call destroy(bpt)

  end subroutine advance_timestep

end module advance_timestep_module
