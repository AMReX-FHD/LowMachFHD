module advance_timestep_module

  use ml_layout_module
  use define_bc_module
  use advance_diffusion_module
  use advance_reaction_module
  use advance_euler_module
  use advance_explicit_midpoint_module
  use multifab_physbc_module
  use bc_module
  use probin_reactdiff_module, only: nspecies, splitting_type, n_bc, reaction_type, &
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

       ! store reactions rates for n_steady in Rn_steady
       ! the input time step does not matter as the reaction_type/use_Poisson_rng settings are
       ! returning an explicit rate in units of number_density/time
       call advance_reaction(mla,n_steady,Rn_steady,dx,dt,the_bc_tower,return_rates_in=.true.)

    end if

    select case(splitting_type)
    case(0)
       ! D + R

       call advance_diffusion(mla,n_old,n_new,dx,dt,the_bc_tower,ext_src_in=Rn_steady)
       call advance_reaction (mla,n_new,n_old,dx,dt,the_bc_tower,ext_src_in=Rn_steady)  ! swap n_new/n_old to avoid calling copy()
       do n=1,nlevs
          call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,n_new(n)%ng) ! make sure n_new contains the new state
       end do

    case(1)
       ! (1/2)R + D + (1/2)R

       call advance_reaction (mla,n_old,n_new,dx,0.5d0*dt,the_bc_tower,ext_src_in=Rn_steady)
       call advance_diffusion(mla,n_new,n_old,dx,dt      ,the_bc_tower,ext_src_in=Rn_steady) ! swap n_new/n_old to avoid calling copy()
       call advance_reaction (mla,n_old,n_new,dx,0.5d0*dt,the_bc_tower,ext_src_in=Rn_steady) ! swap n_new/n_old to avoid calling copy()

    case(2)
       ! (1/2)D + R + (1/2)D

       call advance_diffusion(mla,n_old,n_new,dx,0.5d0*dt,the_bc_tower,ext_src_in=Rn_steady)
       call advance_reaction (mla,n_new,n_old,dx,dt      ,the_bc_tower,ext_src_in=Rn_steady) ! swap n_new/n_old to avoid calling copy()
       call advance_diffusion(mla,n_old,n_new,dx,0.5d0*dt,the_bc_tower,ext_src_in=Rn_steady) ! swap n_new/n_old to avoid calling copy()

    case(3)
       ! no splitting
       ! euler scheme

       call advance_euler(mla,n_old,n_new,dx,dt,the_bc_tower)

    case(4)
       ! no splitting
       ! explicit midpoint predictor-corrector scheme

       call advance_explicit_midpoint(mla,n_old,n_new,dx,dt,the_bc_tower)

    case default
       call bl_error("advance_timestep: invalid splitting_type")
    end select
              
    do n=1,nlevs
       call multifab_destroy(Rn_steady(n))
    end do

    call destroy(bpt)

  end subroutine advance_timestep

end module advance_timestep_module
