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

    type(multifab) :: fz(mla%nlevel)
    type(multifab) :: z(mla%nlevel)

    type(bl_prof_timer),save :: bpt

    call build(bpt,"advance_timestep")

    nlevs = mla%nlevel
    dm = mla%dim

    ! external source term for diffusion/reaction solvers for inhomogeneous bc algorithm
    ! Donev: It seems to me these should only be allocated if needed (splitting_algorithm=3 etc.)
    do n=1,nlevs
       call multifab_build(fz(n),mla%la(n),nspecies,0)
       call setval(fz(n),0.d0,all=.true.)
    end do

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

       do n=1,nlevs
          call multifab_build(z(n),mla%la(n),nspecies,0)
       end do

       ! compute z with old-time boundary conditions
       ! Donev: Actually z should be computed only once in main.f90 or constructed analytically to have a simple gradint
       ! In general the user will know how to solve div D_k grad z_k = 0 manually...
       ! We definitely do NOT want to be solving a Poisson problem every time step -- that is much more expensive than a whole time step of react-diff      

       do n=1,nlevs
          ! temporary hack to test n_bc=1 case
          ! in general we will solve div D_k grad z_k = 0
          call setval(z(n),1.d0,all=.true.)
       end do
       
       ! store reactions rates for z in fz
       ! the input time step does not matter as the reaction_type/use_Poisson_rng settings are
       ! returning an explicit rate in units of number_density/time
       call advance_reaction(mla,z,fz,dx,dt,the_bc_tower,return_rates_in=.true.)

       ! This code should be identical to case=2 just passing in nonzero external sources:
       call advance_diffusion(mla,n_old,n_new,dx,0.5d0*dt,the_bc_tower,ext_src_in=fz)
       call advance_reaction (mla,n_new,n_old,dx,dt      ,the_bc_tower,ext_src_in=fz) ! swap n_new/n_old to avoid calling copy()
       call advance_diffusion(mla,n_old,n_new,dx,0.5d0*dt,the_bc_tower,ext_src_in=fz) ! swap n_new/n_old to avoid calling copy()
              
       do n=1,nlevs
          call multifab_destroy(z(n))
       end do

    case default
       call bl_error("advance_timestep: invalid splitting_type")
    end select

    do n=1,nlevs
       call multifab_destroy(fz(n))
    end do


    call destroy(bpt)

  end subroutine advance_timestep

end module advance_timestep_module
