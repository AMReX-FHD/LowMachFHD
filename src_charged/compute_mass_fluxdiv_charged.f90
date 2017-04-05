module compute_mass_fluxdiv_charged_module

  use multifab_module
  use define_bc_module
  use bc_module
  use div_and_grad_module
  use diffusive_mass_fluxdiv_module
  use stochastic_mass_fluxdiv_module
  use compute_mixture_properties_module
  use external_force_module
  use ml_layout_module
  use mass_flux_utilities_module
  use convert_stag_module
  use compute_mass_fluxdiv_module
  use Epot_mass_fluxdiv_module
  use probin_common_module, only: nspecies
  use probin_charged_module, only: use_charged_fluid

  implicit none

  private

  public :: compute_mass_fluxdiv_charged

contains

  ! compute diffusive, stochastic, and electric potential mass fluxes
  ! includes barodiffusion and thermodiffusion
  subroutine compute_mass_fluxdiv_charged(mla,rho,rhotot,gradp_baro,diff_fluxdiv,stoch_fluxdiv, &
                                          Temp,total_mass_flux,dt,stage_time,dx,weights,the_bc_tower, &
                                          Epot_fluxdiv,charge,grad_Epot,Epot,permittivity, &
                                          flux_diff) ! Donev: Add a logical flag for whether to do electroneutral or not
       
    type(ml_layout), intent(in   )   :: mla
    type(multifab) , intent(inout)   :: rho(:)
    type(multifab) , intent(inout)   :: rhotot(:)
    type(multifab) , intent(in   )   :: gradp_baro(:,:)
    type(multifab) , intent(inout)   :: diff_fluxdiv(:)
    type(multifab) , intent(inout)   :: stoch_fluxdiv(:)
    type(multifab) , intent(in   )   :: Temp(:)
    type(multifab) , intent(inout)   :: total_mass_flux(:,:)
    real(kind=dp_t), intent(in   )   :: dt
    real(kind=dp_t), intent(in   )   :: stage_time 
    real(kind=dp_t), intent(in   )   :: dx(:,:)
    real(kind=dp_t), intent(in   )   :: weights(:) 
    type(bc_tower) , intent(in   )   :: the_bc_tower
    type(multifab) , intent(inout)   :: Epot_fluxdiv(:)
    type(multifab) , intent(inout)   :: charge(:)
    type(multifab) , intent(inout)   :: grad_Epot(:,:)
    type(multifab) , intent(inout)   :: Epot(:)
    type(multifab) , intent(in   )   :: permittivity(:)
    type(multifab) , intent(inout), optional :: flux_diff(:,:)
       
    ! local variables
    type(multifab) :: rhoWchi(mla%nlevel)        ! rho*W*chi*Gama

    integer         :: n,nlevs

    type(bl_prof_timer), save :: bpt

    call build(bpt,"compute_mass_fluxdiv_charged")
    
    ! Donev: I propose the following rewrite:
    ! call compute_mass_fluxdiv() ! Compute F=F_bar+F_tilde using existing routine
    ! if(electroneutral) then
    !    solve Poisson equation with epsilon=0 and compute Epot to return to caller
    !    note no advective fluxes required in this case
    !    project fluxes by adding div(A_Phi grad Phi)
    ! else
    !    call Epot_mass_fluxdiv() ! Add the electrostatic piece using Epot passed in
    !    may be better to compute Epot here by solving the simple Poisson equation though
    !    this way callers don't have to worry about Poisson solves. 
    !    but one needs to check if this will work with all of the existing algorithms
    ! end

    nlevs = mla%nlevel  ! number of levels 
      
    ! build cell-centered multifabs for nspecies and ghost cells contained in rho.
    do n=1,nlevs
       call multifab_build(rhoWchi(n), mla%la(n), nspecies**2, rho(n)%ng)
    end do

    call compute_mass_fluxdiv(mla,rho,rhotot,gradp_baro,diff_fluxdiv,stoch_fluxdiv, &
                              Temp,total_mass_flux,dt,stage_time,dx,weights,the_bc_tower, &
                              flux_diff,rhoWchi)

    ! Donev: Observe that modified densities rho+drho are used when computing charges. Is this what we want?
    ! In particular, in the implicit method, do we want to make sure there is charge everywhere ad-hoc? Probably not
    if (use_charged_fluid) then
       ! compute electric potential mass fluxes
       call Epot_mass_fluxdiv(mla,rho,Epot_fluxdiv,Temp,rhoWchi, &
                              total_mass_flux,dx,the_bc_tower,charge,grad_Epot,Epot, &
                              permittivity)
    end if
      
    ! free the multifab allocated memory
    do n=1,nlevs
       call multifab_destroy(rhoWchi(n))
    end do

    call destroy(bpt)

  end subroutine compute_mass_fluxdiv_charged
  
end module compute_mass_fluxdiv_charged_module
