module diffusive_fluxdiv_module

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
  use electrodiffusive_mass_fluxdiv_module
  use probin_common_module, only: variance_coef_mass, nspecies
  use probin_charged_module, only: use_charged_fluid

  implicit none

  private

  public :: diffusive_fluxdiv

contains

  ! compute diffusive mass and rhoh fluxes
  ! also compute eta and kappa
  subroutine diffusive_fluxdiv(mla,rho,rhotot,gradp_baro,Temp,p0,eta,eta_ed,kappa,
                               diff_mass_fluxdiv,diff_mass_flux,diff_rhoh_fluxdiv, &
                               dx,the_bc_tower)
       
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: rho(:)
    type(multifab) , intent(in   ) :: rhotot(:)
    type(multifab) , intent(in   ) :: gradp_baro(:,:)
    type(multifab) , intent(in   ) :: Temp(:)
    type(multifab) , intent(inout) :: eta(:)
    type(multifab) , intent(inout) :: eta_ed(:)
    type(multifab) , intent(inout) :: kappa(:)
    type(multifab) , intent(inout) :: diff_mass_fluxdiv(:)
    type(multifab) , intent(inout) :: diff_mass_flux(:,:)
    type(multifab) , intent(inout) :: diff_rhoh_fluxdiv(:,:)
    real(kind=dp_t), intent(in   ) :: p0
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    type(multifab) :: conc(mla%nlevel)
    type(multifab) :: molarconc(mla%nlevel)
    type(multifab) :: chi(mla%nlevel)
    type(multifab) :: lambda(mla%nlevel)
    type(multifab) :: rhoWchi(mla%nlevel)        ! rho*W*chi*Gama
    type(multifab) :: Gama(mla%nlevel)           ! Gama-matrix
    type(multifab) :: zeta_over_Temp(mla%nlevel)   ! for Thermo-diffusion 

    integer         :: n,i,dm,nlevs

    type(bl_prof_timer), save :: bpt

    call build(bpt,"diffusive_fluxdiv")

    nlevs = mla%nlevel  ! number of levels 
    dm    = mla%dim     ! dimensionality
      
    do n=1,nlevs
       call multifab_build(conc(n)          , mla%la(n), nspecies   , 1)
       call multifab_build(molarconc(n)     , mla%la(n), nspecies   , 1)
       call multifab_build(chi(n)           , mla%la(n), nspecies**2, 1)
       call multifab_build(lambda(n)        , mla%la(n), 1          , 1)
       call multifab_build(rhoWchi(n)       , mla%la(n), nspecies**2, 1)
       call multifab_build(Gama(n)          , mla%la(n), nspecies**2, 1)
       call multifab_build(zeta_over_Temp(n), mla%la(n), nspecies   , 1)
    end do

    ! compute molmtot, molarconc (primitive variables) for 
    ! each-cell from rho(conserved) 
    call compute_molconc_molmtot(mla,rho,rhotot,molarconc,molmtot)

    ! rho to conc - VALID REGION ONLY
    call convert_rhoc_to_c(mla,rho,rhotot,conc,.true.)

    ! fill conc ghost cells
    call fill_c_ghost_cells(mla,conc,dx,the_bc_tower)

    ! compute transport properties (eta,lambda,kappa,chi,zeta)
    call ideal_mixture_transport_wrapper(mla,rhotot,Temp,p0,conc, &
                                         molarconc,eta,lambda,kappa, &
                                         chi,zeta_over_Temp)

    ! put eta on edges
    if (dm .eq. 2) then
       call average_cc_to_node(nlevs,eta,eta_ed(:,1),1,tran_bc_comp,1,the_bc_level)
    else if (dm .eq. 3) then
       call average_cc_to_edge(nlevs,eta,eta_ed,1,tran_bc_comp,1,the_bc_level)
    end if

    ! set zeta_over_Temp to zeta/Temp
    do n=1,nlevs
       call multifab_copy_c(zeta_over_Temp(n),1,zeta(n),1,nspecies,1)
       do i=1,nspecies
          call multifab_div_div_c(zeta_over_Temp(n),i,Temp(n),1,1,1)
       end do
    end do

    ! set Gama to the identity matrix
    do n=1,nlevs
       call multifab_setval(Gama(n),0.d0,all=.true.)
       j=1
       do i=1,dm
          call multifab_setval_c(Gama(n),1.d0,j,1,all=.true.)
          j=j+dm+1
       end do
    end do

    ! compute rho*W*chi
    call compute_rhoWchi_from_chi(mla,rho,chi,rhoWchi)

    ! compute diffusive mass fluxes, "-F = rho*W*chi*Gamma*grad(x) - ..."
    call diffusive_mass_fluxdiv(mla,rho,rhotot,molarconc,rhoWchi,Gama, &
                                diff_mass_fluxdiv,Temp,zeta_over_Temp,gradp_baro, &
                                diff_mass_flux,dx,the_bc_tower)

    ! compute diffusive rhoh fluxes
    call diffusive_rhoh_fluxdiv(mla,lambda,Temp,diff_mass_flux,rhotot,diff_rhoh_fluxdiv, &
                                dx,time,the_bc_tower)

    do n=1,nlevs
       call multifab_destroy(conc(n))
       call multifab_destroy(molarconc(n))
       call multifab_destroy(chi(n))
       call multifab_destroy(lambda(n))
       call multifab_destroy(rhoWchi(n))
       call multifab_destroy(Gama(n))
       call multifab_destroy(zeta_over_Temp(n))
    end do

    call destroy(bpt)

  end subroutine diffusive_fluxdiv
  
end module diffusive_fluxdiv_module
