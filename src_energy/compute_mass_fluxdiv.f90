module compute_mass_fluxdiv_module

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

  public :: compute_mass_fluxdiv

contains

  ! compute diffusive mass fluxes
  subroutine compute_mass_fluxdiv(mla,rho,rhotot,molarconc,chi,zeta,gradp_baro,Temp, &
                                  diff_mass_fluxdiv,diff_mass_flux,dx,the_bc_tower)
       
    type(ml_layout), intent(in   )   :: mla
    type(multifab) , intent(in   )   :: rho(:)
    type(multifab) , intent(in   )   :: rhotot(:)
    type(multifab) , intent(in   )   :: molarconc(:)
    type(multifab) , intent(in   )   :: gradp_baro(:,:)
    type(multifab) , intent(in   )   :: Temp(:)
    type(multifab) , intent(inout)   :: diff_mass_fluxdiv(:)
    type(multifab) , intent(inout)   :: diff_mass_flux(:,:)
    real(kind=dp_t), intent(in   )   :: dx(:,:)
    type(bc_tower) , intent(in   )   :: the_bc_tower

    ! local variables
    type(multifab) :: rhoWchi(mla%nlevel)        ! rho*W*chi*Gama
    type(multifab) :: Gama(mla%nlevel)           ! Gama-matrix
    type(multifab) :: zeta_by_Temp(mla%nlevel)   ! for Thermo-diffusion 

    integer         :: n,i,dm,nlevs

    type(bl_prof_timer), save :: bpt

    call build(bpt,"compute_mass_fluxdiv")

    nlevs = mla%nlevel  ! number of levels 
    dm    = mla%dim     ! dimensionality
      
    do n=1,nlevs
       call multifab_build(rhoWchi(n),      mla%la(n), nspecies**2, 1)
       call multifab_build(Gama(n),         mla%la(n), nspecies**2, 1)
       call multifab_build(zeta_by_Temp(n), mla%la(n), nspecies,    1)
    end do

    ! set zeta_by_Temp to zeta/Temp
    do n=1,nlevs
       call multifab_copy_c(zeta_by_Temp(n),1,zeta(n),1,nspecies,1)
       do i=1,nspecies
          call multifab_div_div_c(zeta_by_Temp(n),i,Temp(n),1,1,1)
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
    call compute_rhoWchi(mla,rho,chi,rhoWchi)

    do n=1,nlevs
       do i=1,dm
          call setval(mass_flux(n,i),0.d0,all=.true.)
       end do
    end do


    ! compute diffusive mass fluxes, "-F = rho*W*chi*Gamma*grad(x) - ..."
    call diffusive_mass_fluxdiv(mla,rho,rhotot,molarconc,rhoWchi,Gama, &
                                diff_mass_fluxdiv,Temp,zeta_by_Temp,gradp_baro, &
                                diff_mass_flux,dx,the_bc_tower)

    do n=1,nlevs
       call multifab_destroy(rhoWchi(n))
       call multifab_destroy(Gama(n))
       call multifab_destroy(zeta_by_Temp(n))
    end do

    call destroy(bpt)

  end subroutine compute_mass_fluxdiv
  
end module compute_mass_fluxdiv_module
