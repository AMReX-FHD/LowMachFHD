module mass_fluxdiv_energy_module

  use ml_layout_module
  use define_bc_module
  use diffusive_mass_fluxdiv_module
  use mass_flux_utilities_module
  use convert_variables_module
  use probin_multispecies_module, only: nspecies
  use probin_common_module, only: variance_coef_mass

  implicit none

  private

  public :: mass_fluxdiv_energy

contains

  subroutine mass_fluxdiv_energy(mla,rho,gradp_baro, &
                                  mass_fluxdiv, &
                                  Temp, &
                                  dt,stage_time,dx,weights, &
                                  the_bc_tower)
       
    type(ml_layout), intent(in   )   :: mla
    type(multifab) , intent(inout)   :: rho(:)
    type(multifab) , intent(in   )   :: gradp_baro(:,:)
    type(multifab) , intent(inout)   :: mass_fluxdiv(:)
    type(multifab) , intent(in   )   :: Temp(:)
    real(kind=dp_t), intent(in   )   :: dt
    real(kind=dp_t), intent(in   )   :: stage_time 
    real(kind=dp_t), intent(in   )   :: dx(:,:)
    real(kind=dp_t), intent(in   )   :: weights(:) 
    type(bc_tower) , intent(in   )   :: the_bc_tower

    ! local variables
    type(multifab) :: drho(mla%nlevel)           ! correction to rho
    type(multifab) :: rhoWchi(mla%nlevel)        ! rho*W*chi*Gama
    type(multifab) :: rhotot_temp(mla%nlevel)    ! temp storage for rho with drho correction
    type(multifab) :: massfrac(mla%nlevel)       ! mass fractions (concentrations)
    type(multifab) :: molarconc(mla%nlevel)      ! molar concentration
    type(multifab) :: molmtot(mla%nlevel)        ! total molar mass
    type(multifab) :: chi(mla%nlevel)            ! Chi-matrix
    type(multifab) :: Gama(mla%nlevel)           ! Gama-matrix
    type(multifab) :: zeta_by_Temp(mla%nlevel)   ! for Thermo-diffusion 
    type(multifab) :: flux_total(mla%nlevel,mla%dim)

    integer         :: n,i,dm,nlevs

    nlevs = mla%nlevel  ! number of levels 
    dm    = mla%dim     ! dimensionality
      
    ! build cell-centered multifabs for nspecies and ghost cells contained in rho.
    do n=1,nlevs
       call multifab_build(drho(n),         mla%la(n), nspecies,    rho(n)%ng)
       call multifab_build(rhoWchi(n),      mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(rhotot_temp(n),  mla%la(n), 1          , rho(n)%ng)
       call multifab_build(massfrac(n),     mla%la(n), nspecies,    rho(n)%ng)
       call multifab_build(molarconc(n),    mla%la(n), nspecies,    rho(n)%ng)
       call multifab_build(molmtot(n),      mla%la(n), 1,           rho(n)%ng)
       call multifab_build(chi(n),          mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(Gama(n),         mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(zeta_by_Temp(n), mla%la(n), nspecies,    rho(n)%ng)
       do i=1,dm
          call multifab_build_edge(flux_total(n,i),mla%la(n),nspecies,0,i)
          call setval(flux_total(n,i),0.d0,all=.true.)
       end do
    end do
 
    ! modify rho with drho to ensure no mass or mole fraction is zero
    call correct_rho_with_drho(mla,rho,drho)

    ! compute rhotot_temp from corrected rho
    call compute_rhotot(mla,rho,rhotot_temp)

    ! compute mass fractions 
    call convert_rho_to_c(mla,rho,rhotot_temp,massfrac,.true.)

    ! compute mole fractions and molar mass
    call convert_cons_to_prim(mla,rho,rhotot_temp,molarconc,molmtot)
 
    ! use ideal_mixture_transport()
    ! inputs are rho, Temp, P, Y, X
    ! outputs are viscosity (eta), bulk viscosity (kappa), diffusion matrix (chi), 
    ! thermal conductivity (lambda), scaled thermodiffusion coefficients

      

    ! set Gama to the identity matrix


    ! compute rho*W*chi
    call compute_rhoWchi(mla,rho,rhotot_temp,chi,rhoWchi)

    ! reset total flux
    do n=1,nlevs
       do i=1,dm
          call setval(flux_total(n,i),0.d0,all=.true.)
       end do
    end do

    ! compute determinstic mass fluxdiv (interior only), rho contains ghost filled 
    ! in init/end of this code
    call diffusive_mass_fluxdiv(mla,rho,rhotot_temp,molarconc,rhoWchi,Gama, &
                                mass_fluxdiv,Temp,zeta_by_Temp,gradp_baro,flux_total,dx, &
                                the_bc_tower)

    ! revert back rho to it's original form
    do n=1,nlevs
       call saxpy(rho(n),-1.0d0,drho(n))
    end do 
      
    ! free the multifab allocated memory
    do n=1,nlevs
       call multifab_destroy(drho(n))
       call multifab_destroy(rhoWchi(n))
       call multifab_destroy(rhotot_temp(n))
       call multifab_destroy(massfrac(n))
       call multifab_destroy(molarconc(n))
       call multifab_destroy(molmtot(n))
       call multifab_destroy(chi(n))
       call multifab_destroy(Gama(n))
       call multifab_destroy(zeta_by_Temp(n))
       do i=1,dm
          call multifab_destroy(flux_total(n,i))
       end do
    end do

  end subroutine mass_fluxdiv_energy
  
end module mass_fluxdiv_energy_module
