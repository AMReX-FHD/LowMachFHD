module compute_mass_fluxdiv_module

  use multifab_module
  use define_bc_module
  use bc_module
  use div_and_grad_module
  use diffusive_mass_fluxdiv_module
  use fluid_model_module
  use external_force_module
  use ml_layout_module
  use F95_LAPACK
  use stochastic_mass_fluxdiv_module
  use convert_mass_variables_module
  use probin_multispecies_module, only: nspecies, use_stoch

  implicit none

  private

  public :: compute_mass_fluxdiv_wrapper, compute_mass_fluxdiv

contains

  subroutine compute_mass_fluxdiv_wrapper(mla,rho,rho_tot, &
                                          diff_fluxdiv,stoch_fluxdiv,Temp,flux_total, &
                                          dt,stage_time,dx,weights, &
                                          n_rngs,the_bc_level)
       
    type(ml_layout), intent(in   )   :: mla
    type(multifab) , intent(inout)   :: rho(:)
    type(multifab) , intent(inout)   :: rho_tot(:)
    type(multifab) , intent(inout)   :: diff_fluxdiv(:)
    type(multifab) , intent(inout)   :: stoch_fluxdiv(:)
    type(multifab) , intent(in   )   :: Temp(:)
    type(multifab) , intent(inout)   :: flux_total(:,:)
    real(kind=dp_t), intent(in   )   :: dt
    real(kind=dp_t), intent(in   )   :: stage_time 
    real(kind=dp_t), intent(in   )   :: dx(:,:)
    real(kind=dp_t), intent(in   )   :: weights(:) 
    integer,         intent(in   )   :: n_rngs
    type(bc_level) , intent(in   )   :: the_bc_level(:)

    ! local
    integer :: n,nlevs

    type(multifab) :: molarconc(mla%nlevel)      ! molar concentration
    type(multifab) :: molmtot(mla%nlevel)        ! total molar mass
    type(multifab) :: chi(mla%nlevel)            ! Chi-matrix
    type(multifab) :: Gama(mla%nlevel)           ! Gama-matrix
    type(multifab) :: D_MS(mla%nlevel)           ! D_MS-matrix
    type(multifab) :: D_therm(mla%nlevel)        ! DT-matrix
    type(multifab) :: zeta_by_Temp(mla%nlevel)   ! for Thermo-diffusion 

    nlevs = mla%nlevel

    do n=1,nlevs
       call multifab_build(molarconc(n),    mla%la(n), nspecies,    rho(n)%ng)
       call multifab_build(molmtot(n),      mla%la(n), 1,           rho(n)%ng)
       call multifab_build(chi(n),          mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(Gama(n),         mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(D_MS(n),         mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(D_therm(n),      mla%la(n), nspecies,    rho(n)%ng)
       call multifab_build(zeta_by_Temp(n), mla%la(n), nspecies,    rho(n)%ng)
    end do

    call compute_mass_fluxdiv(mla,rho,rho_tot,molarconc,molmtot,chi,Gama,D_MS,&
                              D_therm,diff_fluxdiv,stoch_fluxdiv,Temp,&
                              zeta_by_Temp,flux_total,dt,stage_time,dx,weights,&
                              n_rngs,the_bc_level)

    do n=1,nlevs
       call multifab_destroy(molarconc(n))
       call multifab_destroy(molmtot(n))
       call multifab_destroy(chi(n))
       call multifab_destroy(Gama(n))
       call multifab_destroy(D_MS(n))
       call multifab_destroy(D_therm(n))
       call multifab_destroy(zeta_by_Temp(n))
    end do

  end subroutine compute_mass_fluxdiv_wrapper

  subroutine compute_mass_fluxdiv(mla,rho,rho_tot,molarconc,molmtot,chi,Gama,D_MS,&
                                  D_therm,diff_fluxdiv,stoch_fluxdiv,Temp,&
                                  zeta_by_Temp,flux_total,dt,stage_time,dx,weights,&
                                  n_rngs,the_bc_level)
       
    type(ml_layout), intent(in   )   :: mla
    type(multifab) , intent(inout)   :: rho(:)
    type(multifab) , intent(inout)   :: rho_tot(:)
    type(multifab) , intent(inout)   :: molarconc(:)
    type(multifab) , intent(inout)   :: molmtot(:)
    type(multifab) , intent(inout)   :: chi(:)
    type(multifab) , intent(inout)   :: Gama(:)
    type(multifab) , intent(inout)   :: D_MS(:)
    type(multifab) , intent(inout)   :: D_therm(:)
    type(multifab) , intent(inout)   :: diff_fluxdiv(:)
    type(multifab) , intent(inout)   :: stoch_fluxdiv(:)
    type(multifab) , intent(in   )   :: Temp(:)
    type(multifab) , intent(inout)   :: zeta_by_Temp(:)
    type(multifab) , intent(inout)   :: flux_total(:,:)
    real(kind=dp_t), intent(in   )   :: dt
    real(kind=dp_t), intent(in   )   :: stage_time 
    real(kind=dp_t), intent(in   )   :: dx(:,:)
    real(kind=dp_t), intent(in   )   :: weights(:) 
    integer,         intent(in   )   :: n_rngs
    type(bc_level) , intent(in   )   :: the_bc_level(:)

    ! local variables
    type(multifab)  :: drho(mla%nlevel)  ! correction to rho
    type(multifab)  :: rhoWchi(mla%nlevel)    ! rho*W*chi*Gama

    integer         :: n,i,dm,nlevs

    nlevs = mla%nlevel  ! number of levels 
    dm    = mla%dim     ! dimensionality
      
    ! build cell-centered multifabs for nspecies and ghost cells contained in rho.
    do n=1,nlevs
       call multifab_build(drho(n),    mla%la(n), nspecies,    rho(n)%ng)
       call multifab_build(rhoWchi(n), mla%la(n), nspecies**2, rho(n)%ng)
    end do
 
    ! modify rho with drho to ensure no mass or mole fraction is zero
    call correct_rho_with_drho(mla,rho,drho,the_bc_level)
 
    ! compute molmtot,molarconc & rho_tot (primitive variables) for 
    ! each-cell from rho(conserved) 
    call convert_cons_to_prim(mla,rho,rho_tot,molarconc,molmtot,the_bc_level)
      
    ! populate D_MS and Gama 
    call fluid_model(mla,rho,rho_tot,molarconc,molmtot,D_MS,D_therm,Gama,the_bc_level)

    ! compute chi 
    call compute_chi(mla,rho,rho_tot,molarconc,chi,D_MS,D_therm,Temp,zeta_by_Temp,the_bc_level)
      
    ! compute rho*W*chi
    call compute_rhoWchi(mla,rho,rho_tot,molarconc,molmtot,chi,rhoWchi,the_bc_level)

    ! compute determinstic mass fluxdiv (interior only), rho contains ghost filled 
    ! in init/end of this code
    call diffusive_mass_fluxdiv(mla,rho,rho_tot,molarconc,rhoWchi,Gama,&
                                diff_fluxdiv,Temp,zeta_by_Temp,flux_total,dx,the_bc_level)

    ! compute external forcing for manufactured solution and add to diff_fluxdiv
    call external_source(mla,rho,diff_fluxdiv,dx,stage_time)

    ! compute stochastic fluxdiv 
    if(use_stoch) call stochastic_mass_fluxdiv(mla,rho,rho_tot,molarconc,&
                                               molmtot,chi,Gama,stoch_fluxdiv,flux_total,&
                                               dx,dt,weights,the_bc_level)
      
    ! revert back rho to it's original form
    do n=1,nlevs
       call saxpy(rho(n),-1.0d0,drho(n))
    end do 
      
    ! free the multifab allocated memory
    do n=1,nlevs
       call multifab_destroy(drho(n))
       call multifab_destroy(rhoWchi(n))
    end do

  end subroutine compute_mass_fluxdiv
  
end module compute_mass_fluxdiv_module
