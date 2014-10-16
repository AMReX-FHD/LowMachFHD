module compute_mass_fluxdiv_module

  use multifab_module
  use define_bc_module
  use bc_module
  use div_and_grad_module
  use diffusive_mass_fluxdiv_module
  use compute_mixture_properties_module
  use external_force_module
  use ml_layout_module
  use F95_LAPACK
  use stochastic_mass_fluxdiv_module
  use mass_flux_utilities_module
  use probin_multispecies_module, only: nspecies
  use probin_common_module, only: variance_coef_mass

  implicit none

  private

  public :: compute_mass_fluxdiv_wrapper, compute_mass_fluxdiv

contains

  subroutine compute_mass_fluxdiv_wrapper(mla,rho, &
                                          diff_fluxdiv,stoch_fluxdiv,Temp,flux_total, &
                                          dt,stage_time,dx,weights, &
                                          the_bc_level)
       
    type(ml_layout), intent(in   )   :: mla
    type(multifab) , intent(inout)   :: rho(:)
    type(multifab) , intent(inout)   :: diff_fluxdiv(:)
    type(multifab) , intent(inout)   :: stoch_fluxdiv(:)
    type(multifab) , intent(in   )   :: Temp(:)
    type(multifab) , intent(inout)   :: flux_total(:,:)
    real(kind=dp_t), intent(in   )   :: dt
    real(kind=dp_t), intent(in   )   :: stage_time 
    real(kind=dp_t), intent(in   )   :: dx(:,:)
    real(kind=dp_t), intent(in   )   :: weights(:) 
    type(bc_level) , intent(in   )   :: the_bc_level(:)

    ! local
    integer :: n,nlevs

    type(multifab) :: molarconc(mla%nlevel)      ! molar concentration
    type(multifab) :: molmtot(mla%nlevel)        ! total molar mass
    type(multifab) :: chi(mla%nlevel)            ! Chi-matrix
    type(multifab) :: Hessian(mla%nlevel)        ! Hessian-matrix
    type(multifab) :: Gama(mla%nlevel)           ! Gama-matrix
    type(multifab) :: D_bar(mla%nlevel)          ! D_bar-matrix
    type(multifab) :: D_therm(mla%nlevel)        ! DT-matrix
    type(multifab) :: zeta_by_Temp(mla%nlevel)   ! for Thermo-diffusion 

    nlevs = mla%nlevel

    do n=1,nlevs
       call multifab_build(molarconc(n),    mla%la(n), nspecies,    rho(n)%ng)
       call multifab_build(molmtot(n),      mla%la(n), 1,           rho(n)%ng)
       call multifab_build(chi(n),          mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(Hessian(n),      mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(Gama(n),         mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(D_bar(n),        mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(D_therm(n),      mla%la(n), nspecies,    rho(n)%ng)
       call multifab_build(zeta_by_Temp(n), mla%la(n), nspecies,    rho(n)%ng)
    end do

    call compute_mass_fluxdiv(mla,rho,molarconc,molmtot,chi,Hessian,Gama,D_bar,&
                              D_therm,diff_fluxdiv,stoch_fluxdiv,Temp,&
                              zeta_by_Temp,flux_total,dt,stage_time,dx,weights,&
                              the_bc_level)

    do n=1,nlevs
       call multifab_destroy(molarconc(n))
       call multifab_destroy(molmtot(n))
       call multifab_destroy(chi(n))
       call multifab_destroy(Hessian(n))
       call multifab_destroy(Gama(n))
       call multifab_destroy(D_bar(n))
       call multifab_destroy(D_therm(n))
       call multifab_destroy(zeta_by_Temp(n))
    end do

  end subroutine compute_mass_fluxdiv_wrapper

  subroutine compute_mass_fluxdiv(mla,rho,molarconc,molmtot,chi,Hessian,Gama,D_bar,&
                                  D_therm,diff_fluxdiv,stoch_fluxdiv,Temp,&
                                  zeta_by_Temp,flux_total,dt,stage_time,dx,weights,&
                                  the_bc_level)
       
    type(ml_layout), intent(in   )   :: mla
    type(multifab) , intent(inout)   :: rho(:)
    type(multifab) , intent(inout)   :: molarconc(:)
    type(multifab) , intent(inout)   :: molmtot(:)
    type(multifab) , intent(inout)   :: chi(:)
    type(multifab) , intent(inout)   :: Hessian(:)
    type(multifab) , intent(inout)   :: Gama(:)
    type(multifab) , intent(inout)   :: D_bar(:)
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
    type(bc_level) , intent(in   )   :: the_bc_level(:)

    ! local variables
    type(multifab)  :: drho(mla%nlevel)  ! correction to rho
    type(multifab)  :: rhoWchi(mla%nlevel)    ! rho*W*chi*Gama
    type(multifab)  :: rhotot_temp(mla%nlevel)

    integer         :: n,i,dm,nlevs

    nlevs = mla%nlevel  ! number of levels 
    dm    = mla%dim     ! dimensionality
      
    ! build cell-centered multifabs for nspecies and ghost cells contained in rho.
    do n=1,nlevs
       call multifab_build(drho(n),       mla%la(n), nspecies,    rho(n)%ng)
       call multifab_build(rhoWchi(n),    mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(rhotot_temp(n),mla%la(n), 1          , rho(n)%ng)
    end do
 
    ! modify rho with drho to ensure no mass or mole fraction is zero
    call correct_rho_with_drho(mla,rho,drho,the_bc_level)
 
    ! compute molmtot,molarconc & rhotot_temp (primitive variables) for 
    ! each-cell from rho(conserved) 
    call convert_cons_to_prim(mla,rho,rhotot_temp,molarconc,molmtot,the_bc_level)
      
    ! populate D_bar and Hessian matrix 
    call compute_mixture_properties(mla,rho,rhotot_temp,molarconc,molmtot,D_bar,D_therm, &
                                    Hessian,Temp,the_bc_level)

    ! compute Gama from Hessian
    call compute_Gama(mla,rho,rhotot_temp,molarconc,molmtot,Hessian,Gama,the_bc_level)
   
    ! compute chi 
    call compute_chi(mla,rho,rhotot_temp,molarconc,chi,D_bar,D_therm,Temp,zeta_by_Temp,the_bc_level)
      
    ! compute rho*W*chi
    call compute_rhoWchi(mla,rho,rhotot_temp,molarconc,molmtot,chi,rhoWchi,the_bc_level)

    ! reset total flux
    do n=1,nlevs
       do i=1,dm
          call setval(flux_total(n,i),0.d0,all=.true.)
       end do
    end do

    ! compute determinstic mass fluxdiv (interior only), rho contains ghost filled 
    ! in init/end of this code
    call diffusive_mass_fluxdiv(mla,rho,rhotot_temp,molarconc,rhoWchi,Gama,&
                                diff_fluxdiv,Temp,zeta_by_Temp,flux_total,dx,the_bc_level)

    ! compute external forcing for manufactured solution and add to diff_fluxdiv
    call external_source(mla,rho,diff_fluxdiv,dx,stage_time)

    ! compute stochastic fluxdiv 
    if (variance_coef_mass .ne. 0.d0) then
       call stochastic_mass_fluxdiv(mla,rho,rhotot_temp,molarconc,&
                                    molmtot,chi,Gama,stoch_fluxdiv,flux_total,&
                                    dx,dt,weights,the_bc_level)
    else
       do n=1,nlevs
          call multifab_setval(stoch_fluxdiv(n),0.d0,all=.true.)
       end do
    end if
      
    ! revert back rho to it's original form
    do n=1,nlevs
       call saxpy(rho(n),-1.0d0,drho(n))
    end do 
      
    ! free the multifab allocated memory
    do n=1,nlevs
       call multifab_destroy(drho(n))
       call multifab_destroy(rhoWchi(n))
       call multifab_destroy(rhotot_temp(n))
    end do

  end subroutine compute_mass_fluxdiv
  
end module compute_mass_fluxdiv_module
