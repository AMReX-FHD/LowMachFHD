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
  use probin_common_module, only: variance_coef_mass, nspecies

  implicit none

  private

  public :: compute_mass_fluxdiv

contains

  ! compute diffusive and stochastic mass fluxes
  ! includes barodiffusion and thermodiffusion
  subroutine compute_mass_fluxdiv(mla,rho,rhotot,gradp_baro,diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                  Temp,total_mass_flux,dt,stage_time,dx,weights, &
                                  the_bc_tower,diff_mass_flux,rhoWchi_out)
       
    type(ml_layout), intent(in   )   :: mla
    type(multifab) , intent(inout)   :: rho(:)
    type(multifab) , intent(inout)   :: rhotot(:)
    type(multifab) , intent(in   )   :: gradp_baro(:,:)
    type(multifab) , intent(inout)   :: diff_mass_fluxdiv(:)
    type(multifab) , intent(inout)   :: stoch_mass_fluxdiv(:)
    type(multifab) , intent(in   )   :: Temp(:)
    type(multifab) , intent(inout)   :: total_mass_flux(:,:)
    real(kind=dp_t), intent(in   )   :: dt
    real(kind=dp_t), intent(in   )   :: stage_time 
    real(kind=dp_t), intent(in   )   :: dx(:,:)
    real(kind=dp_t), intent(in   )   :: weights(:) 
    type(bc_tower) , intent(in   )   :: the_bc_tower
    type(multifab) , intent(inout), optional :: diff_mass_flux(:,:)
    type(multifab) , intent(inout), optional :: rhoWchi_out(:)

    ! local variables
    type(multifab) :: drho(mla%nlevel)           ! correction to rho
    type(multifab) :: rhoWchi(mla%nlevel)        ! rho*W*chi*Gama
    type(multifab) :: molarconc(mla%nlevel)      ! molar concentration
    type(multifab) :: molmtot(mla%nlevel)        ! total molar mass
    type(multifab) :: chi(mla%nlevel)            ! Chi-matrix
    type(multifab) :: Hessian(mla%nlevel)        ! Hessian-matrix
    type(multifab) :: Gama(mla%nlevel)           ! Gama-matrix
    type(multifab) :: D_bar(mla%nlevel)          ! D_bar-matrix
    type(multifab) :: D_therm(mla%nlevel)        ! DT-matrix
    type(multifab) :: zeta_by_Temp(mla%nlevel)   ! for Thermo-diffusion 
    type(multifab) :: sqrtLonsager_fc(mla%nlevel,mla%dim) ! cholesky factored Lonsager on face

    integer         :: n,i,dm,nlevs

    type(bl_prof_timer), save :: bpt

    call build(bpt,"compute_mass_fluxdiv")

    nlevs = mla%nlevel  ! number of levels 
    dm    = mla%dim     ! dimensionality
      
    ! build cell-centered multifabs for nspecies and ghost cells contained in rho.
    do n=1,nlevs
       call multifab_build(drho(n),         mla%la(n), nspecies,    rho(n)%ng)
       call multifab_build(rhoWchi(n),      mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(molarconc(n),    mla%la(n), nspecies,    rho(n)%ng)
       call multifab_build(molmtot(n),      mla%la(n), 1,           rho(n)%ng)
       call multifab_build(chi(n),          mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(Hessian(n),      mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(Gama(n),         mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(D_bar(n),        mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(D_therm(n),      mla%la(n), nspecies,    rho(n)%ng)
       call multifab_build(zeta_by_Temp(n), mla%la(n), nspecies,    rho(n)%ng)
       do i=1,dm
          call multifab_build_edge(sqrtLonsager_fc(n,i), mla%la(n), nspecies**2, 0, i)
       end do
    end do
 
    ! modify rho with drho to ensure no mass or mole fraction is zero
    call correct_rho_with_drho(mla,rho,drho)
    call compute_rhotot(mla,rho,rhotot,ghost_cells_in=.true.)
 
    ! compute molmtot, molarconc (primitive variables) for 
    ! each-cell from rho(conserved) 
    call compute_molconc_molmtot(mla,rho,rhotot,molarconc,molmtot)
      
    ! populate D_bar and Hessian matrix 
    call compute_mixture_properties(mla,rho,rhotot,D_bar,D_therm,Hessian)

    ! compute Gama from Hessian
    call compute_Gama(mla,molarconc,Hessian,Gama)
   
    ! compute chi and zeta/Temp
    call compute_chi(mla,rho,rhotot,molarconc,chi,D_bar)
    call compute_zeta_by_Temp(mla,molarconc,D_bar,D_therm,Temp,zeta_by_Temp)

    ! compute rho*W*chi
    call compute_rhoWchi(mla,rho,chi,rhoWchi)

    ! pass out a copy of rhoWchi (needed for the electric potential code)
    if (present(rhoWchi_out)) then
       do n=1,nlevs
          call multifab_copy_c(rhoWchi_out(n),1,rhoWchi(n),1,nspecies**2,rhoWchi_out(n)%ng)
       end do
    end if

    ! reset total flux
    do n=1,nlevs
       do i=1,dm
          call setval(total_mass_flux(n,i),0.d0,all=.true.)
       end do
    end do

    ! compute diffusive mass fluxes, "F = -rho*W*chi*Gamma*grad(x) - ..."
    call diffusive_mass_fluxdiv(mla,rho,rhotot,molarconc,rhoWchi,Gama, &
                                diff_mass_fluxdiv,Temp,zeta_by_Temp,gradp_baro, &
                                total_mass_flux,dx,the_bc_tower)

    ! we need to save only the diffusive fluxes for the implicit potential algorithms
    if (present(diff_mass_flux)) then
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(diff_mass_flux(n,i),1,total_mass_flux(n,i),1,nspecies,0)
          end do
       end do
    end if

    ! compute external forcing for manufactured solution and add to diff_mass_fluxdiv
    call external_source(mla,rho,diff_mass_fluxdiv,dx,stage_time)

    ! revert back rho to it's original form
    do n=1,nlevs
       call saxpy(rho(n),-1.0d0,drho(n),all=.true.)
       call compute_rhotot(mla,rho,rhotot,ghost_cells_in=.true.)
    end do 

    ! compute stochastic fluxdiv 
    if (variance_coef_mass .ne. 0.d0) then

       ! compute face-centered cholesky-factored Lonsager^(1/2)
       call compute_sqrtLonsager_fc(mla,rho,rhotot,sqrtLonsager_fc,dx)

       call stochastic_mass_fluxdiv(mla,rho,rhotot, &
                                    sqrtLonsager_fc,stoch_mass_fluxdiv,total_mass_flux,&
                                    dx,dt,weights,the_bc_tower%bc_tower_array)
    else
       do n=1,nlevs
          call multifab_setval(stoch_mass_fluxdiv(n),0.d0,all=.true.)
       end do
    end if
      
    ! free the multifab allocated memory
    do n=1,nlevs
       call multifab_destroy(drho(n))
       call multifab_destroy(rhoWchi(n))
       call multifab_destroy(molarconc(n))
       call multifab_destroy(molmtot(n))
       call multifab_destroy(chi(n))
       call multifab_destroy(Hessian(n))
       call multifab_destroy(Gama(n))
       call multifab_destroy(D_bar(n))
       call multifab_destroy(D_therm(n))
       call multifab_destroy(zeta_by_Temp(n))
       do i=1,dm
          call multifab_destroy(sqrtLonsager_fc(n,i))
       end do
    end do

    call destroy(bpt)

  end subroutine compute_mass_fluxdiv
  
end module compute_mass_fluxdiv_module
