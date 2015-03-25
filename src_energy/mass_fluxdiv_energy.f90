module mass_fluxdiv_energy_module

  use ml_layout_module
  use define_bc_module
  use diffusive_mass_fluxdiv_module
  use mass_flux_utilities_module
  use convert_variables_module
  use energy_EOS_module
  use probin_multispecies_module, only: nspecies
  use probin_common_module, only: variance_coef_mass

  implicit none

  private

  public :: mass_fluxdiv_energy

contains

  subroutine mass_fluxdiv_energy(mla,rho,rhotot,molefrac,chi,zeta,gradp_baro,mass_fluxdiv, &
                                 Temp,dx,the_bc_tower)
       
    type(ml_layout), intent(in   )   :: mla
    type(multifab) , intent(in   )   :: rho(:)
    type(multifab) , intent(in   )   :: rhotot(:)
    type(multifab) , intent(in   )   :: molefrac(:)
    type(multifab) , intent(in   )   :: chi(:)
    type(multifab) , intent(in   )   :: zeta(:)
    type(multifab) , intent(in   )   :: gradp_baro(:,:)
    type(multifab) , intent(inout)   :: mass_fluxdiv(:)
    type(multifab) , intent(in   )   :: Temp(:)
    real(kind=dp_t), intent(in   )   :: dx(:,:)
    type(bc_tower) , intent(in   )   :: the_bc_tower

    ! local variables
    type(multifab) :: rhoWchi(mla%nlevel)        ! rho*W*chi*Gama
    type(multifab) :: Gama(mla%nlevel)           ! Gama-matrix
    type(multifab) :: zeta_by_Temp(mla%nlevel)   ! for Thermo-diffusion 
    type(multifab) :: flux_total(mla%nlevel,mla%dim)

    integer         :: n,i,j,dm,nlevs

    nlevs = mla%nlevel  ! number of levels 
    dm    = mla%dim     ! dimensionality
      
    ! build cell-centered multifabs for nspecies and ghost cells contained in rho.
    do n=1,nlevs
       call multifab_build(rhoWchi(n),      mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(Gama(n),         mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(zeta_by_Temp(n), mla%la(n), nspecies,    rho(n)%ng)
       do i=1,dm
          call multifab_build_edge(flux_total(n,i),mla%la(n),nspecies,0,i)
          call setval(flux_total(n,i),0.d0,all=.true.)
       end do
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
    call compute_rhoWchi(mla,rho,rhotot,chi,rhoWchi)

    ! compute determinstic mass fluxdiv (interior only), rho contains ghost filled 
    ! in init/end of this code
    call diffusive_mass_fluxdiv(mla,rho,rhotot,molefrac,rhoWchi,Gama, &
                                mass_fluxdiv,Temp,zeta_by_Temp,gradp_baro,flux_total,dx, &
                                the_bc_tower)
      
    ! free the multifab allocated memory
    do n=1,nlevs
       call multifab_destroy(rhoWchi(n))
       call multifab_destroy(Gama(n))
       call multifab_destroy(zeta_by_Temp(n))
       do i=1,dm
          call multifab_destroy(flux_total(n,i))
       end do
    end do

  end subroutine mass_fluxdiv_energy
  
end module mass_fluxdiv_energy_module
