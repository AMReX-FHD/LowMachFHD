module diffusive_fluxdiv_module

  use multifab_module
  use define_bc_module
  use bc_module
  use div_and_grad_module
  use diffusive_flux_module
  use fluid_model_module
  use external_force_module
  use ml_layout_module
  use F95_LAPACK
  use stochastic_fluxdiv_module
  use convert_variables_module
  use probin_multispecies_module

  implicit none

  private

  public :: diffusive_fluxdiv, compute_fluxdiv

contains

  subroutine diffusive_fluxdiv(mla,rho,rho_tot,diff_fluxdiv,molarconc,rhoWchiGama,molmass,dx,the_bc_level)

    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: rho(:)
    type(multifab) , intent(in   )  :: rho_tot(:)
    type(multifab) , intent(inout)  :: diff_fluxdiv(:)
    type(multifab) , intent(inout)  :: molarconc(:)
    type(multifab) , intent(inout)  :: rhoWchiGama(:)
    real(kind=dp_t), intent(in   )  :: molmass(:) 
    real(kind=dp_t), intent(in   )  :: dx(:,:)
    type(bc_level) , intent(in   )  :: the_bc_level(:)

    ! local variables
    integer i,dm,n,nlevs

    ! local array of multifabs for grad and div; one for each direction
    type(multifab) :: flux(mla%nlevel,mla%dim)
    
    nlevs = mla%nlevel  ! number of levels 
    dm    = mla%dim     ! dimensionality
 
    ! build the local multifabs
    do n=1,nlevs
       do i=1,dm
          ! flux(i) is face-centered, has nspecies component, zero ghost cells & nodal in direction i
          call multifab_build_edge(flux(n,i),mla%la(n),nspecies,0,i)
       end do
    end do   
    
    ! compute the face-centered flux (each direction: cells+1 faces while cells contain interior+2 ghost cells) 
    call diffusive_flux(mla,rho,rho_tot,molarconc,rhoWchiGama,flux,dx,the_bc_level)

    ! compute divergence of determinstic flux 
    call compute_div(mla,flux,diff_fluxdiv,dx,1,1,nspecies)
    
    ! destroy the multifab to free the memory
    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(flux(n,i))
       end do
    end do

  end subroutine diffusive_fluxdiv

  subroutine compute_fluxdiv(mla,rho,rho_tot,molarconc,molmtot,molmass,chi,Gama,D_MS,&
                             rhoWchiGama,diff_fluxdiv,stoch_fluxdiv,stoch_W_fc,dt,&
                             stage_time,dx,prob_lo,prob_hi,weights,n_rngs,the_bc_level)
       
    type(ml_layout), intent(in   )   :: mla
    type(multifab) , intent(inout)   :: rho(:)
    type(multifab) , intent(inout)   :: rho_tot(:)
    type(multifab) , intent(inout)   :: molarconc(:)
    type(multifab) , intent(inout)   :: molmtot(:)
    real(kind=dp_t), intent(in   )   :: molmass(nspecies) 
    type(multifab) , intent(inout)   :: chi(:)
    type(multifab) , intent(inout)   :: Gama(:)
    type(multifab) , intent(inout)   :: D_MS(:)
    type(multifab) , intent(inout)   :: rhoWchiGama(:)
    type(multifab) , intent(inout)   :: diff_fluxdiv(:)
    type(multifab) , intent(inout)   :: stoch_fluxdiv(:)
    type(multifab) , intent(inout)   :: stoch_W_fc(:,:,:)
    real(kind=dp_t), intent(in   )   :: dt
    real(kind=dp_t), intent(in   )   :: stage_time 
    real(kind=dp_t), intent(in   )   :: dx(:,:)
    real(kind=dp_t), intent(in   )   :: prob_lo(rho(1)%dim)
    real(kind=dp_t), intent(in   )   :: prob_hi(rho(1)%dim) 
    real(kind=dp_t), intent(in   )   :: weights(:) 
    integer,         intent(in   )   :: n_rngs
    type(bc_level) , intent(in   )   :: the_bc_level(:)

    ! local variables
    type(multifab)  :: drho(mla%nlevel)  ! correction to rho
    integer         :: n,i,dm,nlevs

    nlevs = mla%nlevel  ! number of levels 
    dm    = mla%dim     ! dimensionality
      
    ! build cell-centered multifabs for nspecies and ghost cells contained in rho.
    do n=1,nlevs
       call multifab_build(drho(n),mla%la(n),nspecies,rho(n)%ng)
    enddo
 
    ! modify rho with drho to ensure no mass or mole fraction is zero
    call correct_rho_with_drho(mla,rho,drho,the_bc_level)
 
    ! compute molmtot,molarconc & rho_tot (primitive variables) for each-cell from rho(conserved) 
    call convert_cons_to_prim(mla,rho,rho_tot,molarconc,molmtot,molmass,the_bc_level)
      
    ! populate D_MS and Gama 
    call fluid_model(mla,rho,rho_tot,molarconc,molmtot,D_MS,Gama,the_bc_level)

    ! compute chi 
    call compute_chi(mla,rho,rho_tot,molarconc,molmass,chi,D_MS,the_bc_level)
      
    ! compute rho*W*chi*Gama
    call compute_rhoWchiGama(mla,rho,rho_tot,molarconc,molmass,molmtot,chi,&
                             Gama,rhoWchiGama,the_bc_level)

    ! compute determinstic fluxdiv (interior only), rho contains ghost filled in init/end of this code
    call diffusive_fluxdiv(mla,rho,rho_tot,diff_fluxdiv,molarconc,rhoWchiGama,&
                           molmass,dx,the_bc_level)

    ! compute external forcing for manufactured solution and add to diff_fluxdiv
    call external_source(mla,rho,diff_fluxdiv,prob_lo,prob_hi,dx,stage_time)

    ! compute stochastic fluxdiv 
    if(use_stoch) call stochastic_fluxdiv(mla,rho,rho_tot,molarconc,molmass,molmtot,chi,&
                                          Gama,stoch_W_fc,stoch_fluxdiv,dx,dt,weights,&
                                          the_bc_level)
      
    ! revert back rho to it's original form
    do n=1,nlevs
       call saxpy(rho(n),-1.0d0,drho(n))
    enddo 
      
    ! free the multifab allocated memory
    do n=1,nlevs
       call multifab_destroy(drho(n))
    enddo

  end subroutine compute_fluxdiv
  
end module diffusive_fluxdiv_module
