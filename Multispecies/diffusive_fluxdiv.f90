module diffusive_fluxdiv_module

  use multifab_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use div_and_grad_module
  use diffusive_flux_module
  use ml_layout_module
  use convert_variables_module
  use probin_multispecies_module

  implicit none

  private

  public :: diffusive_fluxdiv

contains

  subroutine diffusive_fluxdiv(mla,rho,fluxdiv,molmtot,Dbar,Gama,mass,dx,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: rho(:)
    type(multifab) , intent(inout) :: fluxdiv(:)
    type(multifab) , intent(inout) :: molmtot(:)
    real(kind=dp_t), intent(in   ) :: Dbar(:,:)
    real(kind=dp_t), intent(in   ) :: Gama(:,:)
    real(kind=dp_t), intent(in   ) :: mass(:) 
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local variables
    integer i,dm,n,nlevs

    ! local array of multifabs for grad and div; one for each direction
    type(multifab) :: flux(mla%nlevel,mla%dim)
    
    ! local array of multifabs for total density, molarconc and BinvGamma 
    ! Donev: The following three arrays should be moved outside this routine (in advance)
    ! since they will be used in several places, including the fluctuating fluxes
    ! For first testing this is OK but later it needs to change
    type(multifab) :: rho_tot(mla%nlevel)
    type(multifab) :: molarconc(mla%nlevel)
    type(multifab) :: BinvGamma(mla%nlevel)
 
    dm    = mla%dim     ! dimensionality
    nlevs = mla%nlevel  ! number of levels 
 
    ! build the local multifabs
    do n=1,nlevs
       ! fluxdiv,rho_tot,molarconc are cell-cented with ghost cells that rho has
       call multifab_build(rho_tot(n),  mla%la(n),1,rho(n)%ng)  ! rho_tot is addition of all component
       call multifab_build(molarconc(n),mla%la(n),nspecies,rho(n)%ng) 
       call multifab_build(BinvGamma(n),mla%la(n),nspecies**2,rho(n)%ng)
       do i=1,dm
          ! flux(i) is face-centered, has nspecies component, zero ghost cells & nodal in direction i
          call multifab_build_edge(flux(n,i),mla%la(n),nspecies,0,i)
       end do
    end do   
    
    ! compute molarconc (primary) and rho_tot (primary) for every cell from rho(1:nspecies) 
    ! Donev: When the code is finalized this call should be outside, in advance, not here
    call convert_cons_to_prim(mla,rho,rho_tot,molarconc,mass,molmtot,the_bc_level)

    ! compute cell-centered B^(-1)*Gamma  
    call compute_BinvGamma(mla,rho,rho_tot,molarconc,BinvGamma,Dbar,Gama,mass,molmtot,the_bc_level)
 
    ! compute the face-centered flux (each direction: cells+1 faces while cells
    ! contain: interior + 2 ghost cells) 
    call diffusive_flux(mla,molarconc,BinvGamma,flux,dx,the_bc_level)

    ! compute divergence of the flux 
    call compute_div(mla,flux,fluxdiv,dx,1,1,nspecies)
    
    ! multiply fluxdiv (having zero ghost cells) with -1 to get -div(-flux).
    do n=1,nlevs
       call multifab_mult_mult_s(fluxdiv(n),-1.0d0,fluxdiv(1)%ng)
    end do
 
    ! destroy the multifab to prevent leakage in memory
    do n=1,nlevs
       call multifab_destroy(rho_tot(n))
       call multifab_destroy(molarconc(n))
       call multifab_destroy(BinvGamma(n))
       do i=1,dm
          call multifab_destroy(flux(n,i))
       end do
    end do

  end subroutine diffusive_fluxdiv
  
end module diffusive_fluxdiv_module
