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

  subroutine diffusive_fluxdiv(mla,rho,rho_tot,fluxdiv,molarconc,molmtot,BinvGamma,chi,&
                               Dbar,Gama,mass,dx,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: rho(:)
    type(multifab) , intent(in)    :: rho_tot(:)
    type(multifab) , intent(inout) :: fluxdiv(:)
    type(multifab) , intent(inout) :: molarconc(:)
    type(multifab) , intent(inout) :: molmtot(:)
    type(multifab) , intent(inout) :: BinvGamma(:)
    type(multifab) , intent(inout) :: chi(:)
    real(kind=dp_t), intent(in   ) :: Dbar(:,:)
    real(kind=dp_t), intent(in   ) :: Gama(:,:)
    real(kind=dp_t), intent(in   ) :: mass(:) 
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local variables
    integer i,dm,n,nlevs

    ! local array of multifabs for grad and div; one for each direction
    type(multifab) :: flux(mla%nlevel,mla%dim)
    
    dm    = mla%dim     ! dimensionality
    nlevs = mla%nlevel  ! number of levels 
 
    ! build the local multifabs
    do n=1,nlevs
       do i=1,dm
          ! flux(i) is face-centered, has nspecies component, zero ghost cells & nodal in direction i
          call multifab_build_edge(flux(n,i),mla%la(n),nspecies,0,i)
       end do
    end do   
    
    ! compute cell-centered B^(-1)*Gamma  
    call compute_BinvGamma(mla,rho,rho_tot,molarconc,BinvGamma,Dbar,Gama,mass,molmtot,the_bc_level)
 
    ! compute the face-centered flux (each direction: cells+1 faces while cells contain interior+2 ghost cells) 
    call diffusive_flux(mla,rho,molarconc,BinvGamma,chi,Gama,flux,dx,the_bc_level)

    ! compute divergence of the flux 
    call compute_div(mla,flux,fluxdiv,dx,1,1,nspecies)
    
    ! multiply fluxdiv (having zero ghost cells) with -1 to get -div(-flux).
    do n=1,nlevs
       call multifab_mult_mult_s(fluxdiv(n),-1.0d0,fluxdiv(1)%ng)
    end do
 
    ! destroy the multifab to free the memory
    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(flux(n,i))
       end do
    end do

  end subroutine diffusive_fluxdiv
  
end module diffusive_fluxdiv_module
