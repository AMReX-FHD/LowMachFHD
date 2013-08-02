module advance_module

  use multifab_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use div_and_grad_module
  use diffusive_flux_module
  use ml_layout_module
  use probin_multispecies_module, only: nspecies

  implicit none

  private

  public :: advance

contains

  ! Donev: molarconc should not be passed to advance, it should be a local multifab
  ! EVERY time you change rho you must recalculate molconc again -- where is this done?
  subroutine advance(mla,rho,molarconc,BinvGama,Dbar,Gama,dx,dt,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rho(:)
    type(multifab) , intent(inout) :: molarconc(:)
    type(multifab) , intent(inout) :: BinvGama(:)
    real(kind=dp_t), intent(in   ) :: Dbar(:,:)
    real(kind=dp_t), intent(in   ) :: Gama(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(in   ) :: dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local variables
    integer i, dm, n, nlevs

    ! an array of multifabs; one for each direction
    type(multifab) :: flux(mla%nlevel,mla%dim)
    type(multifab) :: fluxdiv(mla%nlevel)
 
    dm    = mla%dim     ! dimensionality
    nlevs = mla%nlevel  ! number of levels 

    ! build the flux and div-of-flux multifabs
    do n=1,nlevs
       ! fluxdiv is cell-centered scalar with zero ghost cells 
       call multifab_build(fluxdiv(n),mla%la(n),nspecies,0)
       do i=1,dm
          ! flux(i) is face-centered, has nspecies component, zero ghost cells & nodal in direction i
          call multifab_build_edge(flux(n,i),mla%la(n),nspecies,0,i)
       end do
    end do   
    
    ! compute the face-centered flux in each direction
    call diffusive_flux(mla,rho,molarconc,BinvGama,Dbar,Gama,flux,dx,the_bc_level)
    
    ! compute divergence of the flux 
    call compute_div(mla,flux,fluxdiv,dx,1,1,nspecies)

    ! update rho using forward Euler discretization
    call update_rho(mla,rho,fluxdiv,dt,the_bc_level)

    ! destroy the multifab to prevent leakage in memory
    do n=1,nlevs
       call multifab_destroy(fluxdiv(n))
       do i=1,dm
          call multifab_destroy(flux(n,i))
       end do
    end do

  end subroutine advance

  subroutine update_rho(mla,rho,fluxdiv,dt,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rho(:)
    type(multifab) , intent(inout) :: fluxdiv(mla%nlevel)
    real(kind=dp_t), intent(in   ) :: dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local variables
    integer n, nlevs

    nlevs = mla%nlevel  ! number of levels 
    
    ! Euler explicit time update
    do n=1,nlevs
       ! multiply div-of-flux with dt, starting from component 1 with 0 ghost-cell
       call multifab_mult_mult_s_c(fluxdiv(n),1,dt,nspecies,0) 
    end do 
    
    do n=1,nlevs
       ! add this to rho to advance in time
       call multifab_plus_plus(rho(n),fluxdiv(n))            
    end do 

    do n=1, nlevs
       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(rho(n))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(rho(n),1,scal_bc_comp,nspecies,the_bc_level(n))
    end do   

  end subroutine update_rho

end module advance_module
