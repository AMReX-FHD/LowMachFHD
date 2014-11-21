module baro_mass_fluxdiv_module

  use ml_layout_module
  use multifab_module
  use define_bc_module
  use bc_module
  use convert_stag_module
  use div_and_grad_module
  use probin_multispecies_module, only: nspecies
  implicit none

  private

  public :: baro_mass_fluxdiv

contains

  subroutine baro_mass_fluxdiv(mla,baro_coef,gradp_baro,baro_fluxdiv,flux_total,dx,the_bc_level)

    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: baro_coef(:)
    type(multifab) , intent(in   )  :: gradp_baro(:,:)
    type(multifab) , intent(inout)  :: baro_fluxdiv(:)
    type(multifab) , intent(inout)  :: flux_total(:,:)
    real(kind=dp_t), intent(in   )  :: dx(:,:)
    type(bc_level) , intent(in   )  :: the_bc_level(:)

    ! local variables
    integer i,dm,n,nlevs,comp

    ! local array of multifabs for grad and div; one for each direction
    type(multifab) :: flux(mla%nlevel,mla%dim)
    
    nlevs = mla%nlevel  ! number of levels 
    dm    = mla%dim     ! dimensionality
 
    ! build the local multifabs
    do n=1,nlevs
       do i=1,dm
          ! flux(i) is face-centered, has nspecies component, zero ghost 
          ! cells & nodal in direction i
          call multifab_build_edge(flux(n,i),mla%la(n),nspecies,0,i)
       end do
    end do

    ! average baro_coef to faces, store in 'flux'
    call average_cc_to_face(nlevs, baro_coef, flux, 1, scal_bc_comp, &
                            nspecies, the_bc_level, .false.) 

    ! flux(1:nspecies) = baro_coef(1:nspecies) * gradp_baro
    do n=1,nlevs
       do i=1,dm
          do comp=1,nspecies
             call multifab_mult_mult_c(flux(n,i),comp,gradp_baro(n,i),1,1,0)
          end do
       end do
    end do
    
    ! add fluxes to flux_total
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(flux_total(n,i),1,flux(n,i),1,nspecies,0)
       end do
    end do

    ! compute divergence of flux
    call compute_div(mla,flux,baro_fluxdiv,dx,1,1,nspecies)
    
    ! destroy the multifab to free the memory
    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(flux(n,i))
       end do
    end do

  end subroutine baro_mass_fluxdiv

end module baro_mass_fluxdiv_module
