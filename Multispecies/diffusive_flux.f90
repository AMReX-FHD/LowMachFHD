module diffusive_flux_module

  use multifab_module
  use define_bc_module
  use bc_module
  use div_and_grad_module
  use probin_multispecies_module
  use ml_layout_module
  use convert_stag_module
  use matvec_mul_module
  use F95_LAPACK
  
  implicit none

  private

  public :: diffusive_flux

contains
 
  subroutine diffusive_flux(mla, molarconc, BinvGamma, flux, dx, the_bc_level, & 
                            mol_frac_bc_comp, diff_coeff_bc_comp) 

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: molarconc(:) 
    type(multifab) , intent(in   ) :: BinvGamma(:)  
    type(multifab) , intent(inout) :: flux(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    integer        , intent(in   ) :: mol_frac_bc_comp
    integer        , intent(in   ) :: diff_coeff_bc_comp 

    ! local variables
    integer :: n,i,dm,nlevs
 
    ! local multifab for the face-centered B^(-1)*Gama
    type(multifab) :: BinvGamma_face(mla%nlevel,mla%dim)
   
    dm    = mla%dim     ! dimensionality
    nlevs = mla%nlevel  ! number of levels 

    ! build local face-centered multifab with nspecies^2 component, zero ghost cells 
    ! and nodal in direction i
    do n=1,nlevs
       do i=1,dm
          call multifab_build_edge(BinvGamma_face(n,i), mla%la(n), nspecies**2, 0, i)
       end do
    end do 

    ! calculate face-centrered grad(molarconc) (CHECKED) 
    call compute_grad(mla, molarconc, flux, dx, 1, mol_frac_bc_comp, 1, nspecies, & 
                      the_bc_level, .false.)
   
    ! compute face-centered B^(-1)*Gamma from cell-centered values 
    call average_cc_to_face(nlevs, BinvGamma, BinvGamma_face, 1, diff_coeff_bc_comp, & 
                            nspecies**2, the_bc_level, .false.) 
    
    ! compute flux as B^(-1)*Gama X grad(molarconc). 
    do n=1,nlevs
       do i=1,dm
          call matvec_mul(mla, flux(n,i), BinvGamma_face(n,i))
       end do
    end do    
    
    !Donev: If grad(temperature)
    !call compute_grad(mla,temperature,flux,dx,1,scal_bc_comp+nspecies,nspecies+1,1,the_bc_level)
 
    ! destroy B^(-1)*Gama multifab to prevent leakage in memory
    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(BinvGamma_face(n,i))
       end do
    end do

  end subroutine diffusive_flux

end module diffusive_flux_module
