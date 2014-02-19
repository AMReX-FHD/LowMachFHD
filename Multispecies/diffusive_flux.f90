module diffusive_flux_module

  use multifab_module
  use define_bc_module
  use div_and_grad_module
  use probin_multispecies_module
  use ml_layout_module
  use convert_stag_module
  use matvec_mul_module
  use matmat_mul_module
  use F95_LAPACK
  
  implicit none

  private

  public :: diffusive_flux

contains
 
  subroutine diffusive_flux(mla,rho,rho_tot,molarconc,rhoWchiGama,flux,dx,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: rho(:) 
    type(multifab) , intent(in   ) :: rho_tot(:) 
    type(multifab) , intent(in   ) :: molarconc(:) 
    type(multifab) , intent(in   ) :: rhoWchiGama(:)  
    type(multifab) , intent(inout) :: flux(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local variables
    integer :: n,i,dm,ng,nlevs
 
    ! local face-centered multifabs 
    type(multifab)  :: rhoWchiGama_face(mla%nlevel,mla%dim)
   
    dm    = mla%dim     ! dimensionality
    nlevs = mla%nlevel  ! number of levels 

    ! build local face-centered multifab with nspecies^2 component, zero ghost cells 
    ! and nodal in direction i
    do n=1,nlevs
       do i=1,dm
          call multifab_build_edge(rhoWchiGama_face(n,i),mla%la(n),nspecies**2,0,i)
       enddo
    enddo 

    ! calculate face-centrered grad(molarconc) 
    call compute_grad(mla, molarconc, flux, dx, 1, mol_frac_bc_comp, 1, nspecies, & 
                      the_bc_level, .false.)
   
    ! compute face-centered rhoWchiGama from cell-centered values 
    call average_cc_to_face(nlevs, rhoWchiGama, rhoWchiGama_face, 1, diff_coeff_bc_comp, &
                            nspecies**2, the_bc_level, .false.) 
    
    ! compute flux as B^(-1)*Gama X grad(molarconc). 
    do n=1,nlevs
       do i=1,dm
          call matvec_mul(mla, flux(n,i), rhoWchiGama_face(n,i))
       enddo
    enddo    
     
    ! destroy B^(-1)*Gama multifab to prevent leakage in memory
    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(rhoWchiGama_face(n,i))
       enddo
    enddo

  end subroutine diffusive_flux

end module diffusive_flux_module
