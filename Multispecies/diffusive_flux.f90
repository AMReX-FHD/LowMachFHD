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
 
  subroutine diffusive_flux(mla,rho,molarconc,BinvGamma,chi,Gama,flux,dx,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: rho(:) 
    type(multifab) , intent(in   ) :: molarconc(:) 
    type(multifab) , intent(in   ) :: BinvGamma(:)  
    type(multifab) , intent(in   ) :: chi(:)  
    real(kind=dp_t), intent(in   ) :: Gama(:,:)  
    type(multifab) , intent(inout) :: flux(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local variables
    integer :: n,i,dm,ng,nlevs
 
    ! local cell-centered and face-centered multifabs 
    type(multifab)  :: minusrho(mla%nlevel)
    type(multifab)  :: rhoichiGama(mla%nlevel)
    type(multifab)  :: rhoichiGama_face(mla%nlevel,mla%dim)
    type(multifab)  :: BinvGamma_face(mla%nlevel,  mla%dim)
   
    dm    = mla%dim     ! dimensionality
    nlevs = mla%nlevel  ! number of levels 

    ! build local face-centered multifab with nspecies^2 component, zero ghost cells 
    ! and nodal in direction i
    do n=1,nlevs
       call multifab_build(minusrho(n),   mla%la(n),nspecies,rho(n)%ng)
       call multifab_build(rhoichiGama(n),mla%la(n),nspecies,rho(n)%ng)
       do i=1,dm
          call multifab_build_edge(rhoichiGama_face(n,i), mla%la(n), nspecies**2, 0, i)
          call multifab_build_edge(BinvGamma_face(n,i),   mla%la(n), nspecies**2, 0, i)
       enddo
    enddo 

    ! initialize new multifabs to zero 
    do n=1,nlevs
       call setval(minusrho(n),   0.d0,all=.true.)
       call setval(rhoichiGama(n),0.d0,all=.true.)
       do i=1,dm
          call setval(rhoichiGama_face(n,i),0.d0,all=.true.)
          call setval(BinvGamma_face(n,i),  0.d0,all=.true.)
       enddo
    enddo 

    ! calculate face-centrered grad(molarconc) 
    call compute_grad(mla, molarconc, flux, dx, 1, mol_frac_bc_comp, 1, nspecies, & 
                      the_bc_level, .false.)
   
    ! compute face-centered B^(-1)*Gamma from cell-centered values 
    call average_cc_to_face(nlevs, BinvGamma, BinvGamma_face, 1, diff_coeff_bc_comp, & 
                            nspecies**2, the_bc_level, .false.) 
  
    ! check whether nonideal/ideal mixture 
    if(is_ideal_mixture) then
      do n=1,nlevs
         call saxpy(rhoichiGama(n),1.0d0,chi(n))
      enddo
    else
      do n=1,nlevs
       !call matmat_mul(mla,gama(n),chi(n))  ! this will work once gama is multifab
       !call saxpy(rhoichiGama(n),1.0d0,gama(n))
      enddo
    endif 

    ! multiply with -rhoi
    do n=1,nlevs
       call saxpy(minusrho(n),-1.0d0,rho(n))
       call multifab_mult_mult(rhoichiGama(n),minusrho(n),rho(1)%ng)
    enddo

    ! compute face-centered rhoichiGama from cell-centered values 
    !call average_cc_to_face(nlevs, rhoichiGama, rhoichiGama_face, 1, diff_coeff_bc_comp, &
    !                        nspecies**2, the_bc_level, .false.) 
    
    ! compute flux as B^(-1)*Gama X grad(molarconc). 
    do n=1,nlevs
       do i=1,dm
          call matvec_mul(mla, flux(n,i), BinvGamma_face(n,i))
          !!!!!call matvec_mul(mla, flux(n,i), rhoichiGama_face(n,i))
       enddo
    enddo    
     
    ! destroy B^(-1)*Gama multifab to prevent leakage in memory
    do n=1,nlevs
       call multifab_destroy(rhoichiGama(n))
       call multifab_destroy(minusrho(n))
       do i=1,dm
          call multifab_destroy(BinvGamma_face(n,i))
          call multifab_destroy(rhoichiGama_face(n,i))
       enddo
    enddo

  end subroutine diffusive_flux

end module diffusive_flux_module
