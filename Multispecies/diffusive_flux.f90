module diffusive_flux_module

  use multifab_module
  use define_bc_module
  use div_and_grad_module
  use probin_multispecies_module
  use ml_layout_module
  use convert_stag_module
  use matvec_mul_module
  use matmat_mul_module
  use correction_flux_module
  use F95_LAPACK
  
  implicit none

  private

  public :: diffusive_flux

contains
 
  subroutine diffusive_flux(mla,rho,rho_tot,molarconc,rhoWchi,Gama,Temp,&
                            zeta_by_Temp,flux,dx,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: rho(:) 
    type(multifab) , intent(in   ) :: rho_tot(:) 
    type(multifab) , intent(in   ) :: molarconc(:) 
    type(multifab) , intent(in   ) :: rhoWchi(:)  
    type(multifab) , intent(in   ) :: Gama(:)  
    type(multifab) , intent(in   ) :: Temp(:)  
    type(multifab) , intent(in   ) :: zeta_by_Temp(:)  
    type(multifab) , intent(inout) :: flux(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local variables
    integer :: n,i,dm,ng,nlevs
 
    ! local face-centered multifabs 
    type(multifab)  :: rhoWchi_face(mla%nlevel,mla%dim)
    type(multifab)  :: Gama_face(mla%nlevel,mla%dim)
    type(multifab)  :: zeta_by_Temp_face(mla%nlevel,mla%dim)
    type(multifab)  :: flux_Temp(mla%nlevel,mla%dim)
  
 
    dm    = mla%dim     ! dimensionality
    nlevs = mla%nlevel  ! number of levels 

    ! build local face-centered multifab with nspecies^2 component, zero ghost cells 
    ! and nodal in direction i
    do n=1,nlevs
       do i=1,dm
          call multifab_build_edge(rhoWchi_face(n,i),     mla%la(n),nspecies**2,0,i)
          call multifab_build_edge(Gama_face(n,i),        mla%la(n),nspecies**2,0,i)
          call multifab_build_edge(zeta_by_Temp_face(n,i),mla%la(n),nspecies,   0,i)
          call multifab_build_edge(flux_Temp(n,i),        mla%la(n),1,          0,i)
       end do
    end do 

    ! compute face-centered rhoWchi from cell-centered values 
    call average_cc_to_face(nlevs, rhoWchi, rhoWchi_face, 1, diff_coeff_bc_comp, &
                            nspecies**2, the_bc_level, .false.) 

    !==================================!
    ! compute flux-piece from molarconc
    !==================================! 

    ! calculate face-centrered grad(molarconc) 
    call compute_grad(mla, molarconc, flux, dx, 1, mol_frac_bc_comp, 1, nspecies, & 
                      the_bc_level)

    ! compute face-centered Gama from cell-centered values 
    call average_cc_to_face(nlevs, Gama, Gama_face, 1, diff_coeff_bc_comp, &
                            nspecies**2, the_bc_level, .false.)

    ! compute rhoWchi X Gama (on faces) 
    do n=1,nlevs
       do i=1,dm
          call matmat_mul(mla, Gama_face(n,i), rhoWchi_face(n,i), nspecies)
       end do
    end do    
    
    ! compute flux from molarconc as rhoWchi X Gama X grad(molarconc) 
    do n=1,nlevs
       do i=1,dm
          call matvec_mul(mla, flux(n,i), Gama_face(n,i), nspecies)
       end do
    end do    

    !====================================!
    ! compute flux-piece from Temperature 
    !====================================! 
   
    ! calculate face-centrered grad(Temp) 
    call compute_grad(mla, Temp, flux_Temp, dx, 1, mol_frac_bc_comp, 1, 1, the_bc_level)
    
    ! compute face-centered zeta_by_Temp from cell-centered values 
    call average_cc_to_face(nlevs, zeta_by_Temp, zeta_by_Temp_face, 1, diff_coeff_bc_comp, &
                            nspecies, the_bc_level, .false.) 
 
    ! compute rhoWchi X zeta_by_Temp (on faces) 
    do n=1,nlevs
       do i=1,dm
          call matvec_mul(mla, zeta_by_Temp_face(n,i), rhoWchi_face(n,i), nspecies)
       end do
    end do    
   
    ! compute rhoWchi X zeta_by_Temp X grad(Temp) 
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult(zeta_by_Temp_face(n,i), flux_Temp(n,i), 0)
       end do
    end do  
    
    !===============================!
    ! assemble different flux-pieces 
    !===============================! 
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus(flux(n,i), zeta_by_Temp_face(n,i), 0)
       end do
    end do  

    !correct fluxes to ensure mass conservation to roundoff
    if (correct_flux .and. (nspecies .gt. 1)) then
       !write(*,*) "Checking conservation of deterministic fluxes"
       call correction_flux(mla, rho, rho_tot, flux, the_bc_level)
    end if
    
    ! destroy B^(-1)*Gama multifab to prevent leakage in memory
    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(rhoWchi_face(n,i))
          call multifab_destroy(Gama_face(n,i))
          call multifab_destroy(zeta_by_Temp_face(n,i))
          call multifab_destroy(flux_Temp(n,i))
       end do
    end do

  end subroutine diffusive_flux

end module diffusive_flux_module
