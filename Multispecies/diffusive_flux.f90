module diffusive_flux_module

  use multifab_module
  use define_bc_module
  use bc_module
  use div_and_grad_module
  use probin_multispecies_module
  use ml_layout_module
  use convert_stag_module
  use F95_LAPACK
  
  implicit none

  private

  public :: diffusive_flux

contains
 
  subroutine diffusive_flux(mla,molarconc,BinvGamma,flux,dx,the_bc_level,mol_frac_bc_comp) 

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: molarconc(:) 
    type(multifab) , intent(in   ) :: BinvGamma(:)  
    type(multifab) , intent(inout) :: flux(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    integer,         intent(in   ) :: mol_frac_bc_comp

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
          call multifab_build_edge(BinvGamma_face(n,i),mla%la(n),nspecies**2,0,i)
       end do
    end do 
 
    ! calculate face-centrered grad(molarconc) 
    call compute_grad(mla,molarconc,flux,dx,1,mol_frac_bc_comp,1,nspecies,the_bc_level)
   
    ! compute face-centered B^(-1)*Gamma from cell-centered values 
    do i=1,nspecies**2
       call average_cc_to_face(nlevs,BinvGamma,BinvGamma_face,i,mol_frac_bc_comp,1,the_bc_level) 
    end do
    
    ! compute flux as B^(-1)*Gama X grad(molarconc). 
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_matrixvec_c(flux(n,i),1,BinvGamma_face(n,i),1,nspecies,0)          
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

  subroutine multifab_mult_matrixvec_c(a, targ, b, src, nc, ng)
    
    integer, intent(in)           :: targ, src
    integer, intent(in)           :: nc
    integer, intent(in), optional :: ng
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)    :: b
    real(dp_t), pointer           :: ap(:,:,:,:)
    real(dp_t), pointer           :: bp(:,:,:,:)  ! remember, it's nspecies^2.
    integer                       :: i,lng
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= ng,"not enough ghost cells in multifab_mult_matrixvec_c")
    do i = 1, nlocal(a%la)
       if ( lng > 0 ) then
          ap => dataptr(a, i, grow(get_ibox(a, i),lng), targ, nc)
          bp => dataptr(b, i, grow(get_ibox(b, i),lng), src, nc**2)
       else
          ap => dataptr(a, i, get_ibox(a, i), targ, nc)
          bp => dataptr(b, i, get_ibox(b, i), src, nc**2)
       end if
       call multifab_mult_matrixvec_c_doit(ap, bp)
    end do
  end subroutine multifab_mult_matrixvec_c

   subroutine multifab_mult_matrixvec_c_doit(ap, bp)

    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)

    integer :: i, j, k, n

    !$OMP PARALLEL PRIVATE(i,j,k,n)
    do n=lbound(ap,dim=4), ubound(ap,dim=4)   ! this is 1:nspecies
       !$OMP DO 
       do k = lbound(ap,dim=3), ubound(ap,dim=3)
          do j = lbound(ap,dim=2), ubound(ap,dim=2)
             do i = lbound(ap,dim=1), ubound(ap,dim=1)
                !print*, ap(i,j,k,:)
                call matvec_mul(ap(i,j,k,:), bp(i,j,k,:), n)
             end do
          end do
       end do
       !$OMP END DO NOWAIT
    end do
    !$OMP END PARALLEL

    ! Use contained (internal) subroutine to do the rank conversion 
    contains 
     subroutine matvec_mul(ap_ij, bp_ij, n)
        real(kind=dp_t), dimension(nspecies),          intent(inout) :: ap_ij
        real(kind=dp_t), dimension(nspecies,nspecies), intent(in)    :: bp_ij  
        integer,                                       intent(in)    :: n      
       
        ! local variables
        real(kind=dp_t) :: mvprod
        integer         :: m      
 
        mvprod=0.d0
        do m=1, nspecies
           !print*, bp_ij(m,n), ap_ij(m)
           mvprod = mvprod + bp_ij(n,m)*ap_ij(m)
        enddo
        ap_ij(n) = mvprod
        !print*, ap_ij(n)
     
     end subroutine 

  end subroutine multifab_mult_matrixvec_c_doit

end module diffusive_flux_module
