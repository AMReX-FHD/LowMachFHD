module stochastic_fluxdiv_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use ml_layout_module
  !use ml_restriction_module
  use define_bc_module
  use bc_module
  use convert_variables_module 
  use div_and_grad_module
  use convert_stag_module
  use matvec_mul_module

  use probin_multispecies_module

  implicit none

  private

  public :: stochastic_fluxdiv
  
contains
  
  subroutine stochastic_fluxdiv(mla,rho,rho_tot,molarconc,molmass,molmtot,chi,&
                                Gama,stoch_W_fc,stoch_fluxdiv,dx,dt,weights,the_bc_level)

    ! Donev: For those things that are pure input change intent(inout) to intent(in)
    type(ml_layout), intent(in   )   :: mla
    type(multifab) , intent(inout)   :: rho(:)
    type(multifab) , intent(inout)   :: rho_tot(:)
    type(multifab) , intent(inout)   :: molarconc(:)
    real(kind=dp_t), intent(in   )   :: molmass(nspecies) 
    type(multifab) , intent(inout)   :: molmtot(:)
    type(multifab) , intent(inout)   :: chi(:)
    type(multifab) , intent(inout)   :: Gama(:)
    type(multifab) , intent(inout)   :: stoch_W_fc(:,:,:)
    type(multifab) , intent(inout)   :: stoch_fluxdiv(:) 
    real(kind=dp_t), intent(in   )   :: dx(:,:)
    real(kind=dp_t), intent(in   )   :: dt
    real(kind=dp_t), intent(in)      :: weights(:)         
    type(bc_level) , intent(in   )   :: the_bc_level(:)

    ! Local variables
    type(multifab)   :: Lonsager(mla%nlevel)               ! cholesky factored Lonsager 
    type(multifab)   :: Lonsager_fc(mla%nlevel,mla%dim)    ! cholesky factored Lonsager on face
    type(multifab)   :: stoch_flux_fc(mla%nlevel,mla%dim)  ! face-centered stochastic flux
    integer          :: n,nlevs,i,dm,rng
    real(kind=dp_t)  :: variance

    nlevs = mla%nlevel
    dm    = mla%dim
 
    ! populate the variance (only first level) 
    variance = sqrt(2.d0*k_B/(product(dx(1,1:dm))*dt))

    ! build multifabs 
    do n=1,nlevs
       call multifab_build(Lonsager(n), mla%la(n), nspecies**2, rho(n)%ng)
       do i=1,dm
          call multifab_build_edge(Lonsager_fc(n,i),   mla%la(n), nspecies**2, 0, i)
          call multifab_build_edge(stoch_flux_fc(n,i), mla%la(n), nspecies,    0, i)
       enddo
    enddo
 
    ! set stoch_flux_fc to zero
    do n=1,nlevs
       do i = 1,dm
          call setval(stoch_flux_fc(n,i), 0.d0, all=.true.)   
       enddo   
    enddo   
    
    ! convert stoch_W_fc into stoch_flux_fc
    do n=1,nlevs
       do i = 1,dm
          do rng=1,size(weights)
             call saxpy(stoch_flux_fc(n,i), weights(rng), stoch_W_fc(n,i,rng))
          enddo   
       enddo   
    enddo
    
    ! compute cell-centered cholesky-factored Lonsager^(1/2)
    call compute_Lonsager(mla,rho,rho_tot,molarconc,molmass,molmtot,chi,Gama,Lonsager,the_bc_level)
                  
    ! compute face-centered cholesky factor of cell-centered cholesky factored Lonsager^(1/2)
    call average_cc_to_face(nlevs,Lonsager,Lonsager_fc,1,diff_coeff_bc_comp,nspecies**2,the_bc_level,.false.)

    ! compute sqrt(2k_B) X cholesky-Lonsager-face X W(0,1) 
    do n=1,nlevs
       do i=1,dm
          call matvec_mul(mla, stoch_flux_fc(n,i), Lonsager_fc(n,i))
          call multifab_mult_mult_s(stoch_flux_fc(n,i), variance, 0)
       enddo
    enddo  
     
    ! compute divergence of stochastic flux
    call compute_div(mla,stoch_flux_fc,stoch_fluxdiv,dx,1,1,nspecies)

    ! free the multifab allocated memory
    do n=1,nlevs
       call multifab_destroy(Lonsager(n))
       do i=1,dm
          call multifab_destroy(Lonsager_fc(n,i))
          call multifab_destroy(stoch_flux_fc(n,i))
       enddo
    enddo

  end subroutine stochastic_fluxdiv
   
end module stochastic_fluxdiv_module
