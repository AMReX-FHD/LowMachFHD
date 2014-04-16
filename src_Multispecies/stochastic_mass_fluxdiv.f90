module stochastic_mass_fluxdiv_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use convert_mass_variables_module
  use multifab_fill_random_module
  use div_and_grad_module
  use convert_stag_module
  use matvec_mul_module
  use correction_flux_module
  use probin_common_module
  use probin_multispecies_module

  implicit none

  private

  public :: stochastic_mass_fluxdiv, generate_random_increments, destroy_random_increments
  
contains
  
  subroutine stochastic_mass_fluxdiv(mla,rho,rho_tot,molarconc,molmass,molmtot,chi,&
                                     Gama,stoch_W_fc,stoch_fluxdiv,dx,dt,weights,the_bc_level)

    type(ml_layout), intent(in   )   :: mla
    type(multifab) , intent(in   )   :: rho(:)
    type(multifab) , intent(in   )   :: rho_tot(:)
    type(multifab) , intent(in   )   :: molarconc(:)
    real(kind=dp_t), intent(in   )   :: molmass(nspecies) 
    type(multifab) , intent(in   )   :: molmtot(:)
    type(multifab) , intent(in   )   :: chi(:)
    type(multifab) , intent(in   )   :: Gama(:)
    type(multifab) , intent(in   )   :: stoch_W_fc(:,:,:)
    type(multifab) , intent(inout)   :: stoch_fluxdiv(:) 
    real(kind=dp_t), intent(in   )   :: dx(:,:)
    real(kind=dp_t), intent(in   )   :: dt
    real(kind=dp_t), intent(in   )   :: weights(:)         
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
    variance = sqrt(2.d0*k_B*variance_parameter/(product(dx(1,1:dm))*dt))

    ! build multifabs 
    do n=1,nlevs
       call multifab_build(Lonsager(n), mla%la(n), nspecies**2, rho(n)%ng)
       do i=1,dm
          call multifab_build_edge(Lonsager_fc(n,i),   mla%la(n), nspecies**2, 0, i)
          call multifab_build_edge(stoch_flux_fc(n,i), mla%la(n), nspecies,    0, i)
       end do
    end do
 
    ! set stoch_flux_fc to zero
    do n=1,nlevs
       do i = 1,dm
          call setval(stoch_flux_fc(n,i), 0.d0, all=.true.)   
       end do   
    end do   
    
    ! convert stoch_W_fc into stoch_flux_fc
    do n=1,nlevs
       do i = 1,dm
          do rng=1,size(weights)
             call saxpy(stoch_flux_fc(n,i), weights(rng), stoch_W_fc(n,i,rng))
          end do   
       end do   
    end do
    
    ! compute cell-centered cholesky-factored Lonsager^(1/2)
    call compute_Lonsager(mla,rho,rho_tot,molarconc,molmass,molmtot,chi,Gama,Lonsager,the_bc_level)
                  
    ! compute face-centered cholesky factor of cell-centered cholesky factored Lonsager^(1/2)
    call average_cc_to_face(nlevs,Lonsager,Lonsager_fc,1,diff_coeff_bc_comp,nspecies**2,the_bc_level,.false.)

    ! compute variance X cholesky-Lonsager-face X W(0,1) 
    do n=1,nlevs
       do i=1,dm
          call matvec_mul(mla, stoch_flux_fc(n,i), Lonsager_fc(n,i), nspecies)
          call multifab_mult_mult_s(stoch_flux_fc(n,i), variance, 0)
       end do
    end do  
    
    ! sync the fluxes at the boundaries
    do n=1,nlevs
       do i=1,dm
          call multifab_internal_sync(stoch_flux_fc(n,i))
          call multifab_fill_boundary(stoch_flux_fc(n,i))  
       end do
    end do

    !correct fluxes to ensure mass conservation to roundoff
    if (correct_flux .and. (nspecies .gt. 1)) then
       !write(*,*) "Checking conservation of stochastic fluxes"
       call correction_flux(mla, rho, rho_tot, stoch_flux_fc, the_bc_level)
    end if
 
    ! compute divergence of stochastic flux
    call compute_div(mla,stoch_flux_fc,stoch_fluxdiv,dx,1,1,nspecies)

    ! free the multifab allocated memory
    do n=1,nlevs
       call multifab_destroy(Lonsager(n))
       do i=1,dm
          call multifab_destroy(Lonsager_fc(n,i))
          call multifab_destroy(stoch_flux_fc(n,i))
       end do
    end do

  end subroutine stochastic_mass_fluxdiv

  subroutine generate_random_increments(mla,n_rngs,stoch_W_fc)
  
      type(ml_layout), intent(in   )  :: mla
      integer,         intent(in   )  :: n_rngs   ! how many random numbers to store per time step
      type(multifab),  intent(inout)  :: stoch_W_fc(:,:,:)  

      ! Local variables
      integer :: comp,n,dm,nlevs,box,i,rng
    
      nlevs = mla%nlevel
      dm    = mla%dim    
    
      ! generate and store the stochastic flux (random numbers)
      do rng=1, n_rngs
         do i = 1,dm
            call multifab_fill_random(stoch_W_fc(:,i,rng))
         end do   
      end do   
  
  end subroutine generate_random_increments
  
  subroutine destroy_random_increments(mla,n_rngs,stoch_W_fc)
    
      type(ml_layout), intent(in   )  :: mla
      integer,         intent(in   )  :: n_rngs
      type(multifab),  intent(inout)  :: stoch_W_fc(:,:,:)  

      ! Local variables
      integer :: comp,n,dm,nlevs,box,i,rng

      nlevs = mla%nlevel
      dm    = mla%dim    
  
      ! destroy multifab for stochastic flux
      do n=1, nlevs 
         do rng=1, n_rngs 
            do i = 1,dm
               call multifab_destroy(stoch_W_fc(n,i,rng))
            end do
         end do
      end do

  end subroutine destroy_random_increments
   
end module stochastic_mass_fluxdiv_module
