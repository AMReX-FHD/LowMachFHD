module stochastic_fluxdiv_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use ml_layout_module
  !use ml_restriction_module
  use define_bc_module
  use bc_module
  use BoxLibRNGs
  use multifab_fill_random_module
  use multifab_filter_module
  use probin_multispecies_module

  implicit none

  private

  public :: stochastic_fluxdiv,create_random_increments,destroy_random_increments, &
            add_stochastic_fluxdiv,reuse_stochastic_fluxdiv

contains

  subroutine stochastic_fluxdiv(mla,stoch_flux_fc,stoch_W_fc,stoch_fluxdiv,rho,rho_tot,molarconc,molmass,&
                                molmtot,chi,Gama,Lonsager,dx,weights,the_bc_level)

    type(ml_layout), intent(in   )   :: mla
    type(multifab) , intent(inout)   :: stoch_flux_fc(:,:)
    type(multifab) , intent(inout)   :: stoch_W_fc(:,:,:,:)
    type(multifab) , intent(inout)   :: stoch_fluxdiv(:) 
    type(multifab) , intent(inout)   :: rho(:)
    type(multifab) , intent(inout)   :: rho_tot(:)
    type(multifab) , intent(inout)   :: molarconc(:)
    real(kind=dp_t), intent(in   )   :: molmass(nspecies) 
    type(multifab) , intent(inout)   :: molmtot(:)
    type(multifab) , intent(inout)   :: chi(:)
    type(multifab) , intent(inout)   :: Gama(:)
    type(multifab) , intent(inout)   :: Lonsager(:)
    real(dp_t)     , intent(in   )   :: dx(:,:)
    real(dp_t), intent(in), optional :: weights(:) ! if present, reuse previously-generated rngs
    type(bc_level) , intent(in   )   :: the_bc_level(:)

    ! Local variables
    type(multifab) :: Lonsager_f(mla%nlevel) ! cholesky factor of Lonsager(face) 
    integer        :: i,n,dm,nlevs,box,rng
    logical        :: reuse

    nlevs = mla%nlevel
    dm    = mla%dim
    reuse=.false.

    if(present(weights)) then ! Make a weighted sum of previously-generated random numbers
       if(size(weights)>0) then
         !write(*,*) "REUSING Weiner increments:", weights
          reuse=.true.
          do i = 1,dm
             !if(ntracers>0) then
             !   call multifab_weighted_sum(stoch_W_fc(:,i,:,:), weights)
             !endif   
          enddo
       endif          
    endif
   
    ! convert stoch_W_fc into stoch_flux_fc
    do n=1,nlevs
       do i=1,dm 
          call saxpy(stoch_flux_fc(n,i), stochastic_w1, stoch_W_fc(n,i,:,0))
          call saxpy(stoch_flux_fc(n,i), stochastic_w2, stoch_W_fc(n,i,:,1))
       enddo
    enddo 
 
    ! compute cell-centered cholesky-factored Lonsager
    call compute_Lonsager(mla,rho,rho_tot,molarconc,molmass,molmtot,chi,Gama,Lonsager,the_bc_level)
                  
    ! compute face-centered cholesky factor of cell-centered cholesky factored Lonsager
    call average_cc_to_face(nlevs,Lonsager,Lonsager_f,1,diff_coeff_bc_comp,nspecies**2,the_bc_level,.false.)

    ! compute sqrt(2*k_B) X cholesky-Lonsager-face X W(0,1) 
    do n=1,nlevs
       do i=1,dm
          call matvec_mul(mla, stoch_flux_fc(n,i), Lonsager_f(n,i))
          call multifab_mult_mult_s(stoch_flux_fc(n,i), variance, 0)
       enddo
    enddo  
     
    ! compute divergence of stochastic flux
    call compute_div(mla,stoch_flux_fc,stoch_fluxdiv,dx,1,1,nspecies)

    ! multiply fluxdiv (having zero ghost cells) with -1 to get -div(-flux).
    do n=1,nlevs
       call multifab_mult_mult_s(stoch_fluxdiv(n),-1.0d0,stoch_fluxdiv(1)%ng)
    enddo
 
    ! free the multifab allocated memory
    do n=1,nlevs
       call multifab_destroy(Lonsager_f(n))
    enddo

  end subroutine stochastic_fluxdiv

  subroutine create_random_increments(mla,n_rngs,stoch_W_fc)
  
    type(ml_layout), intent(in   )  :: mla
    integer,         intent(in   )  :: n_rngs  ! how many random numbers to store per time step
                                               ! could be zero if one does not need to store any rngs
    type(multifab),  intent(inout)  :: stoch_W_fc(:,:,:,:)  

    ! Local variables
    integer :: comp,n,dm,nlevs,box,i,rng
    
    nlevs = mla%nlevel
    dm    = mla%dim    
    
    ! generate and store the stochastic diffusive flux (random numbers)
       do rng=1, n_rngs
          do comp=1, nspecies
             do i = 1,dm
                call multifab_fill_random(stoch_W_fc(:,i,comp,rng))
             enddo   
          enddo   
       enddo   
  
  end subroutine create_random_increments
  
  subroutine destroy_random_increments(mla,n_rngs,stoch_W_fc)
    
    type(ml_layout), intent(in   )  :: mla
    integer,         intent(in   )  :: n_rngs
    type(multifab),  intent(inout)  :: stoch_W_fc(:,:,:,:)  

    ! Local variables
    integer :: comp,n,dm,nlevs,box,i,rng

    nlevs = mla%nlevel
    dm    = mla%dim    
  
    ! deallocate stochastic diffusive flux 
    do n = 1, nlevs
       do rng=0, n_rngs 
          do comp=1, nspecies
             do i = 1,dm
                call multifab_destroy(stoch_W_fc(n,i,comp,rng))
             enddo
          enddo
       enddo
    enddo

  end subroutine destroy_random_increments
 
  subroutine add_stochastic_fluxdiv(mla,fluxdiv,stoch_fluxdiv,n_rngs)
 
    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(inout)  :: fluxdiv(:)
    type(multifab) , intent(in   )  :: stoch_fluxdiv(:) 
    integer        , intent(in   )  :: n_rngs 
 
    ! Local variables
    integer :: comp,n,dm,nlevs,box,rng
    
    if(n_rngs==0) then ! There is no weighting in stochastic_fluxdiv, so we do it here
       do n=1,nlevs
          call saxpy(fluxdiv(n), stochastic_w1, stoch_fluxdiv(n))
       enddo
    else ! The weighting was already done in stochastic_fluxdiv, so just add this to the sum
       do n=1,nlevs
          call multifab_plus_plus(fluxdiv(n), stoch_fluxdiv(n))
       enddo
    endif   

  end subroutine add_stochastic_fluxdiv

  subroutine reuse_stochastic_fluxdiv(mla,fluxdiv,stoch_fluxdiv,n_rngs)
    
    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(inout)  :: fluxdiv(:)
    type(multifab) , intent(in   )  :: stoch_fluxdiv(:) 
    integer        , intent(in   )  :: n_rngs 

    ! Local variables
    integer :: comp,n,dm,nlevs,box,rng

    do n=1,nlevs
       call multifab_plus_plus(fluxdiv(n), stoch_fluxdiv(n))
    enddo

  end subroutine reuse_stochastic_fluxdiv

  subroutine multifab_weighted_sum(mfab, weights)
  
    type(multifab), intent(inout) :: mfab(:,0:)
    real(dp_t), intent(in)        :: weights(:)

    integer :: i,n,dm,rng
  
    !--------------------------------------
    ! ghost cells need not be set here since syncs will set them later
    do n=1,size(mfab,1)
       call setval(mfab(n,0), 0.d0)
       ! we want to do mfab(n,0)=sum(weights(i)*mfab(n,i),i=1..n_rngs)
       do rng=1,size(weights)
         call saxpy(mfab(n,0), weights(rng), mfab(n,rng))
       enddo
    enddo   

  end subroutine multifab_weighted_sum

end module stochastic_fluxdiv_module
