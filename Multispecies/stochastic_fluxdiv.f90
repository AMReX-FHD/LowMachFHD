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

  public :: stochastic_fluxdiv
  ! Donev: In this code stochastic_w1 or w2 should never appear, only weights
  ! This code should work the same for any number of nrngs, not just one or two
  
contains

  subroutine stochastic_fluxdiv(mla,stoch_W_fc,stoch_fluxdiv,rho,rho_tot,molarconc,molmass,&
                                molmtot,chi,Gama,Lonsager,dx,weights,the_bc_level)

    type(ml_layout), intent(in   )   :: mla
    type(multifab) , intent(inout)   :: stoch_W_fc(:,:,:)
    type(multifab) , intent(inout)   :: stoch_fluxdiv(:) 
    type(multifab) , intent(inout)   :: rho(:)
    type(multifab) , intent(inout)   :: rho_tot(:)
    type(multifab) , intent(inout)   :: molarconc(:)
    real(kind=dp_t), intent(in   )   :: molmass(nspecies) 
    type(multifab) , intent(inout)   :: molmtot(:)
    type(multifab) , intent(inout)   :: chi(:)
    type(multifab) , intent(inout)   :: Gama(:)
    type(multifab) , intent(inout)   :: Lonsager(:)
    real(kind=dp_t), intent(in   )   :: dx(:,:)
    ! Donev: We always have to pass one weight per component of stoch_W_fc
    ! The old code you are copying sometimes set n_rng=0 which meant to use the same random number in all stages
    ! Do not do this -- simply assume n_rngs>0 always
    real(kind=dp_t), intent(in)      :: weights(:)        ! reuse previously-generated rngs
    type(bc_level) , intent(in   )   :: the_bc_level(:)

    ! Local variables
    type(multifab) :: Lonsager_fc(mla%nlevel,mla%dim)    ! cholesky factor of Lonsager on face
    type(multifab) :: stoch_flux_fc(mla%nlevel,mla%dim)  ! face-centered stochastic flux
    integer        :: i,n,dm,nlevs,box,rng
    logical        :: reuse

    nlevs = mla%nlevel
    dm    = mla%dim
  
    ! build multifabs 
    do n=1,nlevs
       do i=1,dm
          call multifab_build_edge(Lonsager_fc(n,i),    mla%la(n), nspecies, 0, i)
          call multifab_build_edge(stoch_flux_fc(n,i), mla%la(n), nspecies, 0, i)
       enddo
    enddo
 
    ! Donev:
    ! In the code that you copied here, stoch_W_fc(n,i,0), i.e., rng=0
    ! corresponds to your variable stoch_flux_fc
    ! So either remove the variable stoch_flux_fc and use stoch_W_fc(n,i,0) to store those numbers
    ! or make a local multifab stoch_flux_fc and allocate stoch_W_fc to be:
    ! allocate(stoch_W_fc(mla%nlevel,mla%dim,nspecies,1:n_rngs)) ! Not 0:n_rngs
    ! I write the code below assuming we are using stoch_flux_fc

    ! convert stoch_W_fc into stoch_flux_fc
    do i = 1,dm
       do rng=1, size(weights)
          call saxpy(stoch_flux_fc(n,i), weights(rng), stoch_W_fc(n,i,rng))
       enddo   
    enddo
    
    ! compute cell-centered cholesky-factored Lonsager
    call compute_Lonsager(mla,rho,rho_tot,molarconc,molmass,molmtot,chi,Gama,Lonsager,the_bc_level)
                  
    ! compute face-centered cholesky factor of cell-centered cholesky factored Lonsager
    call average_cc_to_face(nlevs,Lonsager,Lonsager_fc,1,diff_coeff_bc_comp,nspecies**2,the_bc_level,.false.)

    ! compute sqrt(2*k_B) X cholesky-Lonsager-face X W(0,1) 
    do n=1,nlevs
       do i=1,dm
          call matvec_mul(mla, stoch_flux_fc(n,i), Lonsager_fc(n,i))
          ! Donev: variance is an undefined variable here
          ! Amit: variance is populated in advance and passed as global varialbe
          ! via probin_multispecies.
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
       do i=1,dm
          call multifab_destroy(Lonsager_fc(n,i))
          call multifab_destroy(stoch_flux_fc(n,i))
       enddo
    enddo

  end subroutine stochastic_fluxdiv
   
end module stochastic_fluxdiv_module
