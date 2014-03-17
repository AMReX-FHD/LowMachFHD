module advance_module

  use multifab_module
  use define_bc_module
  use multifab_physbc_module
  use ml_layout_module
  use diffusive_fluxdiv_module
  use stochastic_fluxdiv_module
  use probin_multispecies_module

  implicit none

  private

  public :: advance

contains

  subroutine advance(mla,rho,molmass,dx,dt,time,prob_lo,prob_hi,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rho(:)
    real(kind=dp_t), intent(in   ) :: molmass(nspecies) 
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(in   ) :: dt,time
    real(kind=dp_t), intent(in   ) :: prob_lo(rho(1)%dim),prob_hi(rho(1)%dim) 
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local variables and array of multifabs
    type(multifab)                 :: rhonew(mla%nlevel),fluxdiv(mla%nlevel),fluxdivnew(mla%nlevel)
    type(multifab)                 :: rho_tot(mla%nlevel)                ! total density
    type(multifab)                 :: molarconc(mla%nlevel)              ! molar concentration
    type(multifab)                 :: molmtot(mla%nlevel)                ! total molar mass
    type(multifab)                 :: chi(mla%nlevel)                    ! Chi-matrix
    type(multifab)                 :: Lonsager(mla%nlevel)               ! Onsager matrix for fluctuations
    type(multifab)                 :: D_MS(mla%nlevel)                   ! D_MS-matrix
    type(multifab)                 :: Gama(mla%nlevel)                   ! Gama-matrix
    type(multifab)                 :: rhoWchiGama(mla%nlevel)            ! rho*W*chi*Gama
    type(multifab)                 :: stoch_fluxdiv(mla%nlevel)          ! stochastic fluxdiv
    type(multifab)                 :: stoch_flux_fc(mla%nlevel,mla%dim)  ! face-centered stochastic flux
    type(multifab),  allocatable   :: stoch_W_fc(:,:,:,:)                ! WA and WB (nlevs,dim,nspecies,n_rngs) 
    real(kind=dp_t), allocatable   :: weights(:)                         ! weights for RNGs       
    real(kind=dp_t)                :: variance,stage_time
    integer                        :: n,i,dm,n_rngs,nlevs
 
    nlevs = mla%nlevel  ! number of levels 
    dm    = mla%dim     ! number of dimensions

    ! build cell-centered multifabs for nspecies and ghost cells contained in rho.
    ! fluxdiv needs no ghost cells as all computations will be in box interior. 
    ! rho_tot is cell-cented with ghost cells that many rho owns.
    do n=1,nlevs
       call multifab_build(rhonew(n),        mla%la(n), nspecies,    rho(n)%ng)
       call multifab_build(rho_tot(n),       mla%la(n), 1,           rho(n)%ng) 
       call multifab_build(fluxdiv(n),       mla%la(n), nspecies,    0) 
       call multifab_build(fluxdivnew(n),    mla%la(n), nspecies,    0)
       call multifab_build(molmtot(n),       mla%la(n), 1,           rho(n)%ng)  
       call multifab_build(molarconc(n),     mla%la(n), nspecies,    rho(n)%ng) 
       call multifab_build(chi(n),           mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(Lonsager(n),      mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(D_MS(n),          mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(Gama(n),          mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(rhoWchiGama(n),   mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(stoch_fluxdiv(n), mla%la(n), nspecies,    0) 
       do i=1,dm
          call multifab_build_edge(stoch_flux_fc(n,i), mla%la(n), nspecies, 0, i)
       enddo
    enddo

    !========================================================
    ! Initialize random number generator for stochastic flux
    !========================================================
    if(use_stoch) then
 
       if (timeinteg_type .eq. 1) then      ! Euler-Maruyama
          n_rngs=0
          stochastic_w1 = 1.d0
          stochastic_w2 = 0.d0
       elseif (timeinteg_type .eq. 2) then  ! Trapezoidal
          n_rngs=1
          stochastic_w1 = 1.d0        !Amit: has to be corrected 
          stochastic_w2 = 0.d0
       else                                 ! Midpoint & RK3
          n_rngs=2  
          stochastic_w1 = 1.d0        !Amit: has to be corrected 
          stochastic_w2 = 0.d0
       endif
 
       ! The zeroth component is used to store the actual stochastic flux
       ! The rest are used to store random numbers that may be reused later
       allocate(stoch_W_fc(mla%nlevel,mla%dim,nspecies,0:n_rngs))
       allocate(weights(n_rngs))
    
       do n=1,nlevs 
          variance = sqrt(2.d0*kT/(product(dx(n,1:dm))*dt))
       enddo
 
       ! initialize stochastic flux on every face W(0,1) 
       call create_random_increments(mla,n_rngs,stoch_W_fc)

    endif
 
   !==================================================================================
    select case(timeinteg_type)
   !==================================================================================
 
      case(1)
      !===================================
      ! Euler time update 
      ! rho(t+dt)  =rho(t) + dt*fluxdiv(t) 
      !===================================
      
      ! compute the total div of flux from rho
      stage_time = time   
      call compute_fluxdiv(mla,rho,stoch_flux_fc,stoch_W_fc,fluxdiv,rho_tot,molarconc,molmtot,molmass,&
           chi,Lonsager,Gama,D_MS,dx,rhoWchiGama,stoch_fluxdiv,stage_time,prob_lo,prob_hi,&
           weights,n_rngs,the_bc_level)

      ! compute rho(t+dt) (only interior) 
      do n=1,nlevs
         call saxpy(rho(n),dt,fluxdiv(n))
      enddo 
    
      case(2)
      !============================================================================
      ! Heun's method: Predictor-Corrector explicit method 
      ! rhonew(t+dt) = rho(t)+dt*fluxdiv(t,rho) 
      !    rho(t+dt) = rho(t)+(dt/2)*[fluxdiv(t,rho) + fluxdivnew(t+1,rhonew(t+1))] 
      !============================================================================
      
      ! store old rho in rhonew (rhonew previously set to zero) 
      do n=1,nlevs
         call saxpy(rhonew(n),1.0d0,rho(n))
      enddo 
      
      !=====================
      ! Euler Predictor step
      !===================== 
      
      ! compute the total div of flux from rho
      stage_time = time 
      call compute_fluxdiv(mla,rho,stoch_flux_fc,stoch_W_fc,fluxdiv,rho_tot,molarconc,molmtot,molmass,&
           chi,Lonsager,Gama,D_MS,dx,rhoWchiGama,stoch_fluxdiv,stage_time,prob_lo,prob_hi,&
           weights,n_rngs,the_bc_level)
      
      ! compute rhonew(t+dt) (only interior) 
      do n=1,nlevs
         call saxpy(rhonew(n),dt,fluxdiv(n))
      enddo 
   
      ! update values of the ghost cells of rhonew
      do n=1,nlevs
         call multifab_fill_boundary(rhonew(n))
        
         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc(rhonew(n),1,rho_part_bc_comp,nspecies,the_bc_level(n),dx(n,:))
      enddo

      !=========================== 
      ! Trapezoidal Corrector step
      !===========================
      
      ! compute the total div of flux from rho
      stage_time = time + dt  
      call compute_fluxdiv(mla,rhonew,stoch_flux_fc,stoch_W_fc,fluxdiv,rho_tot,molarconc,molmtot,molmass,&
           chi,Lonsager,Gama,D_MS,dx,rhoWchiGama,stoch_fluxdiv,stage_time,prob_lo,prob_hi,&
           weights,n_rngs,the_bc_level)

      ! compute rho(t+dt) (only interior)
      do n=1,nlevs
         call saxpy(rho(n),dt*0.5d0,fluxdiv(n))
         call saxpy(rho(n),dt*0.5d0,fluxdivnew(n))
      enddo
 
      case(3)
      !============================================================
      ! Midpoint method (2-stage) 
      ! rhonew(t+dt/2)=rho(t)+(dt/2)*fluxdiv(t,y)                
      !    rho(t+dt)  =rho(t)+ dt*fluxdivnew[t+dt/2,rhonew(t+dt/2)]  
      !============================================================

      ! store old rho in rhonew (rhonew is set zero at the start) 
      do n=1,nlevs
         call saxpy(rhonew(n),1.0d0,rho(n))
      enddo 

      ! compute the total div of flux from rho
      stage_time = time 
      call compute_fluxdiv(mla,rho,stoch_flux_fc,stoch_W_fc,fluxdiv,rho_tot,molarconc,molmtot,molmass,&
           chi,Lonsager,Gama,D_MS,dx,rhoWchiGama,stoch_fluxdiv,stage_time,prob_lo,prob_hi,&
           weights,n_rngs,the_bc_level)
 
      ! compute rhonew(t+dt/2) (only interior) 
      do n=1,nlevs
         call saxpy(rhonew(n),0.5d0*dt,fluxdiv(n))
      enddo 

      ! update values of the ghost cells of rhonew
      do n=1,nlevs
         call multifab_fill_boundary(rhonew(n))
        
         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc(rhonew(n),1,rho_part_bc_comp,nspecies,the_bc_level(n),dx(n,:))
      enddo

      ! compute the total div of flux from rho
      stage_time = time + dt/2.0d0
      call compute_fluxdiv(mla,rhonew,stoch_flux_fc,stoch_W_fc,fluxdiv,rho_tot,molarconc,molmtot,molmass,&
           chi,Lonsager,Gama,D_MS,dx,rhoWchiGama,stoch_fluxdiv,stage_time,prob_lo,prob_hi,&
           weights,n_rngs,the_bc_level)
 
      ! compute rho(t+dt) (only interior) 
      do n=1,nlevs
         call saxpy(rho(n),dt,fluxdivnew(n))
      enddo 
      
      case(4)
      !========================================================================================
      ! 3rd order Runge-Kutta method 
      ! rhonew(t+dt)  =     rho(t)+                       dt*fluxdiv(t)                      
      ! rhonew(t+dt/2)= 3/4*rho(t)+ 1/4*[rhonew(t+dt)   + dt*fluxdivnew(t+dt,rhonew(t+dt))] 
      !    rho(t+dt)  = 1/3*rho(t)+ 2/3*[rhonew(t+dt/2) + dt*fluxdivnew(t+dt/2,rhonew(t+dt/2))] 
      !========================================================================================

      ! store old rho in rhonew (rhonew is set zero at the start) 
      do n=1,nlevs
         call saxpy(rhonew(n),1.0d0,rho(n))
      enddo 
      
      !===========
      ! 1st stage
      !===========

      ! compute the total div of flux from rho
      stage_time = time 
      call compute_fluxdiv(mla,rho,stoch_flux_fc,stoch_W_fc,fluxdiv,rho_tot,molarconc,molmtot,molmass,&
           chi,Lonsager,Gama,D_MS,dx,rhoWchiGama,stoch_fluxdiv,stage_time,prob_lo,prob_hi,&
           weights,n_rngs,the_bc_level)
 
      ! compute rhonew(t+dt) (only interior) 
      do n=1,nlevs
         call saxpy(rhonew(n),dt,fluxdiv(n))
      enddo 
   
      ! update values of the ghost cells of rhonew
      do n=1,nlevs
         call multifab_fill_boundary(rhonew(n))
        
         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc(rhonew(n),1,rho_part_bc_comp,nspecies,the_bc_level(n),dx(n,:))
      enddo

      ! compute the total div of flux from rho
      stage_time = time + dt
      call compute_fluxdiv(mla,rhonew,stoch_flux_fc,stoch_W_fc,fluxdiv,rho_tot,molarconc,molmtot,molmass,&
           chi,Lonsager,Gama,D_MS,dx,rhoWchiGama,stoch_fluxdiv,stage_time,prob_lo,prob_hi,&
           weights,n_rngs,the_bc_level)

      !===========
      ! 2nd stage
      !===========

      ! compute rhonew(t+dt/2) (only interior) 
      do n=1,nlevs
         call multifab_mult_mult_s(rhonew(n),0.25d0)
         call saxpy(rhonew(n), 0.75d0,   rho(n))
         call saxpy(rhonew(n), 0.25d0*dt,fluxdivnew(n))
      enddo 

      ! update values of the ghost cells of rho_star
      do n=1,nlevs
         call multifab_fill_boundary(rhonew(n))
        
         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc(rhonew(n),1,rho_part_bc_comp,nspecies,the_bc_level(n),dx(n,:))
      enddo
      
      ! free up the values of fluxdivnew
      do n=1,nlevs
         call setval(fluxdivnew(n), 0.d0, all=.true.)
      enddo 

      ! compute the total div of flux from rho
      stage_time = time + dt/2.0d0
      call compute_fluxdiv(mla,rhonew,stoch_flux_fc,stoch_W_fc,fluxdiv,rho_tot,molarconc,molmtot,molmass,&
           chi,Lonsager,Gama,D_MS,dx,rhoWchiGama,stoch_fluxdiv,stage_time,prob_lo,prob_hi,&
           weights,n_rngs,the_bc_level)

      !===========
      ! 3rd stage
      !===========

      ! compute rho(t+dt) (only interior) 
      do n=1,nlevs
         call multifab_mult_mult_s(rho(n),(1.0d0/3.0d0))
         call saxpy(rho(n), (2.0d0/3.0d0),   rhonew(n))
         call saxpy(rho(n), (2.0d0/3.0d0)*dt,fluxdivnew(n))
      enddo 

    !=============================================================================================
    end select  
    !=============================================================================================

    ! fill the ghost cells of rho
    do n=1, nlevs
       call multifab_fill_boundary(rho(n))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(rho(n),1,rho_part_bc_comp,nspecies,the_bc_level(n),dx(n,:))
    enddo   

    ! free the multifab allocated memory
    do n=1,nlevs
       call multifab_destroy(rhonew(n))
       call multifab_destroy(fluxdiv(n))
       call multifab_destroy(fluxdivnew(n))
       call multifab_destroy(rho_tot(n))
       call multifab_destroy(molarconc(n))
       call multifab_destroy(molmtot(n))
       call multifab_destroy(chi(n))
       call multifab_destroy(Lonsager(n))
       call multifab_destroy(D_MS(n))
       call multifab_destroy(Gama(n))
       call multifab_destroy(rhoWchiGama(n))
       call multifab_destroy(stoch_fluxdiv(n))
       do i=1,dm
          call multifab_destroy(stoch_flux_fc(n,i))
       enddo
    enddo
    
    if(use_stoch) then 
       call destroy_random_increments(mla,n_rngs)
       deallocate(stoch_W_fc)
       deallocate(weights(n_rngs))
    endif

  end subroutine advance

end module advance_module 
