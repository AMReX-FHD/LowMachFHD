module advance_module

  use multifab_module
  use define_bc_module
  use multifab_physbc_module
  use ml_layout_module
  use diffusive_mass_fluxdiv_module
  use stochastic_mass_fluxdiv_module
  use probin_multispecies_module

  implicit none

  private

  public :: advance

contains

  subroutine advance(mla,rho,molmass,Temp,dx,dt,time,prob_lo,prob_hi,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rho(:)
    real(kind=dp_t), intent(in   ) :: molmass(nspecies) 
    type(multifab) , intent(inout) :: Temp(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(in   ) :: dt,time
    real(kind=dp_t), intent(in   ) :: prob_lo(rho(1)%dim),prob_hi(rho(1)%dim) 
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local variables and array of multifabs
    type(multifab)                 :: rhonew(mla%nlevel)
    type(multifab)                 :: diff_fluxdiv(mla%nlevel)
    type(multifab)                 :: diff_fluxdivnew(mla%nlevel)
    type(multifab)                 :: rho_tot(mla%nlevel)        ! total density
    type(multifab)                 :: molarconc(mla%nlevel)      ! molar concentration
    type(multifab)                 :: molmtot(mla%nlevel)        ! total molar mass
    type(multifab)                 :: chi(mla%nlevel)            ! Chi-matrix
    type(multifab)                 :: D_MS(mla%nlevel)           ! D_MS-matrix
    type(multifab)                 :: Gama(mla%nlevel)           ! Gama-matrix
    type(multifab)                 :: stoch_fluxdiv(mla%nlevel)  ! stochastic fluxdiv
    type(multifab),  allocatable   :: stoch_W_fc(:,:,:)          ! WA and WB (nlevs,dim,n_rngs) 
    real(kind=dp_t), allocatable   :: weights(:)                 ! weights for stoch-time-integrators       
    real(kind=dp_t)                :: stage_time
    integer                        :: n,nlevs,i,dm,rng,n_rngs

    nlevs = mla%nlevel  ! number of levels 
    dm    = mla%dim     ! number of dimensions

    ! build cell-centered multifabs for nspecies and ghost cells contained in rho.
    ! fluxdiv needs no ghost cells as all computations will be in box interior. 
    ! rho_tot is cell-cented with ghost cells that many rho owns.
    do n=1,nlevs
       call multifab_build(rhonew(n),          mla%la(n), nspecies,    rho(n)%ng)
       call multifab_build(rho_tot(n),         mla%la(n), 1,           rho(n)%ng) 
       call multifab_build(diff_fluxdiv(n),    mla%la(n), nspecies,    0) 
       call multifab_build(diff_fluxdivnew(n), mla%la(n), nspecies,    0)
       call multifab_build(molmtot(n),         mla%la(n), 1,           rho(n)%ng)  
       call multifab_build(molarconc(n),       mla%la(n), nspecies,    rho(n)%ng) 
       call multifab_build(chi(n),             mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(D_MS(n),            mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(Gama(n),            mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(stoch_fluxdiv(n),   mla%la(n), nspecies,    0) 
    end do

    !========================================================
    ! initialize random number generator for stochastic flux
    !========================================================
 
    if (timeinteg_type .eq. 1) then       ! Euler-Maruyama
        n_rngs=1
    else if (timeinteg_type .eq. 2) then  ! Trapezoidal
        n_rngs=1
    else                                ! Midpoint & RK3
        n_rngs=2  
    end if
 
    ! allocate the multifabs
    allocate(stoch_W_fc(mla%nlevel,mla%dim,1:n_rngs))
    allocate(weights(n_rngs))
    
    if(use_stoch) then
      ! build the multifabs
       do n=1, nlevs 
          do rng=1, n_rngs 
             do i = 1,dm
                call multifab_build_edge(stoch_W_fc(n,i,rng),mla%la(n),nspecies,0,i)
             end do
          end do
       end do
    
      ! initialize stochastic flux on every face W(0,1) 
      call generate_random_increments(mla,n_rngs,stoch_W_fc)

    endif

   !==================================================================================
    select case(timeinteg_type)
   !==================================================================================
 
      case(1)
      !==================================================
      ! Euler time update 
      ! rho(t+dt)  =rho(t) + dt*fluxdiv(t) + sqrt(dt)K*W(t)
      !==================================================
      
      stage_time = time  
      if(use_stoch) weights(1) = 1.0d0 
      
      ! compute the total div of flux from rho
      call compute_mass_fluxdiv(mla,rho,rho_tot,molarconc,molmtot,molmass,chi,Gama,D_MS,&
                           diff_fluxdiv,stoch_fluxdiv,stoch_W_fc,Temp,dt,&
                           stage_time,dx,prob_lo,prob_hi,weights,n_rngs,the_bc_level)

      ! compute rho(t+dt) (only interior) 
      do n=1,nlevs
                       call saxpy(rho(n),-dt,diff_fluxdiv(n))
         if(use_stoch) call saxpy(rho(n),-dt,stoch_fluxdiv(n))
      end do 
  
      case(2)
      !=========================================================================================
      ! Heun's method: Predictor-Corrector explicit method 
      ! rhonew(t+dt) = rho(t)+ dt*fluxdiv(t,rho)                                 +sqrt(dt)K*W(t) 
      !    rho(t+dt) = rho(t)+(dt/2)*[fluxdiv(t,rho)+fluxdivnew(t+1,rhonew(t+1))]+sqrt(dt)K*W(t)  
      !=========================================================================================
      
      ! store old rho in rhonew after clearing up memory 
      do n=1,nlevs
          call setval(rhonew(n), 0.d0, all=.true.)
      end do
 
      do n=1,nlevs
         call saxpy(rhonew(n),1.0d0,rho(n))
      end do 
      
      !=====================
      ! Euler Predictor step
      !===================== 
      
      stage_time = time 
      if(use_stoch) weights(1) = 1.0d0 
      
      ! compute the total div of flux from rho
      call compute_mass_fluxdiv(mla,rho,rho_tot,molarconc,molmtot,molmass,chi,Gama,D_MS,&
                           diff_fluxdiv,stoch_fluxdiv,stoch_W_fc,Temp,dt,&
                           stage_time,dx,prob_lo,prob_hi,weights,n_rngs,the_bc_level)
      
      ! compute rhonew(t+dt) (only interior) 
      do n=1,nlevs
         call saxpy(rhonew(n),-dt,diff_fluxdiv(n))
         if(use_stoch) call saxpy(rhonew(n),-dt,stoch_fluxdiv(n))
      end do 
   
      ! update values of the ghost cells of rhonew
      do n=1,nlevs
         call multifab_fill_boundary(rhonew(n))
        
         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc(rhonew(n),1,rho_part_bc_comp,nspecies,the_bc_level(n),dx(n,:))
      end do

      !=========================== 
      ! Trapezoidal Corrector step
      !===========================
      
      stage_time = time + dt  
      
      ! compute the total div of flux from rho
      call compute_mass_fluxdiv(mla,rhonew,rho_tot,molarconc,molmtot,molmass,chi,Gama,D_MS,&
                           diff_fluxdivnew,stoch_fluxdiv,stoch_W_fc,Temp,dt,&
                           stage_time,dx,prob_lo,prob_hi,weights,n_rngs,the_bc_level)

      ! compute rho(t+dt) (only interior)
      do n=1,nlevs
                       call saxpy(rho(n),-0.5d0*dt,diff_fluxdiv(n))
                       call saxpy(rho(n),-0.5d0*dt,diff_fluxdivnew(n))
         if(use_stoch) call saxpy(rho(n),      -dt,stoch_fluxdiv(n))
      end do
 
      case(3)
      !========================================================================================
      ! Midpoint method (2-stage) 
      ! rhonew(t+dt/2)=rho(t)+(dt/2)*fluxdiv(t,y)                   + sqrt(dt/2)K*W1(t)                
      !    rho(t+dt)  =rho(t)+ dt*fluxdivnew[t+dt/2,rhonew(t+dt/2)] + sqrt(dt/2)K*(W1(t)+W2(t))
      !========================================================================================
      
      ! store old rho in rhonew after clearing up memory 
      do n=1,nlevs
          call setval(rhonew(n), 0.d0, all=.true.)
      end do
 
      do n=1,nlevs
         call saxpy(rhonew(n),1.0d0,rho(n))
      end do 

      !===========
      ! 1st stage
      !===========

      stage_time = time
      weights(1) = sqrt(0.5d0) 
      weights(2) = 0.0d0 
      
      ! compute the total div of flux from rho
      call compute_mass_fluxdiv(mla,rho,rho_tot,molarconc,molmtot,molmass,chi,Gama,D_MS,&
                           diff_fluxdiv,stoch_fluxdiv,stoch_W_fc,Temp,dt,&
                           stage_time,dx,prob_lo,prob_hi,weights,n_rngs,the_bc_level)
 
      ! compute rhonew(t+dt/2) (only interior) 
      do n=1,nlevs
                       call saxpy(rhonew(n),-0.5d0*dt,diff_fluxdiv(n))
         if(use_stoch) call saxpy(rhonew(n),      -dt,stoch_fluxdiv(n))
      end do 

      ! update values of the ghost cells of rhonew
      do n=1,nlevs
         call multifab_fill_boundary(rhonew(n))
        
         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc(rhonew(n),1,rho_part_bc_comp,nspecies,the_bc_level(n),dx(n,:))
      end do

      !===========
      ! 2nd stage
      !===========

      stage_time = time + dt/2.0d0
      weights(1) = sqrt(0.5d0) 
      weights(2) = sqrt(0.5d0) 

      ! compute the total div of flux from rho
      call compute_mass_fluxdiv(mla,rhonew,rho_tot,molarconc,molmtot,molmass,chi,Gama,D_MS,&
                           diff_fluxdivnew,stoch_fluxdiv,stoch_W_fc,Temp,dt,&
                           stage_time,dx,prob_lo,prob_hi,weights,n_rngs,the_bc_level)
 
      ! compute rho(t+dt) (only interior) 
      do n=1,nlevs
                       call saxpy(rho(n), -dt, diff_fluxdivnew(n))
         if(use_stoch) call saxpy(rho(n), -dt, stoch_fluxdiv(n))
      end do 
      
      case(4)
      !=======================================================================================================================
      ! 3rd order Runge-Kutta method 
      ! rhonew(t+dt)  =     rho(t)+                       dt*fluxdiv(t)                        + sqrt(dt)K*(a1*W1(t)+b1*W2(t)   
      ! rhonew(t+dt/2)= 3/4*rho(t)+ 1/4*[rhonew(t+dt)   + dt*fluxdivnew(t+dt,rhonew(t+dt))     + sqrt(dt)K*(a2*W1(t)+b2*W2(t)] 
      !    rho(t+dt)  = 1/3*rho(t)+ 2/3*[rhonew(t+dt/2) + dt*fluxdivnew(t+dt/2,rhonew(t+dt/2)) + sqrt(dt)K*(a3*W1(t)+b3*W2(t)] 
      !=======================================================================================================================

      ! store old rho in rhonew after clearing up memory 
      do n=1,nlevs
          call setval(rhonew(n), 0.d0, all=.true.)
      end do
 
      do n=1,nlevs
         call saxpy(rhonew(n),1.0d0,rho(n))
      end do 

      !===========
      ! 1st stage
      !===========

      stage_time = time 
      weights(1) = 1.0d0 
      weights(2) = (2*sqrt(2.0d0)+sqrt(3.0d0))/5.0d0 
      
      ! compute the total div of flux from rho
      call compute_mass_fluxdiv(mla,rho,rho_tot,molarconc,molmtot,molmass,chi,Gama,D_MS,&
                           diff_fluxdiv,stoch_fluxdiv,stoch_W_fc,Temp,dt,&
                           stage_time,dx,prob_lo,prob_hi,weights,n_rngs,the_bc_level)
 
      ! compute rhonew(t+dt) (only interior) 
      do n=1,nlevs
                       call saxpy(rhonew(n), -dt, diff_fluxdiv(n))
         if(use_stoch) call saxpy(rhonew(n), -dt, stoch_fluxdiv(n))
      end do 
   
      ! update values of the ghost cells of rhonew
      do n=1,nlevs
         call multifab_fill_boundary(rhonew(n))
        
         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc(rhonew(n),1,rho_part_bc_comp,nspecies,the_bc_level(n),dx(n,:))
      end do

      !===========
      ! 2nd stage
      !===========

      stage_time = time + dt
      weights(1) = 0.25d0 
      weights(2) = (-4*sqrt(2.0d0)+3*sqrt(3.0d0))/5.0d0 

      ! compute the total div of flux from rho
      call compute_mass_fluxdiv(mla,rhonew,rho_tot,molarconc,molmtot,molmass,chi,Gama,D_MS,&
                           diff_fluxdivnew,stoch_fluxdiv,stoch_W_fc,Temp,dt,&
                           stage_time,dx,prob_lo,prob_hi,weights,n_rngs,the_bc_level)


      ! compute rhonew(t+dt/2) (only interior) 
      do n=1,nlevs
                       call multifab_mult_mult_s(rhonew(n),0.25d0)
                       call saxpy(rhonew(n),  0.75d0   , rho(n))
                       call saxpy(rhonew(n), -0.25d0*dt, diff_fluxdivnew(n))
         if(use_stoch) call saxpy(rhonew(n), -0.25d0*dt, stoch_fluxdiv(n))
      end do 

      ! update values of the ghost cells of rho_star
      do n=1,nlevs
         call multifab_fill_boundary(rhonew(n))
        
         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc(rhonew(n),1,rho_part_bc_comp,nspecies,the_bc_level(n),dx(n,:))
      end do
      
      !===========
      ! 3rd stage
      !===========

      ! free up the values of fluxdivnew
      do n=1,nlevs
         call setval(diff_fluxdivnew(n), 0.d0, all=.true.)
      end do 

      stage_time = time + dt/2.0d0
      weights(1) = 2.0d0/3.0d0 
      weights(2) = (sqrt(2.0d0)-2*sqrt(3.0d0))/10.0d0

      ! compute the total div of flux from rho
      call compute_mass_fluxdiv(mla,rhonew,rho_tot,molarconc,molmtot,molmass,chi,Gama,D_MS,&
                           diff_fluxdivnew,stoch_fluxdiv,stoch_W_fc,Temp,dt,&
                           stage_time,dx,prob_lo,prob_hi,weights,n_rngs,the_bc_level)

      ! compute rho(t+dt) (only interior) 
      do n=1,nlevs
                       call multifab_mult_mult_s(rho(n),(1.0d0/3.0d0))
                       call saxpy(rho(n),  (2.0d0/3.0d0)   , rhonew(n))
                       call saxpy(rho(n), -(2.0d0/3.0d0)*dt, diff_fluxdivnew(n))
         if(use_stoch) call saxpy(rho(n), -(2.0d0/3.0d0)*dt, stoch_fluxdiv(n))
      end do 

    !=============================================================================================
    end select  
    !=============================================================================================

    ! fill the ghost cells of rho
    do n=1, nlevs
       call multifab_fill_boundary(rho(n))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(rho(n),1,rho_part_bc_comp,nspecies,the_bc_level(n),dx(n,:))
    end do   

    ! free the multifab allocated memory
    do n=1,nlevs
       call multifab_destroy(rhonew(n))
       call multifab_destroy(diff_fluxdiv(n))
       call multifab_destroy(diff_fluxdivnew(n))
       call multifab_destroy(rho_tot(n))
       call multifab_destroy(molarconc(n))
       call multifab_destroy(molmtot(n))
       call multifab_destroy(chi(n))
       call multifab_destroy(D_MS(n))
       call multifab_destroy(Gama(n))
       call multifab_destroy(stoch_fluxdiv(n))
    end do
    
    if(use_stoch) then
       call destroy_random_increments(mla,n_rngs,stoch_W_fc)
    endif
    deallocate(stoch_W_fc)
    deallocate(weights)
  
  end subroutine advance

end module advance_module 
