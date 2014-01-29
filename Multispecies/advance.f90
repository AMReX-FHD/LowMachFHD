module advance_module

  use multifab_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use ml_layout_module
  use fluid_model_module
  use external_force_module
  use diffusive_fluxdiv_module
  use convert_variables_module
  use probin_multispecies_module

  implicit none

  private

  public :: advance

contains

  subroutine advance(mla,rho,dx,dt,time,prob_lo,prob_hi,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rho(:)
    ! For now that routine can just set them equal to Dbar_in and Gama_in (i.e., constant)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(in   ) :: dt,time
    real(kind=dp_t), intent(in   ) :: prob_lo(rho(1)%dim),prob_hi(rho(1)%dim) 
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local variables and array of multifabs
    type(multifab)               :: rhonew(mla%nlevel),fluxdiv(mla%nlevel),fluxdivnew(mla%nlevel)
    type(multifab)               :: rho_tot(mla%nlevel)     ! total density
    type(multifab)               :: molarconc(mla%nlevel)   ! molar concentration
    type(multifab)               :: molmtot(mla%nlevel)     ! total molar mass
    type(multifab)               :: chi(mla%nlevel)         ! Chi-matrix
    type(multifab)               :: D_MS(mla%nlevel)        ! D_MS-matrix
    type(multifab)               :: Gama(mla%nlevel)        ! Gama-matrix
    type(multifab)               :: rhoWchiGama(mla%nlevel) ! rho*W*chi*Gama
    integer                      :: n,nlevs
    real(kind=dp_t)              :: stage_time  
    real(kind=dp_t), allocatable :: molmass(:) 
    
    nlevs = mla%nlevel  ! number of levels 
 
    ! build cell-centered multifabs for nspecies and ghost cells contained in rho.
    ! fluxdiv needs no ghost cells as all computations will be in box interior. 
    ! rho_tot is cell-cented with ghost cells that many rho owns.
    do n=1,nlevs
       call multifab_build(rhonew(n),      mla%la(n), nspecies,    rho(n)%ng)
       call multifab_build(rho_tot(n),     mla%la(n), 1,           rho(n)%ng) 
       call multifab_build(fluxdiv(n),     mla%la(n), nspecies,    0) 
       call multifab_build(fluxdivnew(n),  mla%la(n), nspecies,    0)
       call multifab_build(molmtot(n),     mla%la(n), 1,           rho(n)%ng)  
       call multifab_build(molarconc(n),   mla%la(n), nspecies,    rho(n)%ng) 
       call multifab_build(chi(n),         mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(D_MS(n),        mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(Gama(n),        mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(rhoWchiGama(n), mla%la(n), nspecies**2, rho(n)%ng)
    enddo
    allocate(molmass(nspecies))
   
    ! free up memory 
    do n=1,nlevs
       call setval(rhonew(n),      0.d0, all=.true.)
       call setval(rho_tot(n),     0.d0, all=.true.)
       call setval(fluxdiv(n),     0.d0, all=.true.)
       call setval(fluxdivnew(n),  0.d0, all=.true.)
       call setval(molmtot(n),     0.d0, all=.true.)
       call setval(molarconc(n),   0.d0, all=.true.)
       call setval(chi(n),         0.d0, all=.true.)
       call setval(D_MS(n),        0.d0, all=.true.)
       call setval(Gama(n),        0.d0, all=.true.)
       call setval(rhoWchiGama(n), 0.d0, all=.true.)
    enddo 
    molmass(1:nspecies) = 1.0d0  

   !==================================================================================
    select case(timeinteg_type)
   !==================================================================================
 
      case(1)
      !===================================
      ! Euler time update 
      ! rho(t+dt)  =rho(t) + dt*fluxdiv(t) 
      !===================================
      
      stage_time = time   

      ! compute molmtot,molarconc & rho_tot (primitive variables) for each-cell from rho(conserved) 
      call convert_cons_to_prim(mla,rho,rho_tot,molarconc,molmtot,molmass,the_bc_level)
      
      ! populate D_MS, Gama and molmass 
      call fluid_model(mla,rho,rho_tot,molarconc,molmtot,molmass,D_MS,Gama,molmass,the_bc_level) 

      ! compute chi and rho*W*chi*Gama
      call compute_chi(mla,rho,rho_tot,molarconc,molmtot,molmass,chi,rhoWchiGama,D_MS,Gama,the_bc_level)
      
      ! compute chi and rho*W*chi*Gama
      call compute_rhoWchiGama(mla,rho,rho_tot,chi,Gama,rhoWchiGama,the_bc_level)
 
      ! compute fluxdiv that contain results in interior only while rho contains ghost values filled in 
      ! init or end of this code
      call diffusive_fluxdiv(mla,rho,fluxdiv,molarconc,rhoWchiGama,molmass,dx,the_bc_level)

      ! compute external forcing for manufactured solution and add to fluxdiv
      ! This routine should not take D_MS as input
      ! It should only take primitive variables as input
      call external_source(mla,rho,fluxdiv,prob_lo,prob_hi,dx,stage_time)
      
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
      
      ! compute fluxdiv 
      call diffusive_fluxdiv(mla,rho,fluxdiv,molarconc,rhoWchiGama,molmass,dx,the_bc_level)
      
      ! compute external forcing for manufactured solution and add to fluxdiv
      call external_source(mla,rho,fluxdiv,prob_lo,prob_hi,dx,stage_time)
      
      ! compute rhonew(t+dt) (only interior) 
      do n=1,nlevs
         call saxpy(rhonew(n),dt,fluxdiv(n))
      enddo 
   
      ! update values of the ghost cells of rhonew
      do n=1,nlevs
         call multifab_fill_boundary(rhonew(n))
        
         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc(rhonew(n),1,rho_part_bc_comp,nspecies,the_bc_level(n),dx(n,:),.false.)
      enddo

      !=========================== 
      ! Trapezoidal Corrector step
      !===========================
      
      ! Donev:
      ! This code is wrong here
      ! The values of 
      ! molarconc,molmtot,chi,rhoWchiGama,D_MS,Gama
      ! need to be updated here because the rho's changed
      ! So you need to repeat the same computation you did at the beginning
      ! again here, passing rhonew instead of rho
      ! For example, you need to call compute_coefficient again!
      ! Similatly in all other temporal schemes. For RK3 you need to do this one more time in 3rd stage
      ! IF THIS IS NOT CLEAR MAKE SURE TO TALK TO ME -- IT IS CRUCIAL YOU UNDERSTAND THIS POINT

      ! compute fluxdiv(t+1,rhonew(t+1)) 
      call diffusive_fluxdiv(mla,rhonew,fluxdivnew,molarconc,rhoWchiGama,molmass,dx,the_bc_level)

      ! compute external forcing for manufactured solution and add to fluxdiv
      stage_time = time + dt  
      call external_source(mla,rhonew,fluxdivnew,prob_lo,prob_hi,dx,stage_time)
      
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
 
      ! compute fluxdiv(t) from rho(t); (interior only) 
      call diffusive_fluxdiv(mla,rho,fluxdiv,molarconc,rhoWchiGama,molmass,dx,the_bc_level)
      
      ! compute external forcing for manufactured solution and add to fluxdiv
      stage_time = time
      call external_source(mla,rho,fluxdiv,prob_lo,prob_hi,dx,stage_time)
      
      ! compute rhonew(t+dt/2) (only interior) 
      do n=1,nlevs
         call saxpy(rhonew(n),0.5d0*dt,fluxdiv(n))
      enddo 

      ! update values of the ghost cells of rhonew
      do n=1,nlevs
         call multifab_fill_boundary(rhonew(n))
        
         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc(rhonew(n),1,rho_part_bc_comp,nspecies,the_bc_level(n),dx(n,:),.false.)
      enddo

      ! compute new div-of-flux 
      call diffusive_fluxdiv(mla,rhonew,fluxdivnew,molarconc,rhoWchiGama,molmass,dx,the_bc_level)

      ! compute external forcing for manufactured solution and add to fluxdiv
      stage_time = time + dt/2.0d0
      call external_source(mla,rhonew,fluxdivnew,prob_lo,prob_hi,dx,stage_time)
      
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

      ! compute fluxdiv(t) from rho(t) (interior only) 
      call diffusive_fluxdiv(mla,rho,fluxdiv,molarconc,rhoWchiGama,molmass,dx,the_bc_level)
      
      ! compute external forcing for manufactured solution and add to fluxdiv
      stage_time = time 
      call external_source(mla,rho,fluxdiv,prob_lo,prob_hi,dx,stage_time)
      
      ! compute rhonew(t+dt) (only interior) 
      do n=1,nlevs
         call saxpy(rhonew(n),dt,fluxdiv(n))
      enddo 
   
      ! update values of the ghost cells of rhonew
      do n=1,nlevs
         call multifab_fill_boundary(rhonew(n))
        
         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc(rhonew(n),1,rho_part_bc_comp,nspecies,the_bc_level(n),dx(n,:),.false.)
      enddo

      ! compute fluxdivnew(t+dt,rhonew(t+dt))
      call diffusive_fluxdiv(mla,rhonew,fluxdivnew,molarconc,rhoWchiGama,molmass,dx,the_bc_level)

      ! compute external forcing for manufactured solution and add to fluxdiv
      stage_time = time + dt
      call external_source(mla,rhonew,fluxdivnew,prob_lo,prob_hi,dx,stage_time)
      
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
         call multifab_physbc(rhonew(n),1,rho_part_bc_comp,nspecies,the_bc_level(n),dx(n,:),.false.)
      enddo
      
      ! free up the values of fluxdivnew
      do n=1,nlevs
         call setval(fluxdivnew(n), 0.d0, all=.true.)
      enddo 
 
      ! compute fluxdiv(t+dt/2,rhonew(t+dt/2)) 
      call diffusive_fluxdiv(mla,rhonew,fluxdivnew,molarconc,rhoWchiGama,molmass,dx,the_bc_level)

      ! compute external forcing for manufactured solution and add to fluxdiv
      stage_time = time + dt/2.0d0
      call external_source(mla,rhonew,fluxdivnew,prob_lo,prob_hi,dx,stage_time)
      
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
       call multifab_physbc(rho(n),1,rho_part_bc_comp,nspecies,the_bc_level(n),dx(n,:),.false.)
    enddo   

    ! free the multifab allocated memory
    do n=1,nlevs
       call multifab_destroy(rhonew(n))
       call multifab_destroy(fluxdiv(n))
       call multifab_destroy(fluxdivnew(n))
       call multifab_destroy(rho_tot(n))
       call multifab_destroy(molmtot(n))
       call multifab_destroy(molarconc(n))
       call multifab_destroy(chi(n))
       call multifab_destroy(D_MS(n))
       call multifab_destroy(Gama(n))
       call multifab_destroy(rhoWchiGama(n))
    enddo
    deallocate(molmass)

  end subroutine advance

end module advance_module 
