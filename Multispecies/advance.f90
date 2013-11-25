module advance_module

  use multifab_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use ml_layout_module
  use external_force_module
  use diffusive_fluxdiv_module
  use probin_multispecies_module

  implicit none

  private

  public :: advance

contains

  subroutine advance(mla,rho,Dbar,Gama,mass,dx,dt,time,prob_lo,prob_hi,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rho(:)
    real(kind=dp_t), intent(in   ) :: Dbar(:,:)
    real(kind=dp_t), intent(in   ) :: Gama(:,:)
    real(kind=dp_t), intent(in   ) :: mass(:) 
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(in   ) :: dt,time
    real(kind=dp_t), intent(in   ) :: prob_lo(rho(1)%dim),prob_hi(rho(1)%dim) 
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local variables
    type(multifab)  :: rhonew(mla%nlevel),fluxdiv(mla%nlevel),fluxdivnew(mla%nlevel)
    integer         :: n,nlevs
    real(kind=dp_t) :: stage_time  
 
    nlevs = mla%nlevel  ! number of levels 
 
    ! build cell-centered multifabs for nspecies and ghost cells contained in rho.
    ! fluxdiv needs no ghost cells as all computations will be in box interior 
    do n=1,nlevs
       call multifab_build(rhonew(n),    mla%la(n),nspecies,rho(n)%ng)
       call multifab_build(fluxdiv(n),   mla%la(n),nspecies,0) 
       call multifab_build(fluxdivnew(n),mla%la(n),nspecies,0)
    enddo
   
    ! initialize new quantities to zero 
      do n=1,nlevs
         call setval(rhonew(n),    0.d0,all=.true.)
         call setval(fluxdiv(n),   0.d0,all=.true.)
         call setval(fluxdivnew(n),0.d0,all=.true.)
      end do 
 
   !==================================================================================
    select case(timeinteg_type)
   !==================================================================================
 
      case(1)
      !===================================
      ! Euler time update 
      ! rho(t+dt)  =rho(t) + dt*fluxdiv(t) 
      !===================================
      
      stage_time = time   
 
      ! compute fluxdiv; fluxdiv contain results in interior only, while rho contains 
      ! ghost values filled in init or end of this code
      call diffusive_fluxdiv(mla,rho,fluxdiv,Dbar,Gama,mass,dx,the_bc_level)

      ! compute external forcing for manufactured solution and add to fluxdiv
      call external_source(mla,rho,fluxdiv,prob_lo,prob_hi,dx,stage_time)
      
      ! compute rho(t+dt) (only interior) 
      do n=1,nlevs
         call saxpy(rho(n),-dt,fluxdiv(n))
      end do 
    
      case(2)
      !============================================================================
      ! Heun's method: Predictor-Corrector explicit method 
      ! rhonew(t+dt) = rho(t)+dt*fluxdiv(t,rho) 
      !    rho(t+dt) = rho(t)+(dt/2)*[fluxdiv(t,rho) + fluxdivnew(t+1,rhonew(t+1))] 
      !============================================================================
      
      ! store old rho in rhonew (rhonew previously set to zero) 
      do n=1,nlevs
         call saxpy(rhonew(n),1.0d0,rho(n))
      end do 
      
      !=====================
      ! Euler Predictor step
      !===================== 
      
      ! compute fluxdiv 
      call diffusive_fluxdiv(mla,rho,fluxdiv,Dbar,Gama,mass,dx,the_bc_level)
      
      ! compute external forcing for manufactured solution and add to fluxdiv
      call external_source(mla,rho,fluxdiv,prob_lo,prob_hi,dx,stage_time)
      
      ! compute rhonew(t+dt) (only interior) 
      do n=1,nlevs
         call saxpy(rhonew(n),-dt,fluxdiv(n))
      end do 
   
      ! update values of the ghost cells of rhonew
      do n=1,nlevs
         call multifab_fill_boundary(rhonew(n))
        
         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc(rhonew(n),1,rho_part_bc_comp,nspecies,the_bc_level(n),dx(n,:),.false.)
      enddo

      ! compute fluxdiv(t+1,rhonew(t+1)) 
      call diffusive_fluxdiv(mla,rhonew,fluxdivnew,Dbar,Gama,mass,dx,the_bc_level)

      ! compute external forcing for manufactured solution and add to fluxdiv
      stage_time = time + dt  
      call external_source(mla,rhonew,fluxdivnew,prob_lo,prob_hi,dx,stage_time)
      
      !=========================== 
      ! Trapezoidal Corrector step
      !===========================

      ! compute rho(t+dt) (only interior)
      do n=1,nlevs
         call saxpy(rho(n),-dt*0.5d0,fluxdiv(n))
         call saxpy(rho(n),-dt*0.5d0,fluxdivnew(n))
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
      end do 
 
      ! compute fluxdiv(t) from rho(t); (interior only) 
      call diffusive_fluxdiv(mla,rho,fluxdiv,Dbar,Gama,mass,dx,the_bc_level)
      
      ! compute external forcing for manufactured solution and add to fluxdiv
      stage_time = time
      call external_source(mla,rho,fluxdiv,prob_lo,prob_hi,dx,stage_time)
      
      ! compute rhonew(t+dt/2) (only interior) 
      do n=1,nlevs
         call saxpy(rhonew(n),-0.5d0*dt,fluxdiv(n))
      end do 

      ! update values of the ghost cells of rhonew
      do n=1,nlevs
         call multifab_fill_boundary(rhonew(n))
        
         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc(rhonew(n),1,rho_part_bc_comp,nspecies,the_bc_level(n),dx(n,:),.false.)
      enddo

      ! compute new div-of-flux 
      call diffusive_fluxdiv(mla,rhonew,fluxdivnew,Dbar,Gama,mass,dx,the_bc_level)

      ! compute external forcing for manufactured solution and add to fluxdiv
      stage_time = time + dt/2.0d0
      call external_source(mla,rhonew,fluxdivnew,prob_lo,prob_hi,dx,stage_time)
      
      ! compute rho(t+dt) (only interior) 
      do n=1,nlevs
         call saxpy(rho(n),-dt,fluxdivnew(n))
      end do 
      
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
      end do 
      
      !===========
      ! 1st stage
      !===========

      ! compute fluxdiv(t) from rho(t) (interior only) 
      call diffusive_fluxdiv(mla,rho,fluxdiv,Dbar,Gama,mass,dx,the_bc_level)
      
      ! compute external forcing for manufactured solution and add to fluxdiv
      stage_time = time 
      call external_source(mla,rho,fluxdiv,prob_lo,prob_hi,dx,stage_time)
      
      ! compute rhonew(t+dt) (only interior) 
      do n=1,nlevs
         call saxpy(rhonew(n),-dt,fluxdiv(n))
      end do 
   
      ! update values of the ghost cells of rhonew
      do n=1,nlevs
         call multifab_fill_boundary(rhonew(n))
        
         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc(rhonew(n),1,rho_part_bc_comp,nspecies,the_bc_level(n),dx(n,:),.false.)
      enddo

      ! compute fluxdivnew(t+dt,rhonew(t+dt))
      call diffusive_fluxdiv(mla,rhonew,fluxdivnew,Dbar,Gama,mass,dx,the_bc_level)

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
         call saxpy(rhonew(n),-0.25d0*dt,fluxdivnew(n))
      end do 

      ! update values of the ghost cells of rho_star
      do n=1,nlevs
         call multifab_fill_boundary(rhonew(n))
        
         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc(rhonew(n),1,rho_part_bc_comp,nspecies,the_bc_level(n),dx(n,:),.false.)
      enddo
      
      ! free up the values of fluxdivnew
      do n=1,nlevs
         call setval(fluxdivnew(n), 0.d0, all=.true.)
      end do 
 
      ! compute fluxdiv(t+dt/2,rhonew(t+dt/2)) 
      call diffusive_fluxdiv(mla,rhonew,fluxdivnew,Dbar,Gama,mass,dx,the_bc_level)

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
         call saxpy(rho(n),-(2.0d0/3.0d0)*dt,fluxdivnew(n))
      end do 

    !=============================================================================================
    end select  
    !=============================================================================================

    ! fill the ghost cells of rho
    do n=1, nlevs
       call multifab_fill_boundary(rho(n))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(rho(n),1,rho_part_bc_comp,nspecies,the_bc_level(n),dx(n,:),.false.)
    end do   

    ! free the multifab allocated memory
    do n=1,nlevs
       call multifab_destroy(rhonew(n))
       call multifab_destroy(fluxdiv(n))
       call multifab_destroy(fluxdivnew(n))
    end do

  end subroutine advance

end module advance_module 
