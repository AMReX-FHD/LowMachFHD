module advance_module

  use multifab_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use ml_layout_module
  use diffusive_fluxdiv_module
  use probin_multispecies_module

  implicit none

  private

  public :: advance

contains

  subroutine advance(mla,rho,Dbar,Gama,mass,dx,dt,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rho(:)
    real(kind=dp_t), intent(in   ) :: Dbar(:,:)
    real(kind=dp_t), intent(in   ) :: Gama(:,:)
    real(kind=dp_t), intent(in   ) :: mass(:) 
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(in   ) :: dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local variables
    type(multifab) :: rho_prev(mla%nlevel), rho_prev1(mla%nlevel)
    type(multifab) :: rho_intm(mla%nlevel), rho_intm1(mla%nlevel)
    type(multifab) :: fluxdiv(mla%nlevel)
    type(multifab) :: fluxdiv_prev(mla%nlevel)
    integer        :: n, nlevs
    
    nlevs = mla%nlevel  ! number of levels 

    ! build cell-centered multifabs for nspecies and 0 ghost cells
    do n=1,nlevs
       call multifab_build(rho_prev(n),mla%la(n),nspecies,0)
       call multifab_build(rho_prev1(n),mla%la(n),nspecies,0)
       call multifab_build(rho_intm(n),mla%la(n),nspecies,0)
       call multifab_build(rho_intm1(n),mla%la(n),nspecies,0)
       call multifab_build(fluxdiv(n),mla%la(n),nspecies,0)
       call multifab_build(fluxdiv_prev(n),mla%la(n),nspecies,0)
    enddo
   
    ! copy rho in fluxdiv for nspecies and 0 ghost cells 
    do n=1,nlevs
       call multifab_copy_c(fluxdiv(n),1,rho(n),1,nspecies)
    end do 
  
    !do n=1,nlevs
    !   call multifab_fill_boundary(fluxdiv(n))
    !end do 
    
    ! compute div-of-flux from rho after filling ghosts also. 
    call diffusive_fluxdiv(mla,fluxdiv,Dbar,Gama,mass,dx,dt,the_bc_level)
 
    select case(timeinteg_type)
 
      ! Donev: Missing calls to fill_boundary for intermediate results
      ! Before I do compute_div etc I need to fill_boundary  
      case(1)
      !==================================================================================
      ! Euler time update 
      !==================================================================================
      
      ! multiply div-of-flux with dt, starting from component 1 with 0 ghost-cell
      do n=1,nlevs
         call multifab_mult_mult_s_c(fluxdiv(n),1,dt,nspecies,0) 
      end do 
    
      ! add this to rho to advance in time
      do n=1,nlevs
         call multifab_plus_plus(rho(n),fluxdiv(n))            
      end do 

      case(2)
      !====================================================================================
      ! Heun's method: Predictor-Corrector explicit method 
      ! rho_predictor(t+1)=rho(t)+dt*fluxdiv(t,rho) 
      ! rho_corrector(t+1)=rho(t)+(dt/2)*[fluxdiv(t,rho) + fluxdiv(t+1,rho_predictor(t+1))] 
      !====================================================================================
      
      ! store the previous values for nspecies and 0 ghost cells
      do n=1,nlevs
         call multifab_copy_c(rho_prev(n),1,rho(n),1,nspecies)
         call multifab_copy_c(fluxdiv_prev(n),1,fluxdiv(n),1,nspecies)
      end do 
   
      !=====================
      ! Euler Predictor step
      !===================== 
      
      do n=1,nlevs
         call multifab_mult_mult_s_c(fluxdiv(n),1,dt,nspecies,0) 
      end do 

      do n=1,nlevs
         call multifab_plus_plus(rho_prev(n),fluxdiv(n))            
      end do 

      ! ghost cells are needed in rho_prev; always call before calling diffusive_fluxdiv
      !do n=1,nlevs
      !   call multifab_fill_boundary(rho_prev(n))
      !enddo
 
      ! compute div-of-flux with Predictor results 
      call diffusive_fluxdiv(mla,rho_prev,Dbar,Gama,mass,dx,dt,the_bc_level)
      
      !===================== 
      ! Corrector step
      !=====================
 
      ! New fluxdiv is stored in returned quantity rho_prev 
      do n=1,nlevs
         call multifab_plus_plus(fluxdiv_prev(n),rho_prev(n))            
      end do 
    
      ! multiply this with dt/2, starting from component 1 with 0 ghost-cell
      do n=1,nlevs
         call multifab_mult_mult_s_c(fluxdiv_prev(n),1,dt*0.5d0,nspecies,0) 
      end do 
    
      ! add this to rho to advance in time
      do n=1,nlevs
         call multifab_plus_plus(rho(n),fluxdiv_prev(n))            
      end do 
    
      case(3)
      !=====================================================================================
      ! Midpoint method (2-stage) 
      ! rho_prev(t+1/2)=rho(t)+(dt/2)*fluxdiv(t,y)                
      ! rho_new(t+1)=rho(t)+ dt*fluxdiv[t+1/2,rho_prev(t+1/2)]  
      !=====================================================================================
      
      ! store the previous values of rho
      do n=1,nlevs
         call multifab_copy_c(rho_prev(n),1,rho(n),1,nspecies)
      end do  
    
      ! multiply fluxdiv with dt/2 
      do n=1,nlevs
         call multifab_mult_mult_s_c(fluxdiv(n),1,dt*0.5d0,nspecies,0) 
      end do 
   
      ! get rho at t+1/2 
      do n=1,nlevs
         call multifab_plus_plus(rho_prev(n),fluxdiv(n))            
      end do

      ! ghost cells are needed in rho_prev; always call before calling diffusive_fluxdiv
      !do n=1,nlevs
      !   call multifab_fill_boundary(rho_prev(n))
      !enddo
      
      call diffusive_fluxdiv(mla,rho_prev,Dbar,Gama,mass,dx,dt*0.5d0,the_bc_level)

      ! multiply that with dt, starting from component 1 with 0 ghost-cell
      do n=1,nlevs
         call multifab_mult_mult_s_c(rho_prev(n),1,dt,nspecies,0) 
      end do 
    
      ! add this to rho to advance in time
      do n=1,nlevs
         call multifab_plus_plus(rho(n),rho_prev(n))            
      end do  

      case(4)
      !=======================================================================================
      ! 3rd order Runge-Kutta method 
      ! rho_prev(t+1)  =     rho(t)+                        dt*fluxdiv(t)                      
      ! rho_intm(t+1/2)= 3/4*rho(t)+ 1/4*[rho_prev(t+1)   + dt*fluxdiv(t+1,  rho_prev(t+1))] 
      ! rho_new(t+1)   = 1/3*rho(t)+ 2/3*[rho_intm(t+1/2) + dt*fluxdiv(t+1/2,rho_intm(t+1/2))] 
      !=======================================================================================

      !===========
      ! 1st stage
      !===========
      do n=1,nlevs
         ! store the previous values for 2nd and 3rd stage
         call multifab_copy_c(rho_prev(n),1,rho(n),1,nspecies)
         call multifab_copy_c(rho_prev1(n),1,rho(n),1,nspecies)
      end do 
 
      do n=1,nlevs
         ! multiply div-of-flux with dt, starting from component 1 with 0 ghost-cell
         call multifab_mult_mult_s_c(fluxdiv(n),1,dt,nspecies,0) 
      end do 
    
      do n=1,nlevs
         ! add this to rho 
         call multifab_plus_plus(rho_prev(n),fluxdiv(n))            
      end do 

      do n=1,nlevs
         ! store this value of rho_prev for 2nd stage
         call multifab_copy_c(rho_intm(n),1,rho_prev(n),1,nspecies)
      end do 

      ! ghost cells are needed in rho_prev; always call before calling diffusive_fluxdiv
      do n=1,nlevs
         call multifab_fill_boundary(rho_prev(n))
      enddo
      
      ! ghost cells are needed in rho_prev; always call before calling diffusive_fluxdiv
      !do n=1,nlevs
      !   call multifab_fill_boundary(rho_prev(n))
      !enddo
 
      !===========
      ! 2nd stage
      !===========
      call diffusive_fluxdiv(mla,rho_prev,Dbar,Gama,mass,dx,dt,the_bc_level)

      ! multiply fluxdiv with dt/4
      do n=1,nlevs
         call multifab_mult_mult_s_c(rho_prev(n),1,0.25d0*dt,nspecies,0) 
      end do 

      ! multiply previous rho with 1/4
      do n=1,nlevs
         call multifab_mult_mult_s_c(rho_intm(n),1,0.25d0,nspecies,0) 
      end do 
    
      ! add this two 
      do n=1,nlevs
         call multifab_plus_plus(rho_intm(n),rho_prev(n))            
      end do 

      ! multiply first saved rho with 3/4
      do n=1,nlevs
         call multifab_mult_mult_s_c(rho_prev1(n),1,0.75d0,nspecies,0) 
      end do 
    
      ! add them all
      do n=1,nlevs
         call multifab_plus_plus(rho_intm(n),rho_prev1(n))            
      end do 

      ! store this intermediate value of rho for 3rd stage
      do n=1,nlevs
         call multifab_copy_c(rho_intm1(n),1,rho_intm(n),1,nspecies)
      end do 
  
      ! ghost cells are needed in rho_intm; always call before calling diffusive_fluxdiv
      !do n=1,nlevs
      !   call multifab_fill_boundary(rho_intm(n))
      !enddo
      
      !===========
      ! 3rd stage
      !===========
      call diffusive_fluxdiv(mla,rho_intm,Dbar,Gama,mass,dx,dt*0.5d0,the_bc_level)

      ! multiply div-of-flux with 2*dt/3
      do n=1,nlevs
         call multifab_mult_mult_s_c(rho_intm(n),1,(2.0d0/3.0d0)*dt,nspecies,0) 
      end do 
    
      ! multiply stored intermediate rho with 2/3 
      do n=1,nlevs
         call multifab_mult_mult_s_c(rho_intm1(n),1,(2.0d0/3.0d0),nspecies,0) 
      end do 

      ! add this to rho to advance in time
      do n=1,nlevs
         call multifab_plus_plus(rho_intm(n),rho_intm1(n))            
      end do 
     
      ! multiply oldest rho with 1/3 
      do n=1,nlevs
         call multifab_mult_mult_s_c(rho(n),1,(1.0d0/3.0d0),nspecies,0) 
      end do 
      
      ! add this to rho to advance in time
      do n=1,nlevs
         call multifab_plus_plus(rho(n),rho_intm(n))            
      end do 

    end select  

    do n=1, nlevs
       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(rho(n))

       ! fill non-periodic domain boundary ghost cells
    !   call multifab_physbc(rho(n),1,rho_part_bc_comp,nspecies,the_bc_level(n), & 
    !                        dx(n,:),.false.)
    end do   

    ! free the multifab allocated memory
    do n=1,nlevs
       call multifab_destroy(rho_prev(n))
       call multifab_destroy(rho_prev1(n))
       call multifab_destroy(rho_intm(n))
       call multifab_destroy(rho_intm1(n))
       call multifab_destroy(fluxdiv(n))
       call multifab_destroy(fluxdiv_prev(n))
    end do

  end subroutine advance

end module advance_module 
