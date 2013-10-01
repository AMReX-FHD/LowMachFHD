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
 
    ! build cell-centered multifabs for nspecies and one ghost cell
    do n=1,nlevs
       call multifab_build(rho_prev(n),mla%la(n),nspecies,1)
       call multifab_build(rho_prev1(n),mla%la(n),nspecies,1)
       call multifab_build(rho_intm(n),mla%la(n),nspecies,1)
       call multifab_build(rho_intm1(n),mla%la(n),nspecies,1)
       call multifab_build(fluxdiv(n),mla%la(n),nspecies,1)
       call multifab_build(fluxdiv_prev(n),mla%la(n),nspecies,1)
    enddo
    
    ! copy rho in fluxdiv for nspecies including ghost cells of rho that has
    ! either been filled in init or at the end of this routine 
    do n=1,nlevs
       call multifab_copy_c(fluxdiv(n),1,rho(n),1,nspecies,fluxdiv(n)%ng)
    end do 
    
    ! compute div-of-flux from rho after filling ghost cells. output will
    ! contain the interior box only, so one has to fill boundary or manually 
    ! fill the ghost cells 
    call diffusive_fluxdiv(mla,fluxdiv,Dbar,Gama,mass,dx,dt,the_bc_level)
   
    ! update the ghost cells 
    do n=1,nlevs
       call multifab_fill_boundary(fluxdiv(n))
       
       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(fluxdiv(n),1,rho_part_bc_comp,nspecies,the_bc_level(n), & 
                            dx(n,:),.false.)
    end do 
 
    select case(timeinteg_type)
 
      case(1)
      !==================================================================================
      ! Euler time update 
      !==================================================================================
      !write(*,*), 'Using Euler explicit method'
      
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
      !write(*,*), 'Using Predictor corrector explicit method'
      
      ! store the previous values for nspecies and 0 ghost cells
      do n=1,nlevs
         call multifab_copy_c(rho_prev(n),1,rho(n),1,nspecies,rho_prev(n)%ng)
         call multifab_copy_c(fluxdiv_prev(n),1,fluxdiv(n),1,nspecies,fluxdiv_prev(n)%ng)
      end do 
   
      !=====================
      ! Euler Predictor step
      !===================== 
     
      ! multiply f(rho,t) with dt (ghosts are copied) 
      do n=1,nlevs
         call multifab_mult_mult_s_c(fluxdiv(n),1,dt,nspecies,0) 
      end do 
      
      ! add this to rho(t) (ghosts added) 
      do n=1,nlevs
         call multifab_plus_plus(rho_prev(n),fluxdiv(n))            
      end do 
      
      ! compute div-of-flux and return results in rho_prev 
      call diffusive_fluxdiv(mla,rho_prev,Dbar,Gama,mass,dx,dt,the_bc_level)
      
      ! fill the ghost cells of rho_prev
      do n=1,nlevs
         call multifab_fill_boundary(rho_prev(n))
        
         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc(rho_prev(n),1,rho_part_bc_comp,nspecies,the_bc_level(n), & 
                              dx(n,:),.false.)
      enddo

      !===================== 
      ! Corrector step
      !=====================
 
      ! add f(rho,t) to that 
      do n=1,nlevs
         call multifab_plus_plus(fluxdiv_prev(n),rho_prev(n))            
      end do 
    
      ! multiply result with dt/2, starting from component 1 with 0 ghost-cell
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
      !write(*,*), 'Using Midpoint explicit method'
      
      ! store the previous values of rho
      do n=1,nlevs
         call multifab_copy_c(rho_prev(n),1,rho(n),1,nspecies,rho_prev(n)%ng)
      end do  
    
      ! multiply fluxdiv with dt/2 
      do n=1,nlevs
         call multifab_mult_mult_s_c(fluxdiv(n),1,dt*0.5d0,nspecies,0) 
      end do 
   
      ! get rho at t+1/2 
      do n=1,nlevs
         call multifab_plus_plus(rho_prev(n),fluxdiv(n))            
      end do

      ! compute fluxdiv[t+1/2,rho_prev(t+1/2)]
      call diffusive_fluxdiv(mla,rho_prev,Dbar,Gama,mass,dx,dt*0.5d0,the_bc_level)
 
      ! fill the ghost cells of rho_prev
      do n=1,nlevs
         call multifab_fill_boundary(rho_prev(n))
        
         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc(rho_prev(n),1,rho_part_bc_comp,nspecies,the_bc_level(n), & 
                              dx(n,:),.false.)
      enddo

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
      !write(*,*), 'Using 3-stage RK explicit method'

      !===========
      ! 1st stage
      !===========
      
      ! store the previous values for 2nd and 3rd stage
      do n=1,nlevs
         call multifab_copy_c(rho_prev(n),1,rho(n),1,nspecies,rho_prev(n)%ng)
         call multifab_copy_c(rho_prev1(n),1,rho(n),1,nspecies,rho_prev1(n)%ng)
      end do 
 
      ! multiply div-of-flux with dt, starting from component 1 with 0 ghost-cell
      do n=1,nlevs
         call multifab_mult_mult_s_c(fluxdiv(n),1,dt,nspecies,0) 
      end do 
    
      ! add this to rho 
      do n=1,nlevs
         call multifab_plus_plus(rho_prev(n),fluxdiv(n))            
      end do 

      ! store this value of rho_prev for 2nd stage
      do n=1,nlevs
         call multifab_copy_c(rho_intm(n),1,rho_prev(n),1,nspecies,rho_intm(n)%ng)
      end do 

      ! compute fluxdiv(t+1,rho_prev(t+1))
      call diffusive_fluxdiv(mla,rho_prev,Dbar,Gama,mass,dx,dt,the_bc_level)
     
      ! fill the ghost cells of rho_prev
      do n=1,nlevs
         call multifab_fill_boundary(rho_prev(n))
        
         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc(rho_prev(n),1,rho_part_bc_comp,nspecies,the_bc_level(n), & 
                              dx(n,:),.false.)
      enddo
 
      !===========
      ! 2nd stage
      !===========

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
         call multifab_copy_c(rho_intm1(n),1,rho_intm(n),1,nspecies,rho_intm1(n)%ng)
      end do 
      
      ! compute fluxdiv(t+1/2,rho_intm(t+1/2))
      call diffusive_fluxdiv(mla,rho_intm,Dbar,Gama,mass,dx,dt*0.5d0,the_bc_level)
      
      ! fill the ghost cells of rho_intm
      do n=1,nlevs
         call multifab_fill_boundary(rho_intm(n))
        
         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc(rho_intm(n),1,rho_part_bc_comp,nspecies,the_bc_level(n), & 
                              dx(n,:),.false.)
      enddo

      !===========
      ! 3rd stage
      !===========

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

    ! fill the ghost cells of rho
    do n=1, nlevs
       call multifab_fill_boundary(rho(n))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(rho(n),1,rho_part_bc_comp,nspecies,the_bc_level(n), & 
                            dx(n,:),.false.)
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
