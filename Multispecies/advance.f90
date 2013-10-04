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
    type(multifab) :: rho_star(mla%nlevel),rho_twostar(mla%nlevel),rho_thrstar(mla%nlevel)
    type(multifab) :: fluxdiv(mla%nlevel),fluxdiv_star(mla%nlevel)
    integer        :: n,nlevs
    
    nlevs = mla%nlevel  ! number of levels 
 
    ! build cell-centered multifabs for nspecies and ghost cells in rho,
    ! fluxdiv needs no ghost cells as all computations will be in box interior 
    do n=1,nlevs
       call multifab_build(rho_star(n),    mla%la(n),nspecies,rho(n)%ng)
       call multifab_build(rho_twostar(n), mla%la(n),nspecies,rho(n)%ng)
       call multifab_build(rho_thrstar(n), mla%la(n),nspecies,rho(n)%ng)
       call multifab_build(fluxdiv(n),     mla%la(n),nspecies,0) 
       call multifab_build(fluxdiv_star(n),mla%la(n),nspecies,0)
    enddo
    
    select case(timeinteg_type)
 
      case(1)
      !==================================================================================
      ! Euler time update 
      !==================================================================================
    
      ! compute fluxdiv; fluxdiv will contain results in interior only while rho contains 
      ! ghost values filled in init or end of this code
      call diffusive_fluxdiv(mla,rho,fluxdiv,Dbar,Gama,mass,dx,the_bc_level)
      
      ! compute rho(t+dt)=rho(t)+dt*fluxdiv(t) (only interior) 
      do n=1,nlevs
         ! Donev: You can just call saxpy here
         call multifab_saxpy_3(rho(n),dt,fluxdiv(n))
      end do 
    
      case(2)
      !====================================================================================
      ! Heun's method: Predictor-Corrector explicit method 
      ! rho_star(t+1)=rho(t)+dt*fluxdiv(t,rho) 
      ! rho_new(t+1) =rho(t)+(dt/2)*[fluxdiv(t,rho) + fluxdiv(t+1,rho_star(t+1))] 
      !====================================================================================
      
      ! store old rho in rho_star 
      ! Donev: Instead of this, use saxpy_4 below (just call saxpy generic)
      do n=1,nlevs
         call multifab_copy_c(rho_star(n),1,rho(n),1,nspecies,rho(n)%ng)
      end do 
      
      !=====================
      ! Euler Predictor step
      !===================== 
      
      ! compute fluxdiv 
      call diffusive_fluxdiv(mla,rho,fluxdiv,Dbar,Gama,mass,dx,the_bc_level)
      
      ! compute rho(t+dt)=rho(t)+dt*fluxdiv(t) (only interior) 
      do n=1,nlevs
         call multifab_saxpy_3(rho_star(n),dt,fluxdiv(n))
      end do 
   
      ! update values of the ghost cells of rho_star
      do n=1,nlevs
         call multifab_fill_boundary(rho_star(n))
        
         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc(rho_star(n),1,rho_part_bc_comp,nspecies,the_bc_level(n), & 
                              dx(n,:),.false.)
      enddo

      ! compute new div-of-flux 
      call diffusive_fluxdiv(mla,rho_star,fluxdiv_star,Dbar,Gama,mass,dx,the_bc_level)

      !=========================== 
      ! Trapezoidal Corrector step
      !===========================
 
      ! add f(rho,t) to that 
      do n=1,nlevs
         call multifab_plus_plus(fluxdiv(n),fluxdiv_star(n))            
      end do 
   
      ! compute rho(t+dt)=rho(t)+(dt/2)*fluxdiv(t) (only interior) 
      do n=1,nlevs
         call multifab_saxpy_3(rho(n),0.5d0*dt,fluxdiv(n))
      end do 
      
      case(3)
      !=====================================================================================
      ! Midpoint method (2-stage) 
      ! rho_star(t+1/2)=rho(t)+(dt/2)*fluxdiv(t,y)                
      ! rho_new(t+1)=rho(t)+ dt*fluxdiv[t+1/2,rho_star(t+1/2)]  
      !=====================================================================================
 
      ! store old rho in rho_star 
      do n=1,nlevs
         call multifab_copy_c(rho_star(n),1,rho(n),1,nspecies,rho(n)%ng)
      end do 

      ! compute fluxdiv(t) from rho(t); fluxdiv will contain results in interior 
      ! only while rho contains ghost values filled in init or end of this code
      call diffusive_fluxdiv(mla,rho,fluxdiv,Dbar,Gama,mass,dx,the_bc_level)
      
      ! compute rho_star(t+dt)=rho(t)+(dt/2)*fluxdiv(t) (only interior) 
      do n=1,nlevs
         call multifab_saxpy_3(rho_star(n),0.5d0*dt,fluxdiv(n))
      end do 

      ! update values of the ghost cells of rho_star
      do n=1,nlevs
         call multifab_fill_boundary(rho_star(n))
        
         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc(rho_star(n),1,rho_part_bc_comp,nspecies,the_bc_level(n), & 
                              dx(n,:),.false.)
      enddo

      ! compute new div-of-flux 
      call diffusive_fluxdiv(mla,rho_star,fluxdiv_star,Dbar,Gama,mass,dx,the_bc_level)

      ! compute rho(t+dt)=rho(t)+dt*fluxdiv(t+dt/2,rho_star(t+1/2)) (only interior) 
      do n=1,nlevs
         call multifab_saxpy_3(rho(n),dt,fluxdiv_star(n))
      end do 
      
      case(4)
      !=======================================================================================
      ! 3rd order Runge-Kutta method 
      ! rho_star(t+dt)     =     rho(t)+                            dt*fluxdiv(t)                      
      ! rho_twostar(t+dt/2)= 3/4*rho(t)+ 1/4*[rho_star(t+dt)      + dt*fluxdiv(t+dt,  rho_star(t+dt))] 
      ! rho_new(t+dt)      = 1/3*rho(t)+ 2/3*[rho_twostar(t+dt/2) + dt*fluxdiv(t+dt/2,rho_twostar(t+dt/2))] 
      !=======================================================================================

      ! store the previous rho for 2nd and 3rd stage
      do n=1,nlevs
         call multifab_copy_c(rho_star(n),   1,rho(n),1,nspecies,rho(n)%ng)
         call multifab_copy_c(rho_twostar(n),1,rho(n),1,nspecies,rho(n)%ng)
         call multifab_copy_c(rho_thrstar(n),1,rho(n),1,nspecies,rho(n)%ng)
      end do 
      
      !===========
      ! 1st stage
      !===========

      ! compute fluxdiv(t) from rho(t); fluxdiv will contain results in interior 
      ! only while rho contains ghost values filled in init or end of this code
      call diffusive_fluxdiv(mla,rho,fluxdiv,Dbar,Gama,mass,dx,the_bc_level)
      
      ! compute rho_star(t+dt)=rho(t)+dt*fluxdiv(t) (only interior) 
      do n=1,nlevs
         call multifab_saxpy_3(rho_star(n),dt,fluxdiv(n))
      end do 
   
      ! update values of the ghost cells of rho_star
      do n=1,nlevs
         call multifab_fill_boundary(rho_star(n))
        
         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc(rho_star(n),1,rho_part_bc_comp,nspecies,the_bc_level(n), & 
                              dx(n,:),.false.)
      enddo

      ! compute fluxdiv_star=fluxdiv(t+dt,rho_star(t+dt))
      call diffusive_fluxdiv(mla,rho_star,fluxdiv_star,Dbar,Gama,mass,dx,the_bc_level)

      !===========
      ! 2nd stage
      !===========

      ! compute rho_star(t+dt)+dt*fluxdiv(t+dt) (only interior) 
      do n=1,nlevs
         call multifab_saxpy_3(rho_star(n),dt,fluxdiv_star(n))
      end do 

      ! multiply old rho with 3/4 (stored in rho_twostar)
      do n=1,nlevs
         call multifab_mult_mult_s_c(rho_twostar(n),1,0.75d0,nspecies,0) 
      end do 

      ! rho_twostar(t+1/2)= 3/4*rho(t)+ 1/4*[rho_star(t+dt) + dt*fluxdiv(t+dt,rho_star(t+dt))] 
      do n=1,nlevs
         call multifab_saxpy_3(rho_twostar(n),0.25d0,rho_star(n))
      end do 

      ! update values of the ghost cells of rho_star
      do n=1,nlevs
         call multifab_fill_boundary(rho_twostar(n))
        
         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc(rho_twostar(n),1,rho_part_bc_comp,nspecies,the_bc_level(n), & 
                              dx(n,:),.false.)
      enddo

      ! compute fluxdiv(t+dt/2,rho_twostar(t+dt/2)) and store in fluxdiv_star  
      call diffusive_fluxdiv(mla,rho_twostar,fluxdiv_star,Dbar,Gama,mass,dx,the_bc_level)

      !===========
      ! 3rd stage
      !===========

      ! compute rho_star(t+dt)=rho(t)+dt*fluxdiv(t) (only interior) 
      do n=1,nlevs
         call multifab_saxpy_3(rho_twostar(n),dt,fluxdiv_star(n))
      end do 

      ! multiply old rho with 1/3
      do n=1,nlevs
         call multifab_mult_mult_s_c(rho(n),1,(1.0d0/3.0d0),nspecies,0) 
      end do 

      ! rho_twostar(t+1/2)= 3/4*rho(t)+ 1/4*rho_star(t+1)
      do n=1,nlevs
         call multifab_saxpy_3(rho(n),(2.0d0/3.0d0),rho_twostar(n))
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
       call multifab_destroy(rho_star(n))
       call multifab_destroy(rho_twostar(n))
       call multifab_destroy(rho_thrstar(n))
       call multifab_destroy(fluxdiv(n))
       call multifab_destroy(fluxdiv_star(n))
    end do

  end subroutine advance

end module advance_module 
