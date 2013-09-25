module update_rho_module

  use multifab_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use ml_layout_module
  use fluxdiv_pred_module
  use probin_multispecies_module

  implicit none

  private

  public :: update_rho

contains

  subroutine update_rho(mla,rho,fluxdiv,Dbar,Gama,mass,dx,dt,the_bc_level,& 
                        rho_part_bc_comp,mol_frac_bc_comp,diff_coeff_bc_comp)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rho(:)
    type(multifab) , intent(inout) :: fluxdiv(mla%nlevel)
    real(kind=dp_t), intent(in   ) :: Dbar(:,:)
    real(kind=dp_t), intent(in   ) :: Gama(:,:)
    real(kind=dp_t), intent(in   ) :: mass(:) 
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(in   ) :: dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    integer,         intent(in   ) :: rho_part_bc_comp,mol_frac_bc_comp
    integer,         intent(in   ) :: diff_coeff_bc_comp 

    ! local variables
    type(multifab) :: rho_prev(mla%nlevel), rho_prev1(mla%nlevel)
    type(multifab) :: rho_intm(mla%nlevel), rho_intm1(mla%nlevel)
    type(multifab) :: fluxdiv_prev(mla%nlevel)
    integer        :: n, nlevs
    
    nlevs = mla%nlevel  ! number of levels 
    !do n=1,nlevs
    !   call setval(rho_prev(n),0.d0,all=.true.)
    !   call setval(rho_prev1(n),0.d0,all=.true.)
    !   call setval(rho_intm(n),0.d0,all=.true.)
    !   call setval(rho_intm1(n),0.d0,all=.true.)
    !   call setval(fluxdiv_prev(n),0.d0,all=.true.)
    !end do
 
    select case(timeinteg_type)
 
    case(1)
    !%%%%%% Euler time update 
    do n=1,nlevs
       ! multiply div-of-flux with dt, starting from component 1 with 0 ghost-cell
       call multifab_mult_mult_s_c(fluxdiv(n),1,dt,nspecies,0) 
    end do 
    
    do n=1,nlevs
       ! add this to rho to advance in time
       call multifab_plus_plus(rho(n),fluxdiv(n))            
    end do 
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%! 

    case(2)
    !%%%%%% Predictor-Corrector explicit method 
    !%%%%%% rho_predictor(t+1)=rho(t)+dt*fluxdiv(t,rho) %%%%%%%! 
    !%%%%%% rho_corrector(t+1)=rho(t)+(dt/2)*[fluxdiv(t,rho) + fluxdiv(t+1,rho_predictor(t+1))] %%%%%%%! 
    
    do n=1,nlevs
       ! store the previous values for nspecies and 0 ghost cells
       call multifab_copy_c(rho_prev(n),1,rho(n),1,nspecies,0)
       call multifab_copy_c(fluxdiv_prev(n),1,fluxdiv(n),1,nspecies,0)
    end do 
   
    ! Euler Predictor step 
    do n=1,nlevs
       call multifab_mult_mult_s_c(fluxdiv(n),1,dt,nspecies,0) 
    end do 

    do n=1,nlevs
       call multifab_plus_plus(rho_prev(n),fluxdiv(n))            
    end do 

    ! compute div-of-flux with Predictor results (almost as calling advance 
    ! routine without the update_rho call which is saved as fluxdiv_pred)
    call fluxdiv_pred(mla,rho_prev,Dbar,Gama,mass,dx,dt,the_bc_level,&
                 rho_part_bc_comp,mol_frac_bc_comp,diff_coeff_bc_comp)

    ! Corrector step
    ! New fluxdiv is stored in returned quantity rho_prev 
    do n=1,nlevs
       call multifab_plus_plus(fluxdiv_prev(n),rho_prev(n))            
    end do 
    
    do n=1,nlevs
       ! multiply div-of-flux with dt, starting from component 1 with 0 ghost-cell
       call multifab_mult_mult_s_c(fluxdiv_prev(n),1,dt*0.5d0,nspecies,0) 
    end do 
    
    do n=1,nlevs
       ! add this to rho to advance in time
       call multifab_plus_plus(rho(n),fluxdiv_prev(n))            
    end do 
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%! 
    
    case(3)
    !%%%%%% Midpoint method (2-stage) 
    !%%%%%% rho_prev(t+1/2)=rho(t)+(dt/2)*fluxdiv(t,y)                %%%! 
    !%%%%%% rho_new(t+1)=rho(t)+ dt*fluxdiv[t+1/2,rho_prev(t+1/2)] %%%! 
    
    do n=1,nlevs
       ! store the previous values of rho
       call multifab_copy_c(rho_prev(n),1,rho(n),1,nspecies,0)
    end do  
    
    do n=1,nlevs
       call multifab_mult_mult_s_c(fluxdiv(n),1,dt*0.5d0,nspecies,0) 
    end do 
    
    do n=1,nlevs
       call multifab_plus_plus(rho_prev(n),fluxdiv(n))            
    end do

    call fluxdiv_pred(mla,rho_prev,Dbar,Gama,mass,dx,dt*0.5d0,the_bc_level,&
                 rho_part_bc_comp,mol_frac_bc_comp,diff_coeff_bc_comp)

    do n=1,nlevs
       ! multiply div-of-flux with dt, starting from component 1 with 0 ghost-cell
       call multifab_mult_mult_s_c(rho_prev(n),1,dt,nspecies,0) 
    end do 
    
    do n=1,nlevs
       ! add this to rho to advance in time
       call multifab_plus_plus(rho(n),rho_prev(n))            
    end do  
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%! 

    case(4)
    !%%%%%% 3rd order Runge-Kutta method 
    !%%%%%% rho_prev(t+1)  =     rho(t)+                        dt*fluxdiv(t)                     %%%! 
    !%%%%%% rho_intm(t+1/2)= 3/4*rho(t)+ 1/4*[rho_prev(t+1)   + dt*fluxdiv(t+1,  rho_prev(t+1))  ]%%%! 
    !%%%%%% rho_new(t+1)   = 1/3*rho(t)+ 2/3*[rho_intm(t+1/2) + dt*fluxdiv(t+1/2,rho_intm(t+1/2))]%%%! 

    ! 1st stage
    do n=1,nlevs
       ! store the previous values for 2nd and 3rd stage
       call multifab_copy_c(rho_prev(n),1,rho(n),1,nspecies,0)
       call multifab_copy_c(rho_prev1(n),1,rho(n),1,nspecies,0)
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
       call multifab_copy_c(rho_intm(n),1,rho_prev(n),1,nspecies,0)
    end do 

    ! 2nd stage
    call fluxdiv_pred(mla,rho_prev,Dbar,Gama,mass,dx,dt,the_bc_level,&
                      rho_part_bc_comp,mol_frac_bc_comp,diff_coeff_bc_comp)

    do n=1,nlevs
       ! multiply fluxdiv with dt/4
       call multifab_mult_mult_s_c(rho_prev(n),1,0.25d0*dt,nspecies,0) 
    end do 

    do n=1,nlevs
       ! multiply previous rho with 1/4
       call multifab_mult_mult_s_c(rho_intm(n),1,0.25d0,nspecies,0) 
    end do 
    
    do n=1,nlevs
       ! add this two 
       call multifab_plus_plus(rho_intm(n),rho_prev(n))            
    end do 

    do n=1,nlevs
       ! multiply first saved rho with 3/4
       call multifab_mult_mult_s_c(rho_prev1(n),1,0.75d0,nspecies,0) 
    end do 
    
    do n=1,nlevs
       ! add them all
       call multifab_plus_plus(rho_intm(n),rho_prev1(n))            
    end do 

    do n=1,nlevs
       ! store this intermediate value of rho for 3rd stage
       call multifab_copy_c(rho_intm1(n),1,rho_intm(n),1,nspecies,0)
    end do 
    
    ! 3rd stage: Amit: After this step code crashing
    call fluxdiv_pred(mla,rho_intm,Dbar,Gama,mass,dx,dt*0.5d0,the_bc_level,&
                      rho_part_bc_comp,mol_frac_bc_comp,diff_coeff_bc_comp)

    do n=1,nlevs
       ! multiply div-of-flux with 2*dt/3
       call multifab_mult_mult_s_c(rho_intm(n),1,(2.0d0/3.0d0)*dt,nspecies,0) 
    end do 
    
    do n=1,nlevs
       ! multiply stored intermediate rho with 2/3 
       call multifab_mult_mult_s_c(rho_intm1(n),1,(2.0d0/3.0d0),nspecies,0) 
    end do 

     do n=1,nlevs
       ! add this to rho to advance in time
       call multifab_plus_plus(rho_intm(n),rho_intm1(n))            
    end do 
     
    do n=1,nlevs
       ! multiply oldest rho with 1/3 
       call multifab_mult_mult_s_c(rho(n),1,(1.0d0/3.0d0),nspecies,0) 
    end do 
    
    do n=1,nlevs
       ! add this to rho to advance in time
       call multifab_plus_plus(rho(n),rho_intm(n))            
    end do 
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%! 

    end select  

    do n=1, nlevs
       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(rho(n))
       !call multifab_fill_boundary(fluxdiv(n))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(rho(n),1,rho_part_bc_comp,nspecies,the_bc_level(n), & 
                            dx(n,:),.false.)
    end do   

    ! destroy the multifab to prevent leakage in memory
    do n=1,nlevs
       call multifab_destroy(rho_prev(n))
       call multifab_destroy(rho_prev1(n))
       call multifab_destroy(rho_intm(n))
       call multifab_destroy(rho_intm1(n))
       call multifab_destroy(fluxdiv_prev(n))
    end do

  end subroutine update_rho

end module update_rho_module 
