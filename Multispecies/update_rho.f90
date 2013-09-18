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
    type(multifab)                 :: rho_predictor(mla%nlevel)
    type(multifab)                 :: fluxdiv_prev(mla%nlevel)
    integer                        :: n, nlevs
    
    nlevs = mla%nlevel  ! number of levels 
   
    select case(timeinteg_type)
 
    case(1)
    !%%%%%%%%%%%%% Euler explicit time update %%%%%%%%! 
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
    !%%%%%% Predictor-Corrector explicit method %%%%%%%! 
    !%%%%%% rho_predictor(t+1)=rho(t)+dt*fluxdiv(t) %%%%%%%! 
    !%%%%%% rho_corrector(t+1)=rho(t)+(dt/2)*[rho(t)+rho_predictor(t+1)] %%%%%%%! 
    do n=1,nlevs
       ! store the previous values of rho
       rho_predictor(n) = rho(n)       
       fluxdiv_prev(n)  = fluxdiv(n) 
    end do 
   
    ! Euler Predictor step 
    do n=1,nlevs
       ! multiply div-of-flux with dt, starting from component 1 with 0 ghost-cell
       call multifab_mult_mult_s_c(fluxdiv(n),1,dt,nspecies,0) 
    end do 

    do n=1,nlevs
       ! add this to rho to advance in time
       call multifab_plus_plus(rho_predictor(n),fluxdiv(n))            
    end do 

    ! calculate new div-of-flux with Predictor results calling advance routine
    ! without the update_rho call which is saved as fluxdiv_pred
    call fluxdiv_pred(mla,rho_predictor,Dbar,Gama,mass,dx,dt,the_bc_level,&
                 rho_part_bc_comp,mol_frac_bc_comp,diff_coeff_bc_comp)

    ! Corrector step
    ! New fluxdiv is stored in returned quantity rho_predictor 
    do n=1,nlevs
       call multifab_plus_plus(fluxdiv_prev(n),rho_predictor(n))            
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

    end select  

    do n=1, nlevs
       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(rho(n))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(rho(n),1,rho_part_bc_comp,nspecies,the_bc_level(n), & 
                            dx(n,:),.false.)
    end do   

  end subroutine update_rho

end module update_rho_module 
