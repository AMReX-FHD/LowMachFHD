module update_rho_module

  use multifab_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use ml_layout_module
  use probin_multispecies_module

  implicit none

  private

  public :: update_rho

contains

  subroutine update_rho(mla,rho,fluxdiv,dt,the_bc_level,dx,rho_part_bc_comp)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rho(:)
    type(multifab) , intent(inout) :: fluxdiv(mla%nlevel)
    real(kind=dp_t), intent(in   ) :: dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    integer,         intent(in   ) :: rho_part_bc_comp 

    ! local variables
    integer n, nlevs

    nlevs = mla%nlevel  ! number of levels 
   
    select case(timeinteg_type)
 
    case(1)
    !%%%%%% Euler explicit time update %%%%%%%%! 
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
    !%%%%%% Option: Predictor-Corrector explicit method %%%%%%%%%%%%%%%! 

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
