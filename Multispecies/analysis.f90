module analysis_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use define_bc_module
  use bc_module
  use init_module
  use probin_common_module
  use probin_multispecies_module
 
  implicit none

  private

  public :: print_errors

  contains

  subroutine print_errors(rho,rho_exact,Dbar,Gama,mass,dx,prob_lo,prob_hi,time,the_bc_level)
  
     type(multifab) , intent(inout)  :: rho(:)            
     type(multifab) , intent(inout)  :: rho_exact(:) 
     real(kind=dp_t), intent(inout)  :: Dbar(:,:)           
     real(kind=dp_t), intent(inout)  :: Gama(:,:)           
     real(kind=dp_t), intent(inout)  :: mass(:)           
     real(kind=dp_t), intent(in   )  :: dx(:,:)           
     real(kind=dp_t), intent(in   )  :: prob_lo(rho(1)%dim)
     real(kind=dp_t), intent(in   )  :: prob_hi(rho(1)%dim)
     real(kind=dp_t), intent(in   )  :: time 
     type(bc_level) , intent(in   )  :: the_bc_level(:)
 
     ! local variables
     integer                         :: i,n,nlevs,n_cell
     real(kind=dp_t)                 :: norm

     nlevs = size(rho,1)

     ! calculate rho_exact defined in init.f90
     call init_rho(rho,rho_exact,Dbar,Gama,mass,dx,prob_lo,prob_hi,time,the_bc_level)
    
     ! substract the values 
     do n=1,nlevs
        call multifab_sub_sub_c(rho_exact(n),1,rho(n),1,nspecies,0)
     end do

     n_cell = multifab_volume(rho_exact(1))/nspecies 
     norm = multifab_norm_inf_c(rho_exact(1),1,nspecies,all=.false.)
    
     ! Linf norm = max(x_i) 
     !if (parallel_IOProcessor()) print*,"Linf-norm rho  =",norm

     ! L1 norm = 1/n_cell*sum(x_i) 
     !norm = multifab_norm_l1_c(rho_exact(1),1,nspecies,all=.false.)/dble(n_cell)
     !if (parallel_IOProcessor()) print*,"L1-norm rho  =",norm

     ! L2 norm = sqrt{1/n_cell*sum(x_i^2)} 
     !norm = multifab_norm_l2_c(rho_exact(1),1,nspecies,all=.false.)/sqrt(dble(n_cell))
     !if (parallel_IOProcessor()) print*,"L2-norm rho  =",norm
     !if (parallel_IOProcessor()) print*,""
     
     if (parallel_IOProcessor()) then 
       if(time .gt. 1) then 
          print*, time, norm
       endif
     endif

  end subroutine print_errors

end module analysis_module
