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

  subroutine print_errors(rho,rho_exact,dx,prob_lo,prob_hi,time,the_bc_level)
  
     type(multifab) , intent(in)     :: rho(:)            
     type(multifab) , intent(inout)  :: rho_exact(:)            
     real(kind=dp_t), intent(in   )  :: dx(:,:)           
     real(kind=dp_t), intent(in   )  :: prob_lo(rho(1)%dim)
     real(kind=dp_t), intent(in   )  :: prob_hi(rho(1)%dim)
     real(kind=dp_t), intent(in   )  :: time 
     type(bc_level) , intent(in   )  :: the_bc_level(:)
 
     ! local variables
     integer                         :: i,n,nlevs,n_cell
     real(kind=dp_t)                 :: norm_inf,norm_l1,norm_l2
     
     nlevs = size(rho,1)

     ! calculate rho_exact
     call init_rho(rho_exact,dx,prob_lo,prob_hi,time,the_bc_level) 
    
     ! substract the values 
     do n=1,nlevs
        call multifab_sub_sub_c(rho_exact(n),1,rho(n),1,nspecies,0)
     end do

     n_cell = multifab_volume(rho_exact(1))/nspecies 
     
     ! Linf norm = max(diff_i) 
     ! Donev: There is a bug here: some value in rho_exact is not initialized
     norm_inf = multifab_norm_inf_c(rho_exact(1),1,nspecies,all=.false.)

     ! L1 norm = 1/n_cell*sum(diff_i) 
     norm_l1 = multifab_norm_l1_c(rho_exact(1),1,nspecies,all=.false.)/dble(n_cell)

     ! L2 norm = sqrt{1/n_cell*sum(diff_i^2)} 
     norm_l2 = multifab_norm_l2_c(rho_exact(1),1,nspecies,all=.false.)/sqrt(dble(n_cell))
     
     ! print the norms
     if(.false.) then
       if (parallel_IOProcessor()) then 
          print*, time, norm_inf, norm_l1, norm_l2
       end if
     end if

  end subroutine print_errors

end module analysis_module
