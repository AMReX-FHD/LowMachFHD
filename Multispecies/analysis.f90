module analysis_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use define_bc_module
  use init_module
  use probin_common_module
  use probin_multispecies_module
 
  implicit none

  private

  public :: print_errors, sum_mass

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
     real(kind=dp_t), dimension(nspecies) :: norm_inf,norm_l1,norm_l2
     real(kind=dp_t)                      :: norm_inf_tot,norm_l1_tot,norm_l2_tot
     integer                              :: i,n,nlevs,n_cell
     
     nlevs = size(rho,1)

     ! calculate rho_exact
     call init_rho(rho_exact,dx,prob_lo,prob_hi,time,the_bc_level) 
    
     ! substract the values 
     do n=1,nlevs
        call multifab_sub_sub_c(rho_exact(n),1,rho(n),1,nspecies,0)
     end do

     n_cell = multifab_volume(rho_exact(1))/nspecies 
     
     ! calculate norms for each species
     do i=1, nspecies 
         
        ! Linf norm = max(diff_i)
        norm_inf(i) = multifab_norm_inf_c(rho_exact(1),i,1,all=.false.)

        ! L1 norm = 1/n_cell*sum(diff_i) 
        norm_l1(i) = multifab_norm_l1_c(rho_exact(1),i,1,all=.false.)/dble(n_cell)

        ! L2 norm = sqrt{1/n_cell*sum(diff_i^2)} 
        norm_l2(i) = multifab_norm_l2_c(rho_exact(1),i,1,all=.false.)/sqrt(dble(n_cell))
     
     end do
     
     ! calculate total norms 
     norm_inf_tot = multifab_norm_inf_c(rho_exact(1),1,nspecies,all=.false.)
     norm_l1_tot = multifab_norm_l1_c(rho_exact(1),1,nspecies,all=.false.)/dble(n_cell)
     norm_l2_tot = multifab_norm_l2_c(rho_exact(1),1,nspecies,all=.false.)/sqrt(dble(n_cell))
     
     ! print the norms
     if(.false.) then
        if (parallel_IOProcessor()) then 
            if(time.gt.2.99d0) then
            !if(time.gt.0.49d0) then
               write(*,*), time, norm_inf(1:nspecies), norm_l1(1:nspecies),norm_l2(1:nspecies),&
                       norm_inf_tot,norm_l1_tot,norm_l2_tot
            end if
        end if
     end if
     
     ! for checking analytic solution with Visit
     call init_rho(rho_exact,dx,prob_lo,prob_hi,time,the_bc_level) 

  end subroutine print_errors

  subroutine sum_mass(rho, step)

    type(multifab), intent(in   ) :: rho(:)
    integer, intent(in) :: step

    ! local
    real(kind=dp_t) :: mass(nspecies)
    integer :: i, n, nlevs

    nlevs = size(rho,1)
    
    do n=1,nlevs
       do i=1,nspecies
          mass(i) = multifab_norm_l1_c(rho(n),i,1,all=.false.)
          if (parallel_IOProcessor()) then
             print*, step, ' sum of rho_i for i=',i,mass(i)
          end if
       end do
       if (parallel_IOProcessor()) then
          print*, step, ' sum of rho=',sum(mass(1:nspecies))
       end if
    end do

  end subroutine sum_mass

end module analysis_module
