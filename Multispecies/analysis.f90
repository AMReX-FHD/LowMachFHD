module analysis_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use define_bc_module
  use ml_layout_module
  use init_module
  use probin_common_module
  use probin_multispecies_module
 
  implicit none

  private

  public :: print_errors, sum_mass, meanvar_w

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
     do n=1,nlevs
        do i=1, nspecies 
         
           ! Linf norm = max(diff_i)
           norm_inf(i) = multifab_norm_inf_c(rho_exact(n),i,1,all=.false.)

           ! L1 norm = 1/n_cell*sum(diff_i) 
           norm_l1(i) = multifab_norm_l1_c(rho_exact(n),i,1,all=.false.)/dble(n_cell)

           ! L2 norm = sqrt{1/n_cell*sum(diff_i^2)} 
           norm_l2(i) = multifab_norm_l2_c(rho_exact(n),i,1,all=.false.)/sqrt(dble(n_cell))
     
        end do
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
    integer,        intent(in   ) :: step

    ! local variables
    real(kind=dp_t) :: mass(nspecies)
    integer         :: i,n,nlevs

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

  subroutine meanvar_w(mla,rho,rho_tot,the_bc_level) 

    type(ml_layout), intent(in) :: mla
    type(multifab),  intent(in) :: rho(:)
    type(multifab),  intent(in) :: rho_tot(:)
    type(bc_level),  intent(in) :: the_bc_level(:)

    ! local variables
    type(multifab)                       :: w(mla%nlevel)
    real(kind=dp_t), dimension(nspecies) :: wavg,wsqrtavg,dwsqavg
    integer                              :: i,n,nlevs,n_cell

    nlevs  = size(rho,1)
    n_cell = multifab_volume(rho(1))/nspecies 
  
    ! build cell-centered multifabs for nspecies and ghost cells contained in rho.
    do n=1,nlevs
       call multifab_build(w(n),mla%la(n),nspecies,rho(n)%ng)
    end do

    ! clear up memory 
    do n=1,nlevs
       call setval(w(n), 0.d0, all=.true.)
    end do

    ! store rho in w 
    do n=1,nlevs
       call saxpy(w(n),1.0d0,rho(n))
    end do

    ! compute w_i=rho_i/rho_tot
    do n=1,nlevs
       call multifab_div_div_c(w(n),1,rho_tot(n),1,nspecies,0)
    end do
 
    do n=1,nlevs
       do i=1,nspecies

          ! wavg = sum(w)/n_cell
          wavg(i) = multifab_norm_l1_c(w(n),i,1,all=.false.)/dble(n_cell)

          ! wsqrtavg = sqrt{sum(w^2)/n_cell} 
          wsqrtavg(i) = multifab_norm_l2_c(w(n),i,1,all=.false.)/sqrt(dble(n_cell))

          ! standard deviation in w
          dwsqavg(i) = wsqrtavg(i)**2 - wavg(i)**2

          if (parallel_IOProcessor()) then
                print*, ' std of w for i=',i, dwsqavg(i)
          end if
    
       end do
    end do

   ! free the multifab allocated memory
   do n=1,nlevs
      call multifab_destroy(w(n))
   end do

  end subroutine meanvar_w

end module analysis_module
