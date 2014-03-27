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

  public :: print_errors, sum_mass, meanvar_W

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

  subroutine meanvar_W(mla,rho,stdW) 

    type(ml_layout), intent(in)    :: mla
    type(multifab),  intent(in)    :: rho(:)
    real(kind=dp_t), intent(inout) :: stdW(nspecies)

    ! local variables
    type(multifab)                       :: W(mla%nlevel)
    real(kind=dp_t), dimension(nspecies) :: wavg,wsqrtavg,dwsqavg
    integer                              :: i,n,nlevs,n_cell

    nlevs  = size(rho,1)
    n_cell = multifab_volume(rho(1))/nspecies 
  
    ! build cell-centered multifabs for nspecies and ghost cells contained in rho
    do n=1,nlevs
       call multifab_build(W(n),mla%la(n),nspecies,rho(n)%ng)
    end do
   
    ! compute W from rho and rho_tot 
    call convert_rho_to_W(mla,rho,W)
  
    do n=1,nlevs
       do i=1,nspecies

          ! wavg = sum(W)/n_cell
          wavg(i) = multifab_norm_l1_c(W(n),i,1,all=.false.)/dble(n_cell)

          ! wsqrtavg = sqrt{sum(W^2)/n_cell} 
          wsqrtavg(i) = multifab_norm_l2_c(W(n),i,1,all=.false.)/sqrt(dble(n_cell))

          ! standard deviation in W
          dwsqavg(i) = wsqrtavg(i)**2 - wavg(i)**2

          ! sum <delW_i> over time
          stdW(i) = stdW(i) + dwsqavg(i)

       end do
    end do

   ! free the multifab allocated memory
   do n=1,nlevs
      call multifab_destroy(W(n))
   end do

  contains 

     subroutine convert_rho_to_W(mla,rho,W)
   
       type(ml_layout), intent(in   )  :: mla
       type(multifab) , intent(in   )  :: rho(:) 
       type(multifab) , intent(inout)  :: W(:) 

       ! local variables
       integer :: lo(rho(1)%dim), hi(rho(1)%dim)
       integer :: n,i,ng,dm,nlevs
 
       ! pointer for rho(nspecies), rho_tot(1), molarconc(nspecies) 
       real(kind=dp_t), pointer        :: dp(:,:,:,:)   ! for rho    
       real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for W 

       dm    = mla%dim     ! dimensionality
       ng    = rho(1)%ng   ! number of ghost cells 
       nlevs = mla%nlevel  ! number of levels 
 
       ! loop over all boxes 
       do n=1,nlevs
          do i=1,nfabs(rho(n))
             dp => dataptr(rho(n),i)
             dp1 => dataptr(W(n),i)
             lo = lwb(get_box(rho(n),i))
             hi = upb(get_box(rho(n),i))
          
             select case(dm)
             case (2)
                call convert_rho_to_W_2d(dp(:,:,1,:),dp1(:,:,1,:),ng,lo,hi) 
             case (3)
                call convert_rho_to_W_3d(dp(:,:,:,:),dp1(:,:,:,:),ng,lo,hi) 
             end select
          end do
       end do

     end subroutine convert_rho_to_W

     subroutine convert_rho_to_W_2d(rho,W,ng,lo,hi)
 
       integer          :: lo(2), hi(2), ng
       real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)   ! density- last dim for #species
       real(kind=dp_t)  ::   W(lo(1)-ng:,lo(2)-ng:,:)   ! W
        
       ! local variables
       integer          :: i,j
    
       ! for specific box, now start loops over alloted cells    
       do j=lo(2)-ng, hi(2)+ng
          do i=lo(1)-ng, hi(1)+ng
             call convert_rho_to_W_local(rho(i,j,:),W(i,j,:))
          end do
       end do
 
     end subroutine convert_rho_to_W_2d

     subroutine convert_rho_to_W_3d(rho,W,ng,lo,hi)
 
       integer          :: lo(3), hi(3), ng
       real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  ! density- last dim for #species
       real(kind=dp_t)  ::   W(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  ! W 
    
       ! local variables
       integer          :: i,j,k
    
       ! for specific box, now start loops over alloted cells    
       do k=lo(3)-ng, hi(3)+ng
          do j=lo(2)-ng, hi(2)+ng
             do i=lo(1)-ng, hi(1)+ng
                call convert_rho_to_W_local(rho(i,j,k,:),W(i,j,k,:))
             end do
          end do
       end do
 
     end subroutine convert_rho_to_W_3d

     subroutine convert_rho_to_W_local(rho,W)
 
       real(kind=dp_t), intent(in)   :: rho(nspecies)  ! density
       real(kind=dp_t), intent(out)  ::   W(nspecies)    ! W 
    
       ! local variables
       integer          :: n
       real(kind=dp_t)  :: rho_tot        ! total density 

       ! calculate total density inside each cell
       rho_tot=0.d0 
       do n=1, nspecies  
          rho_tot = rho_tot + rho(n)
       end do         
  
       ! calculate mass fraction 
       do n=1, nspecies  
          W(n) = rho(n)/rho_tot
       end do

     end subroutine convert_rho_to_W_local 

  end subroutine meanvar_W

end module analysis_module
