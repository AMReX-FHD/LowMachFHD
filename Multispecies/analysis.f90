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

  public :: print_errors, sum_mass, compute_cov, meanvar_W

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

  subroutine compute_cov(mla,rho,covW) 

    type(ml_layout), intent(in)     :: mla
    type(multifab),  intent(in)     :: rho(:)
    real(kind=dp_t), intent(inout)  :: covW(nspecies,nspecies)

    ! local variables
    integer :: lo(rho(1)%dim), hi(rho(1)%dim)
    integer :: n,i,j,ng,dm,nlevs,n_cell
    
    ! pointer for rho 
    real(kind=dp_t), pointer  :: dp(:,:,:,:)  

    dm     = mla%dim     ! dimensionality
    ng     = rho(1)%ng   ! number of ghost cells 
    nlevs  = mla%nlevel 
    n_cell = multifab_volume(rho(1))/nspecies 
  
    ! loop over all boxes 
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp => dataptr(rho(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          
          select case(dm)
          case (2)
             call compute_cov_2d(dp(:,:,1,:),n_cell,covW(:,:),ng,lo,hi) 
          case (3)
             call compute_cov_3d(dp(:,:,:,:),n_cell,covW(:,:),ng,lo,hi) 
          end select
       end do
    end do

    end subroutine compute_cov
     
    subroutine compute_cov_2d(rho,n_cell,covW,ng,lo,hi)
 
       integer         :: lo(2), hi(2), ng, n_cell
       real(kind=dp_t) :: rho(lo(1)-ng:,lo(2)-ng:,:) 
       real(kind=dp_t) :: covW(nspecies,nspecies)  
   
       ! local variables
       real(kind=dp_t), dimension(nspecies)          :: cellW
       real(kind=dp_t), dimension(nspecies,nspecies) :: cellWij
       integer                                       :: i,j
    
       cellW   = 0.d0
       cellWij = 0.d0
       
       ! for specific box, now start loops over alloted cells    
       do j=lo(2)-ng, hi(2)+ng
          do i=lo(1)-ng, hi(1)+ng
             call compute_cov_local(rho(i,j,:),cellW(:),cellWij(:,:))
          end do
       end do
    
       do i=1, nspecies
          !print*, cellW(i)/dble(n_cell)
       end do
       ! average over n_cell and calculate covW
       do i=1,nspecies
          do j=1, nspecies     
             covW(i,j) = covW(i,j) + cellWij(i,j)/dble(n_cell) - cellW(i)*cellW(j)/dble(n_cell**2)
          end do
       end do

     end subroutine compute_cov_2d

     subroutine compute_cov_3d(rho,n_cell,covW,ng,lo,hi)
 
       integer          :: lo(3), hi(3), ng, n_cell
       real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:) 
       real(kind=dp_t)  :: covW(nspecies,nspecies)  
    
       ! local variables
       real(kind=dp_t), dimension(nspecies)          :: cellW
       real(kind=dp_t), dimension(nspecies,nspecies) :: cellWij
       integer                                       :: i,j,k
   
       cellW   = 0.d0
       cellWij = 0.d0
 
       ! for specific box, now start loops over alloted cells    
       do k=lo(3)-ng, hi(3)+ng
          do j=lo(2)-ng, hi(2)+ng
             do i=lo(1)-ng, hi(1)+ng
                call compute_cov_local(rho(i,j,k,:),cellW(:),cellWij(:,:))
             end do
          end do
       end do

       ! average over n_cell and calculate covW
       do i=1,nspecies
          do j=1, nspecies     
             covW(i,j) = covW(i,j) + cellWij(i,j)/dble(n_cell) - cellW(i)*cellW(j)/dble(n_cell**2)
          end do
       end do

     end subroutine compute_cov_3d

     subroutine compute_cov_local(rho,cellW,cellWij)
 
       real(kind=dp_t)  :: rho(nspecies)  ! density
       real(kind=dp_t)  :: cellW(nspecies)    
       real(kind=dp_t)  :: cellWij(nspecies,nspecies) 
    
       ! local variables
       real(kind=dp_t), dimension(nspecies)           :: W
       real(kind=dp_t), dimension(nspecies,nspecies)  :: Wij
       real(kind=dp_t)                                :: rho_tot             
       integer                                        :: i,j

       ! calculate total density inside each cell
       rho_tot=0.d0 
       do i=1, nspecies  
          rho_tot = rho_tot + rho(i)
       end do         
  
       ! calculate mass fraction and sum over cell for each species
       do i=1, nspecies  
          W(i)     = rho(i)/rho_tot
          cellW(i) = cellW(i) + W(i) ! compute this for average
       end do

       ! calculate Wij=wi*wj matrix and <Wij>
       do i=1,nspecies
          do j=1, i-1
                 Wij(i,j) = W(i)*W(j)
                 Wij(j,i) = Wij(i,j) 
             cellWij(i,j) = cellWij(i,j) + Wij(i,j) ! compute for average 
             cellWij(j,i) = cellWij(j,i) + Wij(j,i) 
          end do
              Wij(i,i) = W(i)*W(i)
          cellWij(i,i) = cellWij(i,i) + Wij(i,i) 
       end do 

     end subroutine compute_cov_local 
 
    subroutine meanvar_W(mla,rho,covW) 

    type(ml_layout), intent(in)    :: mla
    type(multifab),  intent(in)    :: rho(:)
    real(kind=dp_t), intent(inout) :: covW(nspecies,nspecies)

    ! local variables
    type(multifab)                       :: W(mla%nlevel)
    type(multifab)                       :: Wij(mla%nlevel)
    real(kind=dp_t), dimension(nspecies) :: wavg,wsqrtavg
    real(kind=dp_t), dimension(nspecies,nspecies) :: wiwjavg
    integer                              :: i,j,n,nlevs,n_cell

    nlevs  = size(rho,1)
    n_cell = multifab_volume(rho(1))/nspecies 
  
    ! build cell-centered multifabs for nspecies and ghost cells contained in rho
    do n=1,nlevs
       call multifab_build(W(n),mla%la(n),nspecies,rho(n)%ng)
       call multifab_build(Wij(n),mla%la(n),nspecies**2,rho(n)%ng)
    end do
   
    ! compute W and Wij from rho and rho_tot 
    call convert_rho_to_W(mla,rho,W,Wij)
  
    do n=1,nlevs
       do i=1,nspecies
          ! wavg = sum(W)/n_cell
          wavg(i) = multifab_norm_l1_c(W(n),i,1,all=.false.)/dble(n_cell)
          !print*, wavg(i)
       end do
       do i=1,nspecies
          do j=1, nspecies     
  
          ! <w_i*w_j>
          wiwjavg(i,j) = multifab_norm_l1_c(Wij(n),i*j,1,all=.false.)/dble(n_cell)
 
          ! sum covariance (W_i W_j) over time
          covW(i,j) = covW(i,j) + wiwjavg(i,j) - wavg(i)*wavg(j)

          end do
       end do
   end do

   ! free the multifab allocated memory
   do n=1,nlevs
      call multifab_destroy(W(n))
      call multifab_destroy(Wij(n))
   end do

  contains 

     subroutine convert_rho_to_W(mla,rho,W,Wij)
   
       type(ml_layout), intent(in   )  :: mla
       type(multifab) , intent(in   )  :: rho(:) 
       type(multifab) , intent(inout)  :: W(:) 
       type(multifab) , intent(inout)  :: Wij(:) 

       ! local variables
       integer :: lo(rho(1)%dim), hi(rho(1)%dim)
       integer :: n,i,ng,dm,nlevs
 
       ! pointer for rho(nspecies), rho_tot(1), molarconc(nspecies) 
       real(kind=dp_t), pointer        :: dp(:,:,:,:)   ! for rho    
       real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for W 
       real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for Wij 

       dm    = mla%dim     ! dimensionality
       ng    = rho(1)%ng   ! number of ghost cells 
       nlevs = mla%nlevel  ! number of levels 
 
       ! loop over all boxes 
       do n=1,nlevs
          do i=1,nfabs(rho(n))
             dp => dataptr(rho(n),i)
             dp1 => dataptr(W(n),i)
             dp2 => dataptr(Wij(n),i)
             lo = lwb(get_box(rho(n),i))
             hi = upb(get_box(rho(n),i))
          
             select case(dm)
             case (2)
                call convert_rho_to_W_2d(dp(:,:,1,:),dp1(:,:,1,:),dp2(:,:,1,:),ng,lo,hi) 
             case (3)
                call convert_rho_to_W_3d(dp(:,:,:,:),dp1(:,:,:,:),dp2(:,:,:,:),ng,lo,hi) 
             end select
          end do
       end do

     end subroutine convert_rho_to_W

     subroutine convert_rho_to_W_2d(rho,W,Wij,ng,lo,hi)
 
       integer          :: lo(2), hi(2), ng
       real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)   ! density- last dim for #species
       real(kind=dp_t)  ::   W(lo(1)-ng:,lo(2)-ng:,:)   ! W
       real(kind=dp_t)  :: Wij(lo(1)-ng:,lo(2)-ng:,:)   ! Wij
        
       ! local variables
       integer          :: i,j
    
       ! for specific box, now start loops over alloted cells    
       do j=lo(2)-ng, hi(2)+ng
          do i=lo(1)-ng, hi(1)+ng
             call convert_rho_to_W_local(rho(i,j,:),W(i,j,:),Wij(i,j,:))
          end do
       end do
 
     end subroutine convert_rho_to_W_2d

     subroutine convert_rho_to_W_3d(rho,W,Wij,ng,lo,hi)
 
       integer          :: lo(3), hi(3), ng
       real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  ! density- last dim for #species
       real(kind=dp_t)  ::   W(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  ! W 
       real(kind=dp_t)  :: Wij(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  ! Wij 
    
       ! local variables
       integer          :: i,j,k
    
       ! for specific box, now start loops over alloted cells    
       do k=lo(3)-ng, hi(3)+ng
          do j=lo(2)-ng, hi(2)+ng
             do i=lo(1)-ng, hi(1)+ng
                call convert_rho_to_W_local(rho(i,j,k,:),W(i,j,k,:),Wij(i,j,k,:))
             end do
          end do
       end do
 
     end subroutine convert_rho_to_W_3d

     subroutine convert_rho_to_W_local(rho,W,Wij)
 
       real(kind=dp_t), intent(in)   :: rho(nspecies)  ! density
       real(kind=dp_t), intent(out)  ::   W(nspecies)    ! Wi (diagonals) 
       real(kind=dp_t), intent(out)  ::   Wij(nspecies,nspecies) ! Wij (full matrix) 
    
       ! local variables
       integer          :: i,j
       real(kind=dp_t)  :: rho_tot        ! total density 

       ! calculate total density inside each cell
       rho_tot=0.d0 
       do i=1, nspecies  
          rho_tot = rho_tot + rho(i)
       end do         
  
       ! calculate mass fraction 
       do i=1, nspecies  
          W(i) = rho(i)/rho_tot
       end do

       do i=1,nspecies
          do j=1, i-1
             Wij(i,j) = W(i)*W(j)
             Wij(j,i) = Wij(i,j) 
          end do
          Wij(i,i) = W(i)*W(i)
       end do 
 
     end subroutine convert_rho_to_W_local 

  end subroutine meanvar_W

end module analysis_module
