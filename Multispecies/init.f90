module init_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use probin_common_module
  use probin_multispecies_module
 
  implicit none

  private

  public :: init_rho
  
  ! init_type: 
  ! 1 = rho in concentric circle (two values inside and outside concentric circular region), 
  ! 2 = constant gradient (constant rho and spatial distortion proportional to x and y),
  ! 3 = gaussian spread with total density constant
  ! 4 = manufactured solution for equal molar mass, gaussian rho & time-independent-space-varying total density
  ! 5 = manufactured solution for unequal molar mass, gaussian rho & time-independent-space-varying total density

contains
  
  subroutine init_rho(rho,dx,prob_lo,prob_hi,time,the_bc_level)

    type(multifab) , intent(inout) :: rho(:)            
    real(kind=dp_t), intent(in   ) :: dx(:,:)           
    real(kind=dp_t), intent(in   ) :: prob_lo(rho(1)%dim)
    real(kind=dp_t), intent(in   ) :: prob_hi(rho(1)%dim)
    real(kind=dp_t), intent(in   ) :: time 
    type(bc_level) , intent(in   ) :: the_bc_level(:)
 
    ! local variables
    integer                        :: lo(rho(1)%dim), hi(rho(1)%dim)
    integer                        :: dm, ng, i, n, nlevs
    real(kind=dp_t), pointer       :: dp(:,:,:,:)   ! pointer for rho (last dim:nspecies)   

    dm = rho(1)%dim
    ng = rho(1)%ng
    nlevs = size(rho,1)

    ! looping over boxes 
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp  => dataptr(rho(n),i)
          lo  = lwb(get_box(rho(n),i))
          hi  = upb(get_box(rho(n),i))
          !print*, lo, hi 
          
          select case(dm)
          case (2)
             call init_rho_2d(dp(:,:,1,:),ng,lo,hi,prob_lo,prob_hi,dx(n,:),time)
          case (3)
             call init_rho_3d(dp(:,:,:,:),ng,lo,hi,prob_lo,prob_hi,dx(n,:),time)
          end select
       end do

       ! fill ghost cells for two adjacent grids including periodic boundary ghost cells
       call multifab_fill_boundary(rho(n))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(rho(n),1,rho_part_bc_comp,nspecies,the_bc_level(n),dx(n,:),.false.)
    end do

  end subroutine init_rho

  subroutine init_rho_2d(rho,ng,lo,hi,prob_lo,prob_hi,dx,time)

    integer          :: lo(2), hi(2), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)  ! last dimension for species
    real(kind=dp_t)  :: prob_lo(2),prob_hi(2)
    real(kind=dp_t)  :: dx(:)
    real(kind=dp_t)  :: time 
 
    ! local varables
    integer          :: i,j
    real(kind=dp_t)  :: x,y,rsq,tau,rhot,L(2)
 
    L(1:2) = prob_hi(1:2)-prob_lo(1:2) ! Domain length
    
    ! for specific box, now start loop over alloted cells     
    select case(init_type) 
    
    case(1) 
    !=========================================================================
    ! Initializing rho's in concentric circle at (Lx/2,Ly/2) with radius^2=0.1
    !=========================================================================
      do j=lo(2),hi(2)
         y = prob_lo(2) + (dble(j)+0.5d0) * dx(2) - 0.5d0
         do i=lo(1),hi(1)
            x = prob_lo(1) + (dble(i)+0.5d0) * dx(1) - 0.5d0
       
            rsq = (x-L(1)*0.5d0)**2 + (y-L(2)*0.5d0)**2
            if (rsq .lt. L(1)*L(2)*0.1d0) then
               rho(i,j,1:nspecies) = rho_in(1,1:nspecies)
            else
               rho(i,j,1:nspecies) = rho_in(2,1:nspecies)
            endif
    
         end do
      end do
  
    case(2) 
    !========================================================
    ! Initializing rho's with constant gradient for 2-species  
    !========================================================
      do j=lo(2),hi(2)
         y = prob_lo(2) + (dble(j)+0.5d0) * dx(2) - 0.5d0
         do i=lo(1),hi(1)
            x = prob_lo(1) + (dble(i)+0.5d0) * dx(1) - 0.5d0
            !rho(i,j,1:nspecies) = rho_in(1,1:nspecies)
            rho(i,j,1) = rho_in(1,1) + 0.001d0*x - 0.002d0*y
            rho(i,j,2) = rho_in(1,2) + 0.003d0*x + 0.001d0*y
            !rho(i,j,3) = rho_in(1,3) + 0.002d0*x - 0.001d0*y
         end do
      end do

    case(3) 
    !========================================================
    ! Initializing rho's in Gaussian so as rho_tot=constant=1.0
    ! Here rho_exact = e^(-r^2/4Dt)/(4piDt)
    !========================================================
  
    do j=lo(2),hi(2)
         y = prob_lo(2) + (dble(j)+0.5d0) * dx(2) - 0.5d0
         do i=lo(1),hi(1)
            x = prob_lo(1) + (dble(i)+0.5d0) * dx(1) - 0.5d0
        
            rsq = (x-L(1)*0.5d0)**2 + (y-L(2)*0.5d0)**2
            rho(i,j,1) = 1.0d0/(4.0d0*M_PI*Dbar_in(1)*time)*dexp(-rsq/(4.0d0*Dbar_in(1)*time))
            rho(i,j,2) = 1.0d0-1.0d0/(4.0d0*M_PI*Dbar_in(1)*time)*dexp(-rsq/(4.0d0*Dbar_in(1)*time))
       
         end do
    end do

    case(4)
    !=======================================================================
    ! Initializing rho1,rho2=Gaussian and rhot=space varying-constant 
    ! in time. Manufactured solution rho1_exact = e^(-r^2/4Dt-t/tau)/(4piDt)
    !=======================================================================
    tau=1.0d0 
 
    do j=lo(2),hi(2)
         y = prob_lo(2) + (dble(j)+0.5d0) * dx(2) - 0.5d0
         do i=lo(1),hi(1)
            x = prob_lo(1) + (dble(i)+0.5d0) * dx(1) - 0.5d0
        
            rsq = (x-L(1)*0.5d0)**2 + (y-L(2)*0.5d0)**2
            rhot = 1.0d0+1.0d0/(4.0d0*M_PI*Dbar_in(1))*dexp(-rsq/(4.0d0*Dbar_in(1)))
           
            rho(i,j,1) = 1.0d0/(4.0d0*M_PI*Dbar_in(1)*time)*dexp(-rsq/(4.0d0*Dbar_in(1)*time)-&
                         time/tau)*rhot

            rho(i,j,2) = (1.0d0-1.0d0/(4.0d0*M_PI*Dbar_in(1)*time)*dexp(-rsq/(4.0d0*Dbar_in(1)*time)-&
                         time/tau))*rhot
           
         end do
    end do

    end select
   
  end subroutine init_rho_2d

  subroutine init_rho_3d(rho,ng,lo,hi,prob_lo,prob_hi,dx,time)
    
    integer          :: lo(3), hi(3), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:) ! Last dimension for species 
    real(kind=dp_t)  :: prob_lo(3),prob_hi(3)
    real(kind=dp_t)  :: dx(:)
    real(kind=dp_t)  :: time 
 
    ! local variables
    integer          :: i,j,k
    real(kind=dp_t)  :: x,y,z,rsq,tau,rhot,L(3)

    L(1:3) = prob_hi(1:3)-prob_lo(1:3) ! Domain length

    ! for specific box, now start loop over alloted cells     
    select case(init_type) 
    
    case(1) 
    !================================================================================
    ! Initializing rho's in concentric circle at (Lx/2,Ly/2,Lz/2) with radius^2=0.001
    !================================================================================
    !$omp parallel private(i,j,k,x,y,z)
    do k=lo(3),hi(3)
       z = prob_lo(3) + (dble(k)+0.5d0) * dx(3) - 0.5d0
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0) * dx(2) - 0.5d0
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0) * dx(1) - 0.5d0
             
             rsq = (x-L(1)*0.5d0)**2 + (y-L(2)*0.5d0)**2 + (z-L(3)*0.5d0)**2
             if (rsq .lt. L(1)*L(2)*L(3)*0.001d0) then
                 rho(i,j,k,1:nspecies) = rho_in(1,1:nspecies)
             else
                 rho(i,j,k,1:nspecies) = rho_in(2,1:nspecies)
             endif
          
          end do
       end do
    end do
    !$omp end parallel do

    case(2) 
    !========================================================
    ! Initializing rho's with constant gradient for 2-species  
    !========================================================
    !$omp parallel private(i,j,k,x,y,z)
    do k=lo(3),hi(3)
       z = prob_lo(3) + (dble(k)+0.5d0) * dx(3) - 0.5d0
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0) * dx(2) - 0.5d0
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0) * dx(1) - 0.5d0
    
            !rho(i,j,k,1:nspecies) = rho_in(1,1:nspecies)
            rho(i,j,k,1) = rho_in(1,1) + 0.001d0*x - 0.002d0*y + 0.003d0*z
            rho(i,j,k,2) = rho_in(1,2) + 0.003d0*x + 0.001d0*y - 0.002d0*z
            !rho(i,j,k,3) = rho_in(1,3) + 0.002d0*x - 0.001d0*y
    
          end do
       end do
     end do
     !$omp end parallel do

     case(3) 
     !===========================================================
     ! Initializing rho's in Gaussian so as rho_tot=constant=1.0. 
     ! Here rho_exact = e^(-r^2/4Dt)/(4piDt)^3/2, For norm, 
     ! sigma/dx >2 (at t=0) & L/sigma < 8 (at t=t)
     !===========================================================
  
     !$omp parallel private(i,j,k,x,y,z)
     do k=lo(3),hi(3)
        z = prob_lo(3) + (dble(k)+0.5d0) * dx(3) - 0.5d0
        do j=lo(2),hi(2)
           y = prob_lo(2) + (dble(j)+0.5d0) * dx(2) - 0.5d0
           do i=lo(1),hi(1)
              x = prob_lo(1) + (dble(i)+0.5d0) * dx(1) - 0.5d0
        
              rsq = (x-L(1)*0.5d0)**2 + (y-L(2)*0.5d0)**2 + (z-L(3)*0.5d0)**2
              rho(i,j,k,1) = dexp(-rsq/(4.0d0*Dbar_in(1)*time))/(4.0d0*M_PI*Dbar_in(1)*time)**1.5d0
              rho(i,j,k,2) = 1.0d0 - dexp(-rsq/(4.0d0*Dbar_in(1)*time))/(4.0d0*M_PI*Dbar_in(1)*time)**1.5d0
       
           end do
        end do
     end do
     !$omp end parallel do

     case(4)
     !=============================================================================
     ! Initializing rho1,rho2=Gaussian and rhot=space varying-constant 
     ! in time. Manufactured solution rho1_exact = e^(-r^2/4Dt-t/tau)/(4piDt)^(3/2)
     !=============================================================================
     tau=1.0d0 
 
     !$omp parallel private(i,j,k,x,y,z)
     do k=lo(3),hi(3)
        z = prob_lo(3) + (dble(k)+0.5d0) * dx(3) - 0.5d0
        do j=lo(2),hi(2)
           y = prob_lo(2) + (dble(j)+0.5d0) * dx(2) - 0.5d0
           do i=lo(1),hi(1)
              x = prob_lo(1) + (dble(i)+0.5d0) * dx(1) - 0.5d0
        
              rsq = (x-L(1)*0.5d0)**2 + (y-L(2)*0.5d0)**2 + (z-L(3)*0.5d0)**2
              rhot = 1.0d0+dexp(-rsq/(4.0d0*Dbar_in(1)))/(4.0d0*M_PI*Dbar_in(1))**1.5d0
           
              rho(i,j,k,1) = 1.0d0/(4.0d0*M_PI*Dbar_in(1)*time)**1.5d0*dexp(-rsq/(4.0d0*Dbar_in(1)*time)-&
                           time/tau)*rhot
              rho(i,j,k,2) = (1.0d0-1.0d0/(4.0d0*M_PI*Dbar_in(1)*time)**1.5d0*dexp(-rsq/(4.0d0*Dbar_in(1)*time)-&
                           time/tau))*rhot

           end do
        end do
     end do
     !$omp end parallel do

    end select
   
  end subroutine init_rho_3d

end module init_module
