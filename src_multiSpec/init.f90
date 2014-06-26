module init_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use multifab_coefbc_module
  use ml_layout_module
  use convert_stag_module
  use probin_common_module, only: prob_lo, prob_hi, visc_coef, prob_type, &
                                  molmass, rhobar, smoothing_width, u_init, n_cells
  use probin_multispecies_module, only: alpha1, beta, delta, sigma, Dbar, &
                                        T_init, rho_init, nspecies, rho_part_bc_comp, rho_init
 
  implicit none

  private

  public :: init_rho_and_umac, init_Temp
  
  ! prob_type: If negative the density of the last species is overwritten to enforce low Mach EOS
  ! 0 = rho constant in space (thermodynamic equilibrium), temperature profile to check thermodiffusion
  ! 1 = rho in concentric circle (two values inside and outside), temperature is distributed similarly 
  ! 2 = constant gradient (spatial distortion proportional to y), temperature is distributed similarly
  ! 3 = gaussian spread with total density constant, 2-species, temperature fixed to 1
  ! 4 = manufactured solution for 2-species equal/unequal molarmass,gaussian-rho,time-independent-space-varying 
  !     totaldensity,temperature fixed to 1 
  ! 5 = manufactured solution for 3-species time-independent-space-varying density,temperature fixed to 1
  ! and more in 2D only, see below

contains
  
  subroutine init_rho_and_umac(rho,umac,dx,time,the_bc_level)

    type(multifab) , intent(inout) :: rho(:)
    type(multifab) , intent(inout) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(in   ) :: time
    type(bc_level) , intent(in   ) :: the_bc_level(:)
 
    ! local variables
    integer                        :: lo(rho(1)%dim), hi(rho(1)%dim)
    integer                        :: dm, ng_r, ng_u, i, n, nlevs
    real(kind=dp_t), pointer       :: dp(:,:,:,:)   ! pointer for rho (last dim:nspecies)   
    real(kind=dp_t), pointer       :: up(:,:,:,:)   ! pointers for mac velocities
    real(kind=dp_t), pointer       :: vp(:,:,:,:)
    real(kind=dp_t), pointer       :: wp(:,:,:,:)

    dm = rho(1)%dim
    ng_r = rho(1)%ng
    ng_u = umac(1,1)%ng
    nlevs = size(rho,1)

    ! assign values of parameters for the gaussian rho, rhototal
    alpha1 = 0.5d0 
    beta   = 0.1d0 
    delta  = 0.5d0 
    sigma  = (prob_hi(1)-prob_lo(1))/10.0d0  ! variance of gaussian distribution

    ! looping over boxes 
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp => dataptr(rho(n),i)
          up => dataptr(umac(n,1),i)
          vp => dataptr(umac(n,2),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          select case(dm)
          case (2)
             call init_rho_and_umac_2d(dp(:,:,1,:),ng_r,up(:,:,1,1),vp(:,:,1,1),ng_u, &
                                       lo,hi,dx(n,:),time)
          case (3)
             wp => dataptr(umac(n,3),i)
             call init_rho_and_umac_3d(dp(:,:,:,:),ng_r,up(:,:,:,1),vp(:,:,:,1),wp(:,:,:,1),ng_u, &
                                       lo,hi,dx(n,:),time)
          end select
       end do

       ! fill ghost cells for two adjacent grids including periodic boundary ghost cells
       call multifab_fill_boundary(rho(n))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(rho(n), 1,rho_part_bc_comp,nspecies,the_bc_level(n),dx(n,:))

       ! Note: not filling umac ghost cells here.

    end do

  end subroutine init_rho_and_umac

  subroutine init_rho_and_umac_2d(rho,ng_r,u,v,ng_u,lo,hi,dx,time)

    integer          :: lo(2), hi(2), ng_r, ng_u
    real(kind=dp_t)  :: rho(lo(1)-ng_r:,lo(2)-ng_r:,:)  ! last dimension for species
    real(kind=dp_t)  ::   u(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t)  ::   v(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t)  :: dx(:)
    real(kind=dp_t)  :: time 
 
    ! local varables
    integer          :: i,j,n
    real(kind=dp_t)  :: x,y,w1,w2,rsq,rhot,L(2),sum,r,y1,rho_loc,rand
 
    L(1:2) = prob_hi(1:2)-prob_lo(1:2) ! Domain length
    
    ! for specific box, now start loop over alloted cells     
    ! if prob_type is negative, we enforce low mach constraint by overwriting final rho_i
    select case (abs(prob_type))

    case(0) 
    !============================================================
    ! Thermodynamic equilibrium
    !============================================================
 
    u = 0.d0
    v = 0.d0

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
 
          rho(i,j,1:nspecies) = rho_init(1,1:nspecies)

       end do
    end do  
    
    case(1) 
    !=============================================================
    ! Initializing rho's in concentric circle with radius^2 = 0.1
    !=============================================================
 
    u = 0.d0
    v = 0.d0

    do j=lo(2),hi(2)
       y = prob_lo(2) + (dble(j)+half)*dx(2) - half*(prob_lo(2)+prob_hi(2))
       do i=lo(1),hi(1)
          x = prob_lo(1) + (dble(i)+half)*dx(1) - half*(prob_lo(1)+prob_hi(1))
       
          rsq = x**2 + y**2
          if (rsq .lt. L(1)*L(2)*0.1d0) then
             rho(i,j,1:nspecies) = rho_init(1,1:nspecies)
          else
             rho(i,j,1:nspecies) = rho_init(2,1:nspecies)
          end if
    
       end do
    end do
  
    case(2) 
    !=========================================================
    ! Initializing rho's with constant gradient 
    !=========================================================
 
    u = 0.d0
    v = 0.d0

    do j=lo(2),hi(2)
       y = prob_lo(2) + (dble(j)+half)*dx(2) 
       do i=lo(1),hi(1)
          x = prob_lo(1) + (dble(i)+half)*dx(1) 
   
            ! linear gradient in rho
            rho(i,j,1:nspecies) = rho_init(1,1:nspecies) + & 
               (rho_init(2,1:nspecies) - rho_init(1,1:nspecies))*(y-prob_lo(2))/L(2)
   
         end do
      end do

    case(3) 
    !===========================================================
    ! Initializing rho's in Gaussian so as rhotot=constant=1.0
    ! Here rho_exact = e^(-r^2/4Dt)/(4piDt)
    !===========================================================
 
    u = 0.d0
    v = 0.d0

    do j=lo(2),hi(2)
         y = prob_lo(2) + (dble(j)+half) * dx(2) - half*(prob_lo(2)+prob_hi(2))
         do i=lo(1),hi(1)
            x = prob_lo(1) + (dble(i)+half) * dx(1) - half*(prob_lo(1)+prob_hi(1))
        
            rsq = x**2 + y**2
            rho(i,j,1) = 1.0d0/(4.0d0*M_PI*Dbar(1)*time)*dexp(-rsq/(4.0d0*Dbar(1)*time))
            rho(i,j,2) = 1.0d0-1.0d0/(4.0d0*M_PI*Dbar(1)*time)*dexp(-rsq/(4.0d0*Dbar(1)*time))
       
         end do
      end do

    case(4)
    !==================================================================================
    ! Initializing rho1,rho2=Gaussian and rhototal=1+alpha*exp(-r^2/4D)/(4piD) (no-time 
    ! dependence). Manufactured solution rho1_exact = exp(-r^2/4Dt-beta*t)/(4piDt)
    !==================================================================================
 
    u = 0.d0
    v = 0.d0

    do j=lo(2),hi(2)
         y = prob_lo(2) + (dble(j)+half) * dx(2) - half
         do i=lo(1),hi(1)
            x = prob_lo(1) + (dble(i)+half) * dx(1) - half
        
            rsq = (x-L(1)*half)**2 + (y-L(2)*half)**2
            rhot = 1.0d0 + alpha1/(4.0d0*M_PI*Dbar(1))*dexp(-rsq/(4.0d0*Dbar(1)))
            rho(i,j,1) = 1.0d0/(4.0d0*M_PI*Dbar(1)*time)*dexp(-rsq/(4.0d0*Dbar(1)*time)-&
                         beta*time)*rhot
            rho(i,j,2) = rhot - rho(i,j,1)

         end do
      end do

    case(5)
    !==================================================================================
    ! Initializing m2=m3, D12=D13 where Dbar(1)=D12, Dbar(2)=D13, 
    ! Dbar(3)=D23, Grad(w2)=0, manufactured solution for rho1 and rho2 
    ! (to benchmark eqn1) Initializing rho1, rho2=Gaussian and rhototal has no-time dependence.
    !==================================================================================
 
    u = 0.d0
    v = 0.d0

    do j=lo(2),hi(2)
         y = prob_lo(2) + (dble(j)+half) * dx(2) - half
         do i=lo(1),hi(1)
            x = prob_lo(1) + (dble(i)+half) * dx(1) - half
        
            rsq = (x-L(1)*half)**2 + (y-L(2)*half)**2
            w1  = alpha1*dexp(-rsq/(2.0d0*sigma**2))
            w2  =  delta*dexp(-beta*time)
            rhot = 1.0d0 + (molmass(2)*Dbar(3)/(molmass(1)*Dbar(1))-1.0d0)*w1
            rho(i,j,1) = rhot*w1
            rho(i,j,2) = rhot*w2 
            rho(i,j,3) = rhot-rho(i,j,1)-rho(i,j,2)
           
         end do
    end do

    case(6) 
    !=========================================================
    ! Test of thermodiffusion steady-state for 2 species 
    !=========================================================
 
    u = 0.d0
    v = 0.d0

    do j=lo(2),hi(2)
       y = prob_lo(2) + (dble(j)+half)*dx(2) 
       do i=lo(1),hi(1)
          x = prob_lo(1) + (dble(i)+half)*dx(1) 
   
            ! Solution to diff(c(y),y)=-K*c(y)*(1-c(y))
            ! Here K=grad(T)*S_T=0.15
            ! Height of domain H=32
            ! And average <rho1>=.4830852506
            rho(i,j,1) = 1.0d0/(1.0d0+0.1d0*exp(0.15d0*y))
            rho(i,j,2) = 1.0d0 - rho(i,j,1) 
   
         end do
      end do

   case(7)

    !=============================================================
    ! smoothed circle
    !=============================================================
 
    u = 0.d0
    v = 0.d0

    do j=lo(2),hi(2)
       y = prob_lo(2) + (dble(j)+half)*dx(2) - half*(prob_lo(2)+prob_hi(2))
       do i=lo(1),hi(1)
          x = prob_lo(1) + (dble(i)+half)*dx(1) - half*(prob_lo(1)+prob_hi(1))
       
          r = sqrt(x**2 + y**2)

          rho(i,j,1:nspecies-1) = rho_init(1,1:nspecies-1) + &
               0.5d0*(rho_init(2,1:nspecies-1) - rho_init(1,1:nspecies-1))* &
                  (1.d0 + tanh((r-15.d0)/2.d0))

       end do
    end do

    case(8)

    !=============================================================
    ! 4-species, 4-stripes
    !=============================================================
 
    u = 0.d0
    v = 0.d0

    rho = 0.d0
           
    do j=lo(2),hi(2)
       y = prob_lo(2) + (dble(j)+half)*dx(2) 
       do i=lo(1),hi(1)
          if (y .le. 0.75d0*prob_lo(2)+0.25d0*prob_hi(2)) then
             rho(i,j,1) = rho_init(1,1)
             rho(i,j,2) = 0.1d0*rhobar(2)
             rho(i,j,3) = 0.1d0*rhobar(3)
             rho(i,j,4) = 0.1d0*rhobar(4)
          else if (y .le. 0.5d0*prob_lo(2)+0.5d0*prob_hi(2)) then
             rho(i,j,1) = 0.1d0*rhobar(1)
             rho(i,j,2) = rho_init(1,2)
             rho(i,j,3) = 0.1d0*rhobar(3)
             rho(i,j,4) = 0.1d0*rhobar(4)
          else if (y .le. 0.25d0*prob_lo(2)+0.75d0*prob_hi(2)) then
             rho(i,j,1) = 0.1d0*rhobar(1)
             rho(i,j,2) = 0.1d0*rhobar(2)
             rho(i,j,3) = rho_init(1,3)
             rho(i,j,4) = 0.1d0*rhobar(4)
          else
             rho(i,j,1) = 0.1d0*rhobar(1)
             rho(i,j,2) = 0.1d0*rhobar(2)
             rho(i,j,3) = 0.1d0*rhobar(3)
             rho(i,j,4) = rho_init(1,4)
          end if
       end do
    end do

    case(9)

    !=============================================================
    ! one fluid on top of another
    ! rho1 = rho_init(1) in lower half of domain (in y)
    ! rho1 = rho_init(2) in upper half
    ! use the EOS to compute rho2 by setting prob_type negative
    !=============================================================
 
    u = 0.d0
    v = 0.d0

    ! middle of domain
    y1 = (prob_lo(2)+prob_hi(2)) / 2.d0

    if(abs(smoothing_width)>epsilon(1.d0)) then

       ! smoothed version
       do j=lo(2),hi(2)
          y = prob_lo(2) + dx(2)*(dble(j)+0.5d0) - y1
          
          rho_loc = rho_init(1,1) + (rho_init(1,2)-rho_init(1,1))*0.5d0*(tanh(y/(smoothing_width*dx(2)))+1.d0)
          rho(lo(1):hi(1),j,1) = rho_loc

       end do

    else

       ! discontinuous version
       do j=lo(2),hi(2)
          y = prob_lo(2) + (j+0.5d0)*dx(2)
          if (y .lt. y1) then
             rho(lo(1):hi(1),j,1) = rho_init(1,1)
          else
             rho(lo(1):hi(1),j,1) = rho_init(1,2)
          end if
       end do

    end if

    case (10)

    !=============================================================
    ! low Mach Kelvin-Helmholtz comparison to binary version
    ! one fluid on top of another
    ! discontinuous interface, but with random density perturbation added 
    ! in a 1-cell thick transition region
    !=============================================================

    v = 0.d0

    ! middle of domain
    y1 = (prob_lo(2)+prob_hi(2)) / 2.d0

    ! rho1 = rho_init(1,1) in lower half of domain (in y)
    ! rho1 = rho_init(2,1) in upper half
    ! random perturbation below centerline

    do j=lo(2),hi(2)
       y = prob_lo(2) + (j+0.5d0)*dx(2)
          
       if (y .lt. y1) then
          rho_loc = rho_init(1,1)
       else
          rho_loc = rho_init(2,1)
       end if

       rho(lo(1):hi(1),j,1) = rho_loc
       rho(lo(1):hi(1),j,2) = (1.d0 - rho_loc/rhobar(1))*rhobar(2)

       ! add random perturbation below centerline
       if (j .eq. n_cells(2)/2-1) then
          do i=lo(1),hi(1)
             call random_number(rand)
             rho_loc = rand*rho_init(1,1) + (1.d0-rand)*rho_init(2,1)
             rho(i,j,1) = rho_loc
             rho(i,j,2) = (1.d0 - rho_loc/rhobar(1))*rhobar(2)
          end do
       end if
          
    end do
       
    ! velocity = u_init(1) below centerline
    !            u_init(2) above centerline
    do j=lo(2),hi(2)
       y = prob_lo(2) + (j+0.5d0)*dx(2)

       if (y .lt. y1) then
          u(:,j) = u_init(1)
       else
          u(:,j) = u_init(2)
       end if

    end do

    case (11)

    !=============================================================
    ! 1 fluid on top of another
    ! rho(:) = rho_init(1,:) on bottom
    ! rho(:) = rho_init(2,:) on top
    !=============================================================

    u = 0.d0
    v = 0.d0

    ! middle of domain
    y1 = (prob_lo(2)+prob_hi(2)) / 2.d0

    ! rho1 = rho_init(1,1) in lower half of domain (in y)
    ! rho1 = rho_init(2,1) in upper half
    ! random perturbation below centerline

    do j=lo(2),hi(2)
       y = prob_lo(2) + (j+0.5d0)*dx(2)
          
       if (y .lt. y1) then
          do i=lo(1),hi(1)
             rho(i,j,1:nspecies) = rho_init(1,1:nspecies)
          end do
       else
          do i=lo(1),hi(1)
             rho(i,j,1:nspecies) = rho_init(2,1:nspecies)
          end do
       end if

    end do

    case default
      
      call bl_error("Desired prob_type not supported in 3D")
      
    end select

    if (prob_type .lt. 0) then

       ! enforce low mach constraint by overwriting final rho_i
       do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          sum = 0
          do n=1,nspecies-1
             sum = sum + rho(i,j,n)/rhobar(n)
          end do
          rho(i,j,nspecies) = rhobar(nspecies)*(1.d0 - sum)

       end do
       end do          

    end if
   
  end subroutine init_rho_and_umac_2d

  subroutine init_rho_and_umac_3d(rho,ng_r,u,v,w,ng_u,lo,hi,dx,time)
    
    integer          :: lo(3), hi(3), ng_r, ng_u
    real(kind=dp_t)  :: rho(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:,:) ! Last dimension for species 
    real(kind=dp_t)  ::   u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t)  ::   v(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t)  ::   w(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t)  :: dx(:)
    real(kind=dp_t)  :: time 
 
    ! local variables
    integer          :: i,j,k,n
    real(kind=dp_t)  :: x,y,z,rsq,w1,w2,rhot,L(3),sum,rand,rho_loc,y1,r

    L(1:3) = prob_hi(1:3)-prob_lo(1:3) ! Domain length

    ! for specific box, now start loop over alloted cells
    ! if prob_type is negative, we enforce low mach constraint by overwriting final rho_i
    select case (abs(prob_type))
    
    case(0) 
    !================================================================================
    ! Thermodynamic equilibrium
    !================================================================================
 
    u = 0.d0
    v = 0.d0
    w = 0.d0

    !$omp parallel do private(i,j,k)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             rho(i,j,k,1:nspecies) = rho_init(1,1:nspecies)

          end do
       end do
    end do
    !$omp end parallel do
    
    case(1) 
    !================================================================================
    ! Initializing rho's in concentric circle 
    !================================================================================
 
    u = 0.d0
    v = 0.d0
    w = 0.d0
  
    !$omp parallel do private(i,j,k,x,y,z,rsq)
    do k=lo(3),hi(3)
       z = prob_lo(3) + (dble(k)+half)*dx(3) - half*(prob_lo(3)+prob_hi(3))
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+half)*dx(2) - half*(prob_lo(2)+prob_hi(2))
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+half)*dx(1) - half*(prob_lo(1)+prob_hi(1))

               rsq = x**2 + y**2 + z**2
               if (rsq .lt. L(1)*L(2)*L(3)*0.001d0) then
                  rho(i,j,k,1:nspecies) = rho_init(1,1:nspecies)
               else
                  rho(i,j,k,1:nspecies) = rho_init(2,1:nspecies)
               end if
          
          end do
       end do
    end do
    !$omp end parallel do

    case(2) 
    !========================================================
    ! Initializing rho's with constant gradient
    !========================================================
 
    u = 0.d0
    v = 0.d0
    w = 0.d0
 
    !$omp parallel do private(i,j,k,x,y,z)
    do k=lo(3),hi(3)
       z = prob_lo(3) + (dble(k)+half)*dx(3) 
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+half)*dx(2)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+half)*dx(1)

               rho(i,j,k,1:nspecies) = rho_init(1,1:nspecies) + &
                   (rho_init(2,1:nspecies) - rho_init(1,1:nspecies))*(y-prob_lo(2))/L(2)

          end do
       end do
     end do
     !$omp end parallel do

     case(3) 
     !================================================================================
     ! Initializing rho's in Gaussian so as rhotot=constant=1.0. Here rho_exact = 
     ! e^(-r^2/4Dt)/(4piDt)^3/2, For norm, sigma/dx >2 (at t=0) & L/sigma < 8 (at t=t)
     !================================================================================
 
     u = 0.d0
     v = 0.d0
     w = 0.d0
  
     !$omp parallel do private(i,j,k,x,y,z,rsq)
     do k=lo(3),hi(3)
        z = prob_lo(3) + (dble(k)+half)*dx(3) - half*(prob_lo(3)+prob_hi(3))
        do j=lo(2),hi(2)
           y = prob_lo(2) + (dble(j)+half)*dx(2) - half*(prob_lo(2)+prob_hi(2))
           do i=lo(1),hi(1)
              x = prob_lo(1) + (dble(i)+half)*dx(1) - half*(prob_lo(1)+prob_hi(1))
        
              rsq = x**2 + y**2 + z**2
              rho(i,j,k,1) = dexp(-rsq/(4.0d0*Dbar(1)*time))/(4.0d0*M_PI*&
                             Dbar(1)*time)**1.5d0
              rho(i,j,k,2) = 1.0d0 - dexp(-rsq/(4.0d0*Dbar(1)*time))/(4.0d0*&
                             M_PI*Dbar(1)*time)**1.5d0
       
           end do
        end do
     end do
     !$omp end parallel do

     case(4)
     !==============================================================================
     ! Initializing rho1,rho2=Gaussian and rhot=space varying-constant 
     ! in time. Manufactured solution rho1_exact = e^(-r^2/4Dt-t/tau)/(4piDt)^(3/2)
     !==============================================================================
 
     u = 0.d0
     v = 0.d0
     w = 0.d0
     
     !$omp parallel do private(i,j,k,x,y,z,rsq,rhot)
     do k=lo(3),hi(3)
        z = prob_lo(3) + (dble(k)+half) * dx(3) - half
        do j=lo(2),hi(2)
           y = prob_lo(2) + (dble(j)+half) * dx(2) - half
           do i=lo(1),hi(1)
              x = prob_lo(1) + (dble(i)+half) * dx(1) - half
        
              rsq = (x-L(1)*half)**2 + (y-L(2)*half)**2 + (z-L(3)*half)**2
              rhot = 1.0d0 + alpha1*dexp(-rsq/(4.0d0*Dbar(1)))/(4.0d0*M_PI*&
                     Dbar(1))**1.5d0
           
              rho(i,j,k,1) = 1.0d0/(4.0d0*M_PI*Dbar(1)*time)**1.5d0*dexp(-rsq/&
                             (4.0d0*Dbar(1)*time) - time*beta)*rhot
              rho(i,j,k,2) = rhot - rho(i,j,k,1) 

           end do
        end do
     end do
     !$omp end parallel do

     case(5)
     !==================================================================================
     ! Initializing m2=m3, D12=D13 where Dbar(1)=D12, Dbar(2)=D13, Dbar(3)=D23, 
     ! Grad(w2)=0, manufactured solution for rho1 and rho2 
     !==================================================================================
 
     u = 0.d0
     v = 0.d0
     w = 0.d0

     !$omp parallel do private(i,j,k,x,y,z,rsq,w1,w2,rhot)
     do k=lo(3),hi(3)
        z = prob_lo(3) + (dble(k)+half) * dx(3) - half
        do j=lo(2),hi(2)
           y = prob_lo(2) + (dble(j)+half) * dx(2) - half
           do i=lo(1),hi(1)
              x = prob_lo(1) + (dble(i)+half) * dx(1) - half

              rsq = (x-L(1)*half)**2 + (y-L(2)*half)**2 + (z-L(3)*half)**2
              w1  = alpha1*dexp(-rsq/(2.0d0*sigma**2))
              w2  =  delta*dexp(-beta*time)
              rhot = 1.0d0 + (molmass(2)*Dbar(3)/(molmass(1)*Dbar(1))-1.0d0)*w1
              rho(i,j,k,1) = rhot*w1
              rho(i,j,k,2) = rhot*w2
              rho(i,j,k,3) = rhot-rho(i,j,k,1)-rho(i,j,k,2)
           
              if(rho(i,j,k,1).lt.0.d0 .or. rho(i,j,k,2).lt.0.d0 .or. rho(i,j,k,3).lt.0.d0) then 
                 write(*,*), "rho1 / rho2 / rho3 is negative: STOP"
                 write(*,*), i, j, " w1=", w1, " w2=", w2, " rho1=",rho(i,j,k,1)," rho2=",&
                             rho(i,j,k,2), " rho3=",rho(i,j,k,3), " rhot=",rhot
              end if
 
          end do
       end do
    end do
    !$omp end parallel do

   case(7)

    !=============================================================
    ! smoothed circle
    !=============================================================
 
    u = 0.d0
    v = 0.d0
    w = 0.d0

    do k=lo(3),hi(3)
       z = prob_lo(3) + (dble(k)+half)*dx(3) - half*(prob_lo(3)+prob_hi(3))
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+half)*dx(2) - half*(prob_lo(2)+prob_hi(2))
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+half)*dx(1) - half*(prob_lo(1)+prob_hi(1))
       
             r = sqrt(x**2 + y**2 + z**2)

             rho(i,j,k,1:nspecies-1) = rho_init(1,1:nspecies-1) + &
                  0.5d0*(rho_init(2,1:nspecies-1) - rho_init(1,1:nspecies-1))* &
                  (1.d0 + tanh((r-15.d0)/2.d0))

          end do
       end do
    end do

    case(8)

    !=============================================================
    ! 4-species, 4-stripes
    !=============================================================
 
    u = 0.d0
    v = 0.d0
    w = 0.d0

    rho = 0.d0
       
    do k=lo(3),hi(3)
    do j=lo(2),hi(2)
       y = prob_lo(2) + (dble(j)+half)*dx(2) 
       do i=lo(1),hi(1)
          if (y .le. 0.75d0*prob_lo(2)+0.25d0*prob_hi(2)) then
             rho(i,j,k,1) = rho_init(1,1)
             rho(i,j,k,2) = 0.1d0*rhobar(2)
             rho(i,j,k,3) = 0.1d0*rhobar(3)
             rho(i,j,k,4) = 0.1d0*rhobar(4)
          else if (y .le. 0.5d0*prob_lo(2)+0.5d0*prob_hi(2)) then
             rho(i,j,k,1) = 0.1d0*rhobar(1)
             rho(i,j,k,2) = rho_init(1,2)
             rho(i,j,k,3) = 0.1d0*rhobar(3)
             rho(i,j,k,4) = 0.1d0*rhobar(4)
          else if (y .le. 0.25d0*prob_lo(2)+0.75d0*prob_hi(2)) then
             rho(i,j,k,1) = 0.1d0*rhobar(1)
             rho(i,j,k,2) = 0.1d0*rhobar(2)
             rho(i,j,k,3) = rho_init(1,3)
             rho(i,j,k,4) = 0.1d0*rhobar(4)
          else
             rho(i,j,k,1) = 0.1d0*rhobar(1)
             rho(i,j,k,2) = 0.1d0*rhobar(2)
             rho(i,j,k,3) = 0.1d0*rhobar(3)
             rho(i,j,k,4) = rho_init(1,4)
          end if
       end do
    end do
    end do

    case(9)

    !=============================================================
    ! one fluid on top of another
    ! rho1 = rho_init(1) in lower half of domain (in y)
    ! rho1 = rho_init(2) in upper half
    ! use the EOS to compute rho2 by setting prob_type negative
    !=============================================================
 
    u = 0.d0
    v = 0.d0
    w = 0.d0

    ! middle of domain
    y1 = (prob_lo(2)+prob_hi(2)) / 2.d0

    if(abs(smoothing_width)>epsilon(1.d0)) then

       ! smoothed version
       do j=lo(2),hi(2)
          y = prob_lo(2) + dx(2)*(dble(j)+0.5d0) - y1
          
          rho_loc = rho_init(1,1) + (rho_init(1,2)-rho_init(1,1))*0.5d0*(tanh(y/(smoothing_width*dx(2)))+1.d0)
          rho(lo(1):hi(1),j,lo(3):hi(3),1) = rho_loc

       end do

    else

       ! discontinuous version
       do j=lo(2),hi(2)
          y = prob_lo(2) + (j+0.5d0)*dx(2)
          if (y .lt. y1) then
             rho(lo(1):hi(1),j,lo(3):hi(3),1) = rho_init(1,1)
          else
             rho(lo(1):hi(1),j,lo(3):hi(3),1) = rho_init(1,2)
          end if
       end do

    end if

    case (10)

    !=============================================================
    ! low Mach Kelvin-Helmholtz comparison to binary version
    ! one fluid on top of another
    ! discontinuous interface, but with random density perturbation added 
    ! in a 1-cell thick transition region
    !=============================================================

    v = 0.d0
    w = 0.d0

    ! middle of domain
    y1 = (prob_lo(2)+prob_hi(2)) / 2.d0

    ! rho1 = rho_init(1,1) in lower half of domain (in y)
    ! rho1 = rho_init(2,1) in upper half
    ! random perturbation below centerline

    do j=lo(2),hi(2)
       y = prob_lo(2) + (j+0.5d0)*dx(2)
          
       if (y .lt. y1) then
          rho_loc = rho_init(1,1)
       else
          rho_loc = rho_init(2,1)
       end if

       rho(lo(1):hi(1),j,lo(3):hi(3),1) = rho_loc
       rho(lo(1):hi(1),j,lo(3):hi(3),2) = (1.d0 - rho_loc/rhobar(1))*rhobar(2)

       ! add random perturbation below centerline
       if (j .eq. n_cells(2)/2-1) then
          do k=lo(3),hi(3)
          do i=lo(1),hi(1)
             call random_number(rand)
             rho_loc = rand*rho_init(1,1) + (1.d0-rand)*rho_init(2,1)
             rho(i,j,k,1) = rho_loc
             rho(i,j,k,2) = (1.d0 - rho_loc/rhobar(1))*rhobar(2)
          end do
          end do
       end if
          
    end do
       
    ! velocity = u_init(1) below centerline
    !            u_init(2) above centerline
    do j=lo(2),hi(2)
       y = prob_lo(2) + (j+0.5d0)*dx(2)

       if (y .lt. y1) then
          u(:,j,:) = u_init(1)
       else
          u(:,j,:) = u_init(2)
       end if

    end do

    case (11)

    !=============================================================
    ! 1 fluid on top of another
    ! rho(:) = rho_init(1,:) on bottom
    ! rho(:) = rho_init(2,:) on top
    !=============================================================

    u = 0.d0
    v = 0.d0
    w = 0.d0

    ! middle of domain
    y1 = (prob_lo(2)+prob_hi(2)) / 2.d0

    ! rho1 = rho_init(1,1) in lower half of domain (in y)
    ! rho1 = rho_init(2,1) in upper half
    ! random perturbation below centerline

    do j=lo(2),hi(2)
       y = prob_lo(2) + (j+0.5d0)*dx(2)
          
       if (y .lt. y1) then
          do k=lo(3),hi(3)
          do i=lo(1),hi(1)
             rho(i,j,k,1:nspecies) = rho_init(1,1:nspecies)
          end do
          end do
       else
          do k=lo(3),hi(3)
          do i=lo(1),hi(1)
             rho(i,j,k,1:nspecies) = rho_init(2,1:nspecies)
          end do
          end do
       end if

    end do

    case default
      
      call bl_error("Desired prob_type not supported in 3D")
      
    end select

    if (prob_type .lt. 0) then

       ! enforce low mach constraint by overwriting final rho_i
       do k=lo(3),hi(3)
       do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          sum = 0
          do n=1,nspecies-1
             sum = sum + rho(i,j,k,n)/rhobar(n)
          end do
          rho(i,j,k,nspecies) = rhobar(nspecies)*(1.d0 - sum)

       end do
       end do
       end do

    end if
   
  end subroutine init_rho_and_umac_3d

  subroutine init_Temp(Temp,dx,time,the_bc_level)

    type(multifab) , intent(inout) :: Temp(:)            
    real(kind=dp_t), intent(in   ) :: dx(:,:)           
    real(kind=dp_t), intent(in   ) :: time 
    type(bc_level) , intent(in   ) :: the_bc_level(:)
 
    ! local variables
    integer                        :: lo(Temp(1)%dim), hi(Temp(1)%dim)
    integer                        :: dm, ng, i, n, nlevs
    real(kind=dp_t), pointer       :: dp1(:,:,:,:)  ! pointer for Temp 

    dm = Temp(1)%dim
    ng = Temp(1)%ng
    nlevs = size(Temp,1)

    sigma  = (prob_hi(1)-prob_lo(1))/10.0d0  ! variance of gaussian distribution

    ! looping over boxes 
    do n=1,nlevs
       do i=1,nfabs(Temp(n))
          dp1 => dataptr(Temp(n),i)
          lo  = lwb(get_box(Temp(n),i))
          hi  = upb(get_box(Temp(n),i))
          
          select case(dm)
          case (2)
             call init_Temp_2d(dp1(:,:,1,1),ng,lo,hi,dx(n,:),time)
          case (3)
             call init_Temp_3d(dp1(:,:,:,1),ng,lo,hi,dx(n,:),time)
          end select
       end do

       ! fill ghost cells for two adjacent grids including periodic boundary ghost cells
       call multifab_fill_boundary(Temp(n))

       ! fill non-periodic domain boundary ghost cells
       call multifab_coefbc(Temp(n),1,1,the_bc_level(n))

    end do

  end subroutine init_Temp

  subroutine init_Temp_2d(Temp,ng,lo,hi,dx,time)

    integer          :: lo(2), hi(2), ng
    real(kind=dp_t)  :: Temp(lo(1)-ng:,lo(2)-ng:)  
    real(kind=dp_t)  :: dx(:)
    real(kind=dp_t)  :: time 
 
    ! local varables
    integer          :: i,j
    real(kind=dp_t)  :: x,y,rsq,L(2)
 
    L(1:2) = prob_hi(1:2)-prob_lo(1:2) ! Domain length
    
    ! for specific box, now start loop over alloted cells     
    select case (abs(prob_type)) 

    case(0) 
    !=========================================================================
    ! Thermodynamic equilibrium
    !=========================================================================
    do j=lo(2)-ng,hi(2)+ng
       y = prob_lo(2) + (dble(j)+half)*dx(2) 
       do i=lo(1)-ng,hi(1)+ng
          x = prob_lo(1) + (dble(i)+half)*dx(1) 
 
            Temp(i,j) = T_init(1)
            ! These were used for testing thermodiffusion only
            if(.false.) Temp(i,j) = 1.0d0 + 0.01d0*cos(2.0d0*M_PI*x/L(1))*sin(2.0d0*M_PI*y/L(1))
            if(.false.) Temp(i,j) = 1.0d0 + 0.01d0*sin(2.0d0*M_PI*y/L(1))

       end do
    end do  
    
    case(1) 
    !=========================================================================
    ! Initializing T in concentric circle at (Lx/2,Ly/2) with radius^2=0.1*L(1)*L(2)
    !=========================================================================
    do j=lo(2)-ng,hi(2)+ng
       y = prob_lo(2) + (dble(j)+half)*dx(2) - half*(prob_lo(2)+prob_hi(2))
       do i=lo(1)-ng,hi(1)+ng
          x = prob_lo(1) + (dble(i)+half)*dx(1) - half*(prob_lo(1)+prob_hi(1))
      
          ! temperature distribution follows the density 
          rsq = x**2 + y**2
          if (rsq .lt. L(1)*L(2)*0.1d0) then
              Temp(i,j) = T_init(1)
          else
              Temp(i,j) = T_init(2)
          end if
    
        end do
    end do
  
    case(2,6) 
    !========================================================
    ! Initializing T with constant gradient along y axes
    !========================================================
    do j=lo(2)-ng,hi(2)+ng
       y = prob_lo(2) + (dble(j)+half)*dx(2) 
       do i=lo(1)-ng,hi(1)+ng
          x = prob_lo(1) + (dble(i)+half)*dx(1) 
      
          ! linear gradient in y direction
          Temp(i,j) = T_init(1) + (T_init(2) - T_init(1))*(y-prob_lo(2))/L(2)

       end do
    end do

    case default

    do j=lo(2)-ng,hi(2)+ng
         y = prob_lo(2) + (dble(j)+half) * dx(2) - half
         do i=lo(1)-ng,hi(1)+ng
            x = prob_lo(1) + (dble(i)+half) * dx(1) - half
        
            Temp(i,j)  = T_init(1) 

         end do
    end do

   end select
   
  end subroutine init_Temp_2d

  subroutine init_Temp_3d(Temp,ng,lo,hi,dx,time)
    
    integer          :: lo(3), hi(3), ng
    real(kind=dp_t)  :: Temp(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)  
    real(kind=dp_t)  :: dx(:)
    real(kind=dp_t)  :: time 
 
    ! local variables
    integer          :: i,j,k
    real(kind=dp_t)  :: x,y,z,rsq,L(3)

    L(1:3) = prob_hi(1:3)-prob_lo(1:3) ! Domain length

    ! for specific box, now start loop over alloted cells     
    select case (abs(prob_type)) 
    
    case(0) 
    !================================================================================
    ! Thermodynamic equilibrium
    !================================================================================
 
    !$omp parallel do private(i,j,k,x,y,z)
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             
             Temp(i,j,k) = T_init(1)

          end do
       end do
    end do
    !$omp end parallel do
    
    case(1) 
    !================================================================================
    ! Initializing temperature in concentric circle 
    !================================================================================
  
    !$omp parallel do private(i,j,k,x,y,z,rsq)
    do k=lo(3)-ng,hi(3)+ng
       z = prob_lo(3) + (dble(k)+half)*dx(3) - half*(prob_lo(3)+prob_hi(3))
       do j=lo(2)-ng,hi(2)+ng
          y = prob_lo(2) + (dble(j)+half)*dx(2) - half*(prob_lo(2)+prob_hi(2))
          do i=lo(1)-ng,hi(1)+ng
             x = prob_lo(1) + (dble(i)+half)*dx(1) - half*(prob_lo(1)+prob_hi(1))
             
             !temperature distribution follows the density 
             rsq = x**2 + y**2 + z**2
             if (rsq .lt. L(1)*L(2)*L(3)*0.001d0) then
                Temp(i,j,k) = T_init(1)
             else
                Temp(i,j,k) = T_init(2)
             end if
          
          end do
       end do
    end do
    !$omp end parallel do

    case(2,6) 
    !========================================================
    ! Initializing T with constant gradient along y direction
    !========================================================
 
    !$omp parallel do private(i,j,k,x,y,z)
    do k=lo(3)-ng,hi(3)+ng
       z = prob_lo(3) + (dble(k)+half)*dx(3) 
       do j=lo(2)-ng,hi(2)+ng
          y = prob_lo(2) + (dble(j)+half)*dx(2) 
          do i=lo(1)-ng,hi(1)+ng
             x = prob_lo(1) + (dble(i)+half)*dx(1) 

             ! linear gradient in y direction for temperature 
             Temp(i,j,k) = T_init(1) + (T_init(2) - T_init(1))*(y-prob_lo(2))/L(2)
 
          end do
       end do
    end do
    !$omp end parallel do

    case default

    !$omp parallel do private(i,j,k)
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng

             Temp(i,j,k) = T_init(1)

          end do
       end do
    end do
    !$omp end parallel do

   end select
   
  end subroutine init_Temp_3d

end module init_module
