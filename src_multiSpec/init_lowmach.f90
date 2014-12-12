module init_lowmach_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use multifab_coefbc_module
  use ml_layout_module
  use convert_stag_module
  use convert_variables_module
  use probin_common_module, only: prob_lo, prob_hi, prob_type, k_B, grav, &
                                  molmass, rhobar, smoothing_width, u_init, n_cells
  use probin_multispecies_module, only: alpha1, beta, delta, sigma, Dbar, Dtherm, &
                                        rho_init, nspecies, rho_part_bc_comp, &
                                        T_init
 
  implicit none

  private

  public :: init_rho_and_umac  ! used in low Mach code; initialize c first then convert to rho

  ! IMPORTANT: In the diffusion only code (init_rho), rho_init specifies initial values for DENSITY
  ! In the low-Mach code (init_rho_and_umac), rho_init specifies initial MASS FRACTIONS
  ! (should sum to unity!... but we overwrite the final concentration so sum(c_i)=1 before computing rho)
  ! The density follows from the EOS in the LM case so it cannot be specified
  ! Same applies to boundary conditions

  ! prob_types codes for init_lowmach:

  !=============================================================
  ! case 1:
  ! bubble with radius = 1/4 of domain in x
  ! c=rho_init(1,:) inside, c=rho_init(2,:) outside
  ! can be discontinous or smooth depending on smoothing_width

  !=========================================================
  ! case 2:
  ! constant concentration gradient along y
  ! c=rho_init(1,:) on bottom, c=rho_init(2,:) on top

  !=========================================================
  ! case 3:
  ! 1 fluid on top of another
  ! c = rho_init(1,:) on bottom; c = rho_init(2,:) on top
  ! smoothing_width > 0 is a tanh smoothed interface where smoothing width is approx the # of grid 
  !   cells and then c = rand*rho_init(1,:) + (1-rand)*rho_init(2,:)
  ! smoothing_width between 0 and -1 is random perturbation where rand = abs(smoothing_width)*rand()
  ! smoothing width of -2 is a sinusoidal perturbation
  ! x-vel = u_init(1) below centerline, u_init(2) above centerline

  !=========================================================
  ! case 4:
  ! not defined

  !=========================================================
  ! case 5:
  ! not defined

  !=========================================================
  ! case 6:
  ! not defined

  !=========================================================
  ! case 7:
  ! not defined

  !=========================================================
  ! case 8:
  ! not defined

  !=========================================================
  ! case 9:
  ! not defined

  !=========================================================
  ! case 10:
  ! not defined

  !=========================================================
  ! case 11:
  ! not defined

  !=========================================================
  ! case 12:

  !=========================================================
  ! case 13:

  !=========================================================
  ! case 14:


contains

  subroutine init_rho_and_umac(mla,rho,umac,dx,time,the_bc_level)

    ! initialize rho_i and umac in the valid region
    ! we first initialize c_i in the valid region
    ! then enforce that sum(c_i)=1 by overwriting the final concentration,
    ! and then use the EOS to compute rho_i

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rho(:)
    type(multifab) , intent(inout) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(in   ) :: time
    type(bc_level) , intent(in   ) :: the_bc_level(:)
 
    ! local variables
    integer                        :: lo(mla%dim), hi(mla%dim)
    integer                        :: dm, ng_c, ng_u, ng_r, i, n, nlevs
    real(kind=dp_t), pointer       :: dp(:,:,:,:)   ! pointer for rho (last dim:nspecies)   
    real(kind=dp_t), pointer       :: up(:,:,:,:)   ! pointers for mac velocities
    real(kind=dp_t), pointer       :: vp(:,:,:,:)
    real(kind=dp_t), pointer       :: wp(:,:,:,:)
    real(kind=dp_t), pointer       :: rp(:,:,:,:)

    type(multifab) :: conc(mla%dim)

    dm = mla%dim
    nlevs = mla%nlevel

    do n=1,nlevs
       call multifab_build(conc(n),mla%la(n),nspecies,rho(n)%ng)
    end do

    ng_u = umac(1,1)%ng
    ng_c = conc(1)%ng
    ng_r = rho(1)%ng

    ! assign values of parameters for the gaussian rho, rhototal
    alpha1 = 0.5d0 
    beta   = 0.1d0 
    delta  = 0.5d0 
    sigma  = (prob_hi(1)-prob_lo(1))/10.0d0  ! variance of gaussian distribution

    ! looping over boxes 
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp => dataptr(conc(n),i)
          rp => dataptr(rho(n),i)
          up => dataptr(umac(n,1),i)
          vp => dataptr(umac(n,2),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          select case (dm)
          case (2)
             call init_rho_and_umac_2d(dp(:,:,1,:),ng_c,rp(:,:,1,:),ng_r,up(:,:,1,1),vp(:,:,1,1),ng_u, &
                                     lo,hi,dx(n,:),time)
          case (3)
             wp => dataptr(umac(n,3),i)
             call init_rho_and_umac_3d(dp(:,:,:,:),ng_c,rp(:,:,:,:),ng_r,up(:,:,:,1),vp(:,:,:,1),wp(:,:,:,1),ng_u, &
                                     lo,hi,dx(n,:),time)
          end select
       end do
    end do

    do n=1,nlevs
       call multifab_destroy(conc(n))
    end do

  end subroutine init_rho_and_umac

  subroutine init_rho_and_umac_2d(c,ng_c,rho,ng_r,u,v,ng_u,lo,hi,dx,time)

    integer          :: lo(2), hi(2), ng_c, ng_u, ng_r
    real(kind=dp_t)  ::   c(lo(1)-ng_c:,lo(2)-ng_c:,:)
    real(kind=dp_t)  :: rho(lo(1)-ng_r:,lo(2)-ng_r:,:)
    real(kind=dp_t)  ::   u(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t)  ::   v(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t)  :: dx(:)
    real(kind=dp_t)  :: time 
 
    ! local varables
    integer          :: i,j,n,seed
    real(kind=dp_t)  :: x,y,w1,w2,rad,rsq,rhot,L(2),sum,r,y1,c_loc,random
    real(kind=dp_t)  :: gradToverT,m_e
 
    L(1:2) = prob_hi(1:2)-prob_lo(1:2) ! Domain length

    seed = n_cells(1)*lo(2) + lo(1)
    call srand(seed)

    select case (abs(prob_type))
    
    case (1)

       !=============================================================
       ! bubble with radius = 1/4 of domain in x
       ! c=rho_init(1,:) inside, c=rho_init(2,:) outside
       ! can be discontinous or smooth depending on smoothing_width
       !=============================================================
 
       u = 0.d0
       v = 0.d0

       rad = L(1)/4.d0

       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+half)*dx(2) - half*(prob_lo(2)+prob_hi(2))
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+half)*dx(1) - half*(prob_lo(1)+prob_hi(1))
       
             r = sqrt(x**2 + y**2)

             if (smoothing_width .eq. 0) then

                ! discontinuous interface
                if (r .lt. rad) then
                   c(i,j,1:nspecies) = rho_init(1,1:nspecies)
                else
                   c(i,j,1:nspecies) = rho_init(2,1:nspecies)
                end if

             else

                ! smooth interface
                c(i,j,1:nspecies-1) = rho_init(1,1:nspecies-1) + &
                     (rho_init(2,1:nspecies-1) - rho_init(1,1:nspecies-1))* &
                     0.5d0*(1.d0 + tanh((r-rad)/(smoothing_width*dx(1))))

             end if
    
          end do
       end do

    case (2) 

       !=========================================================
       ! constant concentration gradient along y
       ! c=rho_init(1,:) on bottom, c=rho_init(2,:) on top
       !=========================================================

       u = 0.d0
       v = 0.d0

       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+half)*dx(2) 
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+half)*dx(1) 

             ! linear gradient in mass fractions
             c(i,j,1:nspecies) = rho_init(1,1:nspecies) + & 
                  (rho_init(2,1:nspecies) - rho_init(1,1:nspecies))*(y-prob_lo(2))/L(2)

          end do
       end do

    case (3)

       !=============================================================
       ! 1 fluid on top of another
       ! c = rho_init(1,:) on bottom; c = rho_init(2,:) on top
       ! smoothing_width > 0 is a tanh smoothed interface where smoothing width is approx the # of grid 
       !   cells and then c = rand*rho_init(1,:) + (1-rand)*rho_init(2,:)
       ! smoothing_width between 0 and -1 is random perturbation where rand = abs(smoothing_width)*rand()
       ! smoothing width of -2 is a sinusoidal perturbation
       ! x-vel = u_init(1) below centerline, u_init(2) above centerline
       !=============================================================

       v = 0.d0

       ! x-velocity = u_init(1) below centerline
       !              u_init(2) above centerline
       do j=lo(2),hi(2)
          y = prob_lo(2) + (j+0.5d0)*dx(2)
          if (y .lt. y1) then
             u(:,j) = u_init(1)
          else
             u(:,j) = u_init(2)
          end if
       end do

       ! middle of domain
       y1 = (prob_lo(2)+prob_hi(2)) / 2.d0

       if (smoothing_width .le. 0.d0 .and. smoothing_width .ge. -1.d0) then

          ! discontinuous version with random perturbation
          do j=lo(2),hi(2)
             y = prob_lo(2) + (j+0.5d0)*dx(2)
             if (y .lt. y1) then
                do n=1,nspecies
                   c(lo(1):hi(1),j,n) = rho_init(1,n)
                end do
             else
                do n=1,nspecies
                   c(lo(1):hi(1),j,n) = rho_init(2,n)
                end do
             end if

             if (j .eq. n_cells(2)/2) then
                do i=lo(1),hi(1)
                   x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                   c_loc = abs(smoothing_width)*rand()
                   do n=1,nspecies
                      c(i,j,n) = c_loc*(rho_init(1,n)) + (1.d0-c_loc)*rho_init(2,n)
                   end do
                end do
             end if

          end do

       else if (smoothing_width .gt. 0.d0) then

          ! smoothed version
          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2)*(dble(j)+0.5d0) - y1
             do n=1,nspecies
                c_loc = rho_init(1,n) + (rho_init(2,n)-rho_init(1,n))*0.5d0*(tanh(y/(smoothing_width*dx(2)))+1.d0)
                c(lo(1):hi(1),j,n) = c_loc
             end do
          end do

       else if (smoothing_width .eq. -2.d0) then

          ! discontinuous version with sinusoidal perturbation
          do j=lo(2),hi(2)
             y = prob_lo(2) + (j+0.5d0)*dx(2)
             if (y .lt. y1) then
                do n=1,nspecies
                   c(lo(1):hi(1),j,n) = rho_init(1,n)
                end do
             else
                do n=1,nspecies
                   c(lo(1):hi(1),j,n) = rho_init(2,n)
                end do
             end if

             if (j .eq. n_cells(2)/2) then
                do i=lo(1),hi(1)
                   x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                   c_loc = 0.5d0*(cos(4.d0*M_PI*x/L(1))+1.d0)
                   do n=1,nspecies
                      c(i,j,n) = c_loc*(rho_init(1,n)) + (1.d0-c_loc)*rho_init(2,n)
                   end do
                end do
             end if

          end do

       else

          call bl_error("init_rho_and_umac_2d: smoothing_width not compatible with prob_type")

       end if

    case (12)

       ! Gaussian bubble
       ! centered in domain
       ! if smoothing_width = 0, this is a discontinuous square in the central 25% of domain

       u = 0.d0
       v = 0.d0

       do j=lo(2),hi(2)
          y = prob_lo(2) + dx(2) * (dble(j)+0.5d0) -  0.5d0*(prob_lo(2)+prob_hi(2))
          do i=lo(1),hi(1)
             x = prob_lo(1) + dx(1) * (dble(i)+0.5d0) - 0.5d0*(prob_lo(1)+prob_hi(1))

             r = sqrt (x**2 + y**2)

             if (smoothing_width .eq. 0.d0) then

                ! initialize c to a square region
                if (i .ge. n_cells(1)/4 .and. i .le. 3*n_cells(1)/4-1 .and. &
                    j .ge. n_cells(2)/4 .and. j .le. 3*n_cells(2)/4-1) then
                   c(i,j,1:nspecies) = rho_init(1,1:nspecies)
                else
                   c(i,j,1:nspecies) = rho_init(2,1:nspecies)
                end if

             else
                ! set c using Gaussian bump
                c(i,j,1) = rho_init(1,1)*exp(-75.d0*r**2)
                c(i,j,2) = rho_init(1,2)*exp(-75.d0*r**2)
             end if

          enddo
       enddo

    case (13)

       ! stratified multispecies due to barodiffusion
       ! assumes the final species is the light solvent

       u = 0.d0
       v = 0.d0

       do n=1,nspecies-1
          m_e = (rhobar(nspecies)/rhobar(n) - 1.d0)*molmass(n)

          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2)*(dble(j)+0.5d0)
             do i=lo(1),hi(1)
                c(i,j,n) = rho_init(1,n)*exp(-m_e*grav(2)*y/(k_B*T_init(1)))
             enddo
          enddo
       enddo

    case (14)

       ! stratified multispecies due to thermodiffusion
       ! assumes the final species is the light solvent
       ! assume ternary for now

       if (nspecies .ne. 3) then
          call bl_error("prob_type=14 requires nspecies=3")
       end if

       gradToverT = (T_init(2)-T_init(1))/(T_init(1)*(prob_hi(2)-prob_lo(2)))

       u = 0.d0
       v = 0.d0

       do j=lo(2),hi(2)
          y = prob_lo(2) + dx(2)*(dble(j)+0.5d0)
          do i=lo(1),hi(1)
             
             c(i,j,1) = rho_init(1,1)*exp((Dtherm(1)-Dtherm(3))*gradToverT*y/Dbar(2))
             c(i,j,2) = rho_init(1,2)*exp((Dtherm(2)-Dtherm(3))*gradToverT*y/Dbar(3))

          enddo
       enddo

    case default

       call bl_error("Desired prob_type not supported in 2D")

    end select

    ! set final c_i such that sum(c_i) = 1 to within roundoff
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          sum = 0
          do n=1,nspecies-1
             sum = sum + c(i,j,n)
          end do
          c(i,j,nspecies) = 1.d0 - sum

       end do
    end do

    ! compute rho using the eos
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          sum = 0.d0
          do n=1,nspecies
             ! sum represents rhoinv
             sum = sum + c(i,j,n)/rhobar(n)
          end do
          rho(i,j,1:nspecies) = c(i,j,1:nspecies)/sum

       end do
    end do

  end subroutine init_rho_and_umac_2d

  subroutine init_rho_and_umac_3d(c,ng_c,rho,ng_r,u,v,w,ng_u,lo,hi,dx,time)
    
    integer          :: lo(3), hi(3), ng_c, ng_u, ng_r
    real(kind=dp_t)  ::   c(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:,:)
    real(kind=dp_t)  :: rho(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:,:)
    real(kind=dp_t)  ::   u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t)  ::   v(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t)  ::   w(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t)  :: dx(:)
    real(kind=dp_t)  :: time 
 
    ! local variables
    integer          :: i,j,k,n,seed
    real(kind=dp_t)  :: x,y,z,rsq,w1,w2,rhot,L(3),sum,random,c_loc,y1,r

    L(1:3) = prob_hi(1:3)-prob_lo(1:3) ! Domain length
    
    seed = n_cells(1)*n_cells(2)*lo(3) + n_cells(1)*lo(2) + lo(1)
    call srand(seed)

    select case (abs(prob_type))

    case (1) 
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
                   c(i,j,k,1:nspecies) = rho_init(1,1:nspecies)
                else
                   c(i,j,k,1:nspecies) = rho_init(2,1:nspecies)
                end if

             end do
          end do
       end do
       !$omp end parallel do

    case (2) 
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

                c(i,j,k,1:nspecies) = rho_init(1,1:nspecies) + &
                     (rho_init(2,1:nspecies) - rho_init(1,1:nspecies))*(y-prob_lo(2))/L(2)

             end do
          end do
       end do
       !$omp end parallel do

    case (3) 

       !=============================================================
       ! 1 fluid on top of another
       ! c(:) = rho_init(1,:) on bottom
       ! c(:) = rho_init(2,:) on top
       !=============================================================

       ! for prob_type=3, smoothing_width is interpreted as:
       ! positive is a tanh smoothed interface where smoothing width is approx the # of grid cells
       ! between 0 and -1 is random perturbation where rand = abs(smoothing_width)*rand() and then
       !  c = rand*rho_init(1,:) + (1-rand)*rho_init(2,:)
       ! -2 is sinusoidal

       u = 0.d0
       v = 0.d0
       w = 0.d0

       ! middle of domain
       y1 = (prob_lo(2)+prob_hi(2)) / 2.d0

       if (smoothing_width .le. 0.d0 .and. smoothing_width .ge. -1.d0) then

          ! discontinuous version with random perturbation
          do j=lo(2),hi(2)

             y = prob_lo(2) + (j+0.5d0)*dx(2)

             if (y .lt. y1) then
                do n=1,nspecies
                   c(lo(1):hi(1),j,lo(3):hi(3),n) = rho_init(1,n)
                end do
             else
                do n=1,nspecies
                   c(lo(1):hi(1),j,lo(3):hi(3),n) = rho_init(2,n)
                end do
             end if

             if (j .eq. n_cells(2)/2) then
                do k=lo(3),hi(3)
                   do i=lo(1),hi(1)
                      c_loc = abs(smoothing_width)*rand()
                      do n=1,nspecies
                         c(i,j,k,n) = c_loc*(rho_init(1,n)) + (1.d0-c_loc)*rho_init(2,n)
                      end do
                   end do
                end do
             end if

          end do

       else if (smoothing_width .gt. 0.d0) then

          ! smoothed version
          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2)*(dble(j)+0.5d0) - y1
             do n=1,nspecies
                c_loc = rho_init(1,n) + (rho_init(2,n)-rho_init(1,n))*0.5d0*(tanh(y/(smoothing_width*dx(2)))+1.d0)
                c(lo(1):hi(1),j,lo(3):hi(3),n) = c_loc
             end do
          end do

       else if (smoothing_width .eq. -2.d0) then

          ! discontinuous version with sinusoidal perturbation
          do j=lo(2),hi(2)

             y = prob_lo(2) + (j+0.5d0)*dx(2)

             if (y .lt. y1) then
                do n=1,nspecies
                   c(lo(1):hi(1),j,lo(3):hi(3),n) = rho_init(1,n)
                end do
             else
                do n=1,nspecies
                   c(lo(1):hi(1),j,lo(3):hi(3),n) = rho_init(2,n)
                end do
             end if

             if (j .eq. n_cells(2)/2) then
                do k=lo(3),hi(3)
                   do i=lo(1),hi(1)
                      x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                      c_loc = 0.5d0*(cos(4.d0*M_PI*x/L(1))+1.d0)
                      do n=1,nspecies
                         c(i,j,k,n) = c_loc*(rho_init(1,n)) + (1.d0-c_loc)*rho_init(2,n)
                      end do
                   end do
                end do
             end if

          end do

       else

          call bl_error("init_rho_and_umac_3d: smoothing_width not compatible with prob_type")

       end if

    case (4)
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

                c(i,j,k,1) = 1.0d0/(4.0d0*M_PI*Dbar(1)*time)**1.5d0*dexp(-rsq/&
                     (4.0d0*Dbar(1)*time) - time*beta)*rhot
                c(i,j,k,2) = rhot - c(i,j,k,1) 

             end do
          end do
       end do
       !$omp end parallel do

    case (5)
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
                c(i,j,k,1) = rhot*w1
                c(i,j,k,2) = rhot*w2
                c(i,j,k,3) = rhot-c(i,j,k,1)-c(i,j,k,2)

                if(c(i,j,k,1).lt.0.d0 .or. c(i,j,k,2).lt.0.d0 .or. c(i,j,k,3).lt.0.d0) then 
                   write(*,*), "rho1 / rho2 / rho3 is negative: STOP"
                   write(*,*), i, j, " w1=", w1, " w2=", w2, " rho1=",c(i,j,k,1)," rho2=",&
                        c(i,j,k,2), " rho3=",c(i,j,k,3), " rhot=",rhot
                end if

             end do
          end do
       end do
       !$omp end parallel do

    case (7)

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

                c(i,j,k,1:nspecies-1) = rho_init(1,1:nspecies-1) + &
                     0.5d0*(rho_init(2,1:nspecies-1) - rho_init(1,1:nspecies-1))* &
                     (1.d0 + tanh((r-15.d0)/2.d0))

             end do
          end do
       end do

    case (8)

       !=============================================================
       ! 4-species, 4-stripes
       !=============================================================

       u = 0.d0
       v = 0.d0
       w = 0.d0

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+half)*dx(2) 
             do i=lo(1),hi(1)
                if (y .le. 0.75d0*prob_lo(2)+0.25d0*prob_hi(2)) then
                   c(i,j,k,1) = rho_init(1,1)
                   c(i,j,k,2) = 0.1d0*rhobar(2)
                   c(i,j,k,3) = 0.1d0*rhobar(3)
                   c(i,j,k,4) = 0.1d0*rhobar(4)
                else if (y .le. 0.5d0*prob_lo(2)+0.5d0*prob_hi(2)) then
                   c(i,j,k,1) = 0.1d0*rhobar(1)
                   c(i,j,k,2) = rho_init(1,2)
                   c(i,j,k,3) = 0.1d0*rhobar(3)
                   c(i,j,k,4) = 0.1d0*rhobar(4)
                else if (y .le. 0.25d0*prob_lo(2)+0.75d0*prob_hi(2)) then
                   c(i,j,k,1) = 0.1d0*rhobar(1)
                   c(i,j,k,2) = 0.1d0*rhobar(2)
                   c(i,j,k,3) = rho_init(1,3)
                   c(i,j,k,4) = 0.1d0*rhobar(4)
                else
                   c(i,j,k,1) = 0.1d0*rhobar(1)
                   c(i,j,k,2) = 0.1d0*rhobar(2)
                   c(i,j,k,3) = 0.1d0*rhobar(3)
                   c(i,j,k,4) = rho_init(1,4)
                end if
             end do
          end do
       end do

    case (9)

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

             c_loc = rho_init(1,1) + (rho_init(1,2)-rho_init(1,1))*0.5d0*(tanh(y/(smoothing_width*dx(2)))+1.d0)
             c(lo(1):hi(1),j,lo(3):hi(3),1) = c_loc

          end do

       else

          ! discontinuous version
          do j=lo(2),hi(2)
             y = prob_lo(2) + (j+0.5d0)*dx(2)
             if (y .lt. y1) then
                c(lo(1):hi(1),j,lo(3):hi(3),1) = rho_init(1,1)
             else
                c(lo(1):hi(1),j,lo(3):hi(3),1) = rho_init(1,2)
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

       ! c1 = rho_init(1,1) in lower half of domain (in y)
       ! c1 = rho_init(2,1) in upper half
       ! random perturbation below centerline

       do j=lo(2),hi(2)
          y = prob_lo(2) + (j+0.5d0)*dx(2)

          if (y .lt. y1) then
             c_loc = rho_init(1,1)
          else
             c_loc = rho_init(2,1)
          end if

          c(lo(1):hi(1),j,lo(3):hi(3),1) = c_loc
          c(lo(1):hi(1),j,lo(3):hi(3),2) = 1.d0 - c_loc

          ! add random perturbation below centerline
          if (j .eq. n_cells(2)/2-1) then
             do k=lo(3),hi(3)
                do i=lo(1),hi(1)
                   random = rand()
                   c_loc = random*rho_init(1,1) + (1.d0-random)*rho_init(2,1)
                   c(i,j,k,1) = c_loc
                   c(i,j,k,2) = 1.d0 - c_loc
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
       ! c(:) = rho_init(1,:) on bottom
       ! c(:) = rho_init(2,:) on top
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
                   c(i,j,k,1:nspecies) = rho_init(1,1:nspecies)
                end do
             end do
          else
             do k=lo(3),hi(3)
                do i=lo(1),hi(1)
                   c(i,j,k,1:nspecies) = rho_init(2,1:nspecies)
                end do
             end do
          end if

       end do

    case (12)

       ! Gaussian bubble
       ! centered in domain
       ! if smoothing_width = 0, this is a discontinuous square in the central 25% of domain

       u = 0.d0
       v = 0.d0
       w = 0.d0

       do k=lo(3),hi(3)
          z = prob_lo(3) + dx(3) * (dble(k)+0.5d0) -  0.5d0*(prob_lo(3)+prob_hi(3))
          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2) * (dble(j)+0.5d0) -  0.5d0*(prob_lo(2)+prob_hi(2))
             do i=lo(1),hi(1)
                x = prob_lo(1) + dx(1) * (dble(i)+0.5d0) - 0.5d0*(prob_lo(1)+prob_hi(1))

                r = sqrt (x**2 + y**2 + z**2)

                if (smoothing_width .eq. 0.d0) then

                   ! initialize c to a square region
                   if (i .ge. n_cells(1)/4 .and. i .le. 3*n_cells(1)/4-1 .and. &
                       j .ge. n_cells(2)/4 .and. j .le. 3*n_cells(2)/4-1 .and. &
                       k .ge. n_cells(3)/4 .and. k .le. 3*n_cells(3)/4-1) then
                      c(i,j,k,1:nspecies) = rho_init(1,1:nspecies)
                   else
                      c(i,j,k,1:nspecies) = rho_init(2,1:nspecies)
                   end if

                else
                   ! set c using Gaussian bump
                   c(i,j,k,1) = rho_init(1,1)*exp(-75.d0*r**2)
                   c(i,j,k,2) = rho_init(1,2)*exp(-75.d0*r**2)
                end if

             enddo
          enddo
       enddo

    case default

       call bl_error("Desired prob_type not supported in 3D")

    end select

    ! set final c_i such that sum(c_i) = 1
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             sum = 0
             do n=1,nspecies-1
                sum = sum + c(i,j,k,n)
             end do
             c(i,j,k,nspecies) = 1.d0 - sum

          end do
       end do
    end do

    ! compute rho using the eos
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             sum = 0.d0
             do n=1,nspecies
                ! sum represents rhoinv
                sum = sum + c(i,j,k,n)/rhobar(n)
             end do
             rho(i,j,k,1:nspecies) = c(i,j,k,1:nspecies)/sum

          end do
       end do
    end do

  end subroutine init_rho_and_umac_3d

end module init_lowmach_module
