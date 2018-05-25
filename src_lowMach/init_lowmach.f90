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
  use BoxLibRNGs
  use bl_rng_module
  use probin_common_module, only: prob_lo, prob_hi, prob_type, k_B, grav, &
                                  molmass, rhobar, rho0, smoothing_width, u_init, n_cells, &
                                  use_bl_rng, nspecies, algorithm_type, initial_variance_mass
  use probin_multispecies_module, only: alpha1, beta, delta, sigma, Dbar, Dtherm, &
                                        c_init, T_init, temp_type, sigma, is_ideal_mixture
  use probin_charged_module, only: charge_per_mass, use_charged_fluid, electroneutral
 
  implicit none

  private

  public :: init_rho_and_umac, & ! used in low Mach code; initialize c first then convert to rho
            init_temp

  ! IMPORTANT: In the diffusion only code (init_rho), c_init specifies initial values for DENSITY
  ! In the low-Mach code (init_rho_and_umac), c_init specifies initial MASS FRACTIONS
  ! (should sum to unity!... but we overwrite the final concentration so sum(c_i)=1 before computing rho)
  ! The density follows from the EOS in the LM case so it cannot be specified
  ! Same applies to boundary conditions

  ! prob_types codes for init_lowmach:

  !=============================================================
  ! case 1:
  ! bubble with radius = 1/4 of domain in x
  ! c=c_init(1,:) inside, c=c_init(2,:) outside
  ! can be discontinous or smooth depending on smoothing_width

  !=========================================================
  ! case 2:
  ! constant concentration gradient along y
  ! c=c_init(1,:) on bottom, c=c_init(2,:) on top

  !=========================================================
  ! case 3:
  ! one fluid on top of another
  ! c = c_init(1,:) on bottom; c = c_init(2,:) on top
  ! smoothing_width > 0 is a tanh smoothed interface where smoothing width is approx the # of grid 
  !   cells and then c = rand*c_init(1,:) + (1-rand)*c_init(2,:)
  ! smoothing_width between 0 and -1 is random perturbation where rand = abs(smoothing_width)*rand()
  ! smoothing width of -2 is a sinusoidal perturbation
  ! x-vel = u_init(1) below centerline, u_init(2) above centerline

  !=========================================================
  ! case 4:
  ! not defined
  ! external_source has an analytic solution, not sure what it was for
  ! I think it's a remnant from test_diffusion

  !=========================================================
  ! case 5:
  ! not defined
  ! external_source has an analytic solution, not sure what it was for
  ! I think it's a remnant from test_diffusion

  !=========================================================
  ! case 6:
  ! Two Gaussian bubbles with different centers
  ! c(1) peak is c_init(1,1) 
  ! c(2) peak is c_init(2,2) 

  !=========================================================
  ! case 7:
  ! 1D sin wave of c1

  !=========================================================
  ! case 8:
  ! not defined

  !=========================================================
  ! case 9:
  ! not defined
  ! mixture_properties_mass_local computes Dbar's using water/glycerol
  ! compute_eta uses water/glycerol

  !=========================================================
  ! case 10:
  ! not defined

  !=========================================================
  ! case 11:
  ! Discontinuous square in the central 25% of domain
  ! c=c_init(1,:) inside; c=c_init(2,:) outside

  !=========================================================
  ! case 12:
  ! Gaussian bubble centered in domain
  ! c=c_init(1,:) inside; c=c_init(2,:) outside
  ! lo- and hi-y walls move with prescribed velocity,
  ! see inhomogeneous_bc_val.f90
  ! compute_eta uses linear profile in rho if prob_type < 0

  !=========================================================
  ! case 13:
  ! stratified multispecies due to barodiffusion
  ! approximate analytical steady solution
  ! assumes the final species is the dominant component

  !=========================================================
  ! case 14:
  ! stratified multispecies due to thermodiffusion
  ! approximate analytical steady solution
  ! assumes nspecies=3
  ! assumes the final species is the dominant component

  !=========================================================
  ! case 15:
  ! Discontinuous band in central 50% of domain
  ! c=c_init(1,:) inside; c=c_init(2,:) outside
  ! if prob_type=-15, add another tanh along other dimension for last two species (two stripes crossing)

  !=========================================================
  ! case 16:
  ! Same as case 15 but now the central stripe is 1/3 of the domain

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

    type(bl_prof_timer), save :: bpt

    ! currently, adding initial mass fluctuations is only available for ideal mixture
    if (abs(initial_variance_mass) .gt. 0.d0) then
      if (.not. is_ideal_mixture) then
        call bl_error("currently, nonzero initial_variance_mass is supported only for ideal mixture")
      end if
    end if

    call build(bpt, "init_rho_and_umac")

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
             call init_rho_and_umac_2d(dp(:,:,1,:),ng_c,rp(:,:,1,:),ng_r, &
                                       up(:,:,1,1),vp(:,:,1,1),ng_u, &
                                       lo,hi,dx(n,:),time)
          case (3)
             wp => dataptr(umac(n,3),i)
             call init_rho_and_umac_3d(dp(:,:,:,:),ng_c,rp(:,:,:,:),ng_r, &
                                       up(:,:,:,1),vp(:,:,:,1),wp(:,:,:,1),ng_u, &
                                       lo,hi,dx(n,:),time)
          end select
       end do
    end do

    do n=1,nlevs
       call multifab_destroy(conc(n))
    end do

    call destroy(bpt)

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
    integer          :: i,j,n
    real(kind=dp_t)  :: x,y,rad,L(2),sum,r,r1,r2,y1,y2,c_loc,x1,x2,coeff
    real(kind=dp_t)  :: gradToverT,m_e
 
    real(kind=dp_t)  :: random

    real(kind=dp_t)  :: rho_total

    L(1:2) = prob_hi(1:2)-prob_lo(1:2) ! Domain length

    select case (abs(prob_type))
    
    case (1)

       !=============================================================
       ! bubble with radius = 1/4 of domain in x
       ! c=c_init(1,:) inside, c=c_init(2,:) outside
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
                   c(i,j,1:nspecies) = c_init(1,1:nspecies)
                else
                   c(i,j,1:nspecies) = c_init(2,1:nspecies)
                end if

             else

                ! smooth interface
                c(i,j,1:nspecies-1) = c_init(1,1:nspecies-1) + &
                     (c_init(2,1:nspecies-1) - c_init(1,1:nspecies-1))* &
                     0.5d0*(1.d0 + tanh((r-rad)/(smoothing_width*dx(1))))

             end if
    
          end do
       end do

    case (2) 

       !=========================================================
       ! constant concentration gradient along y
       ! c=c_init(1,:) on bottom, c=c_init(2,:) on top
       !=========================================================

       u = 0.d0
       v = 0.d0

       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+half)*dx(2) 
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+half)*dx(1) 

             ! linear gradient in mass fractions
             c(i,j,1:nspecies) = c_init(1,1:nspecies) + & 
                  (c_init(2,1:nspecies) - c_init(1,1:nspecies))*(y-prob_lo(2))/L(2)

          end do
       end do

    case (3)

       !=============================================================
       ! 1 fluid on top of another
       ! c = c_init(1,:) on bottom; c = c_init(2,:) on top
       ! smoothing_width > 0 is a tanh smoothed interface where smoothing width is approx the # of grid 
       !   cells and then c = rand*c_init(1,:) + (1-rand)*c_init(2,:)
       ! smoothing_width between 0 and -1 is random perturbation where rand = abs(smoothing_width)*rand()
       ! smoothing width of -2 is a sinusoidal perturbation
       ! x-vel = u_init(1) below centerline, u_init(2) above centerline
       !=============================================================

       u = 0.d0
       v = 0.d0

       ! middle of domain
       y1 = (prob_lo(2)+prob_hi(2)) / 2.d0

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

       if (smoothing_width .le. 0.d0 .and. smoothing_width .ge. -1.d0) then
          ! discontinuous version with random perturbation
          do j=lo(2),hi(2)
             y = prob_lo(2) + (j+0.5d0)*dx(2)
             if (y .lt. y1) then
                do n=1,nspecies
                   c(lo(1):hi(1),j,n) = c_init(1,n)
                end do
             else
                do n=1,nspecies
                   c(lo(1):hi(1),j,n) = c_init(2,n)
                end do
             end if

             if (j .eq. n_cells(2)/2) then
                do i=lo(1),hi(1)
                   x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                   if (use_bl_rng) then
                      call UniformRNG(random,engine=rng_eng_init_mass%p)
                   else
                      call UniformRNG(random)
                   end if
                   c_loc = abs(smoothing_width)*random
                   do n=1,nspecies
                      c(i,j,n) = c_loc*(c_init(1,n)) + (1.d0-c_loc)*c_init(2,n)
                   end do
                end do
             end if

          end do

       else if (smoothing_width .gt. 0.d0) then

          ! smoothed version
          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2)*(dble(j)+0.5d0) - y1
             do n=1,nspecies                
                c_loc = c_init(1,n) + (c_init(2,n)-c_init(1,n))*0.5d0*(tanh(y/(smoothing_width*dx(2)))+1.d0)
                c(lo(1):hi(1),j,n) = c_loc
             end do
          end do
          

       else if (smoothing_width .eq. -2.d0) then

          ! discontinuous version with sinusoidal perturbation
          do j=lo(2),hi(2)
             y = prob_lo(2) + (j+0.5d0)*dx(2)
             if (y .lt. y1) then
                do n=1,nspecies
                   c(lo(1):hi(1),j,n) = c_init(1,n)
                end do
             else
                do n=1,nspecies
                   c(lo(1):hi(1),j,n) = c_init(2,n)
                end do
             end if

             if (j .eq. n_cells(2)/2) then
                do i=lo(1),hi(1)
                   x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                   c_loc = 0.5d0*(cos(4.d0*M_PI*x/L(1))+1.d0)
                   do n=1,nspecies
                      c(i,j,n) = c_loc*(c_init(1,n)) + (1.d0-c_loc)*c_init(2,n)
                   end do
                end do
             end if

          end do

       else

          call bl_error("init_rho_and_umac_2d: smoothing_width not compatible with prob_type")

       end if

    case (6)

       !=============================================================
       ! Two Gaussian bubbles with different centers
       ! c(1) peak is c_init(1,1) 
       ! c(2) peak is c_init(2,2) 
       !=============================================================

       u = 0.d0
       v = 0.d0

       do j=lo(2),hi(2)
          y = prob_lo(2) + dx(2) * (dble(j)+0.5d0)
          do i=lo(1),hi(1)
             x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)

             r1 = sqrt ((x-(0.5d0*prob_lo(1)+0.5d0*prob_hi(1)))**2 + (y-(0.6d0*prob_lo(2)+0.4d0*prob_hi(2)))**2)
             r2 = sqrt ((x-(0.5d0*prob_lo(1)+0.5d0*prob_hi(1)))**2 + (y-(0.4d0*prob_lo(2)+0.6d0*prob_hi(2)))**2)
             

             ! set c using Gaussian bump
             c(i,j,1) = c_init(1,1)*exp(-75.d0*r1**2)
             c(i,j,2) = c_init(2,2)*exp(-75.d0*r2**2)

          enddo
       enddo

    case (7)

       !=============================================================
       ! 1D sin wave of c1
       !=============================================================

       u = 0.d0
       v = 0.d0

       ! note: c(:,:,3) will be computed below to enforce sum(c)=1
       c(:,:,2) = 0.05d0

       do j=lo(2),hi(2)
          y = prob_lo(2) + dx(2) * (dble(j)+0.5d0)
          do i=lo(1),hi(1)
             x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)

             c(i,j,1) = 0.05d0 + 0.0005d0*sin(2.d0*M_PI*x/L(1))

          enddo
       enddo

    case (11)

       !=============================================================
       ! Discontinuous square in the central 25% of domain
       ! c=c_init(1,:) inside; c=c_init(2,:) outside
       !=============================================================

       u = 0.d0
       v = 0.d0

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ! initialize c to a square region
             if (i .ge. n_cells(1)/4 .and. i .le. 3*n_cells(1)/4-1 .and. &
                 j .ge. n_cells(2)/4 .and. j .le. 3*n_cells(2)/4-1) then
                c(i,j,1:nspecies) = c_init(1,1:nspecies)
             else
                c(i,j,1:nspecies) = c_init(2,1:nspecies)
             end if

          enddo
       enddo

    case (12)

       !=============================================================
       ! Gaussian bubble centered in domain
       ! c=c_init(1,:) inside; c=c_init(2,:) outside
       ! lo- and hi-y walls move with prescribed velocity,
       ! see inhomogeneous_bc_val.f90
       !=============================================================

       u = 0.d0
       v = 0.d0

       do j=lo(2),hi(2)
          y = prob_lo(2) + dx(2) * (dble(j)+0.5d0) - 0.5d0*(prob_lo(2)+prob_hi(2))
          do i=lo(1),hi(1)
             x = prob_lo(1) + dx(1) * (dble(i)+0.5d0) - 0.5d0*(prob_lo(1)+prob_hi(1))

             r = sqrt (x**2 + y**2)

             ! set c using Gaussian bump
             c(i,j,1:nspecies-1) = c_init(1,1:nspecies-1)*exp(-75.d0*r**2)

          enddo
       enddo

    case (13)

       !=============================================================
       ! stratified multispecies due to barodiffusion
       ! approximate analytical steady solution
       ! assumes the final species is the dominant component
       !=============================================================

       u = 0.d0
       v = 0.d0

       do n=1,nspecies-1
          m_e = (rhobar(nspecies)/rhobar(n) - 1.d0)*molmass(n)

          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2)*(dble(j)+0.5d0)
             do i=lo(1),hi(1)
                c(i,j,n) = c_init(1,n)*exp(-m_e*grav(2)*y/(k_B*T_init(1)))
             enddo
          enddo
       enddo

    case (14)

       !=============================================================
       ! stratified multispecies due to thermodiffusion
       ! approximate analytical steady solution
       ! assumes nspecies=3
       ! assumes the final species is the dominant component
       !=============================================================

       if (nspecies .ne. 3) then
          call bl_error("prob_type=14 requires nspecies=3")
       end if

       gradToverT = (T_init(2)-T_init(1))/(T_init(1)*(prob_hi(2)-prob_lo(2)))

       u = 0.d0
       v = 0.d0

       do j=lo(2),hi(2)
          y = prob_lo(2) + dx(2)*(dble(j)+0.5d0)
          do i=lo(1),hi(1)
             
             c(i,j,1) = c_init(1,1)*exp((Dtherm(1)-Dtherm(3))*gradToverT*y/Dbar(2))
             c(i,j,2) = c_init(1,2)*exp((Dtherm(2)-Dtherm(3))*gradToverT*y/Dbar(3))

          enddo
       enddo

    case (15,16)

       !=========================================================
       ! Discontinuous band in central 50% of domain
       ! c=c_init(1,:) inside; c=c_init(2,:) outside
       !=============================================================

       u = 0.d0
       v = 0.d0

       if(prob_type==15) then
          ! first quarter of domain
          y1 = (3*prob_lo(2) + prob_hi(2)) / 4.d0
          x1 = (3*prob_lo(1) + prob_hi(1)) / 4.d0

          ! last quarter of domain
          y2 = (prob_lo(2) + 3*prob_hi(2)) / 4.d0
          x2 = (prob_lo(1) + 3*prob_hi(1)) / 4.d0       
       else
          ! first third of domain
          y1 = (2*prob_lo(2) + prob_hi(2)) / 3.d0
          x1 = (2*prob_lo(1) + prob_hi(1)) / 3.d0
       
          ! last third of domain
          y2 = (prob_lo(2) + 2*prob_hi(2)) / 3.d0
          x2 = (prob_lo(1) + 2*prob_hi(1)) / 3.d0
       end if

       ! x-velocity = u_init(1) below centerline
       !              u_init(2) above centerline
       do j=lo(2),hi(2)
          y = prob_lo(2) + (j+0.5d0)*dx(2)
          if (y .lt. y2 .and. y .gt. y1) then
             u(:,j) = u_init(1)
          else
             u(:,j) = u_init(2)
          end if
       end do

       if (smoothing_width .le. 0.d0 .and. smoothing_width .ge. -1.d0) then
          ! discontinuous version with random perturbation
          do j=lo(2),hi(2)
             y = prob_lo(2) + (j+0.5d0)*dx(2)
             if (y .lt. y2 .and. y .gt. y1) then
                do n=1,nspecies
                   c(lo(1):hi(1),j,n) = c_init(1,n)
                end do
             else
                do n=1,nspecies
                   c(lo(1):hi(1),j,n) = c_init(2,n)
                end do
             end if

             if (j .eq. n_cells(2)/4) then
                do i=lo(1),hi(1)
                   x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                   if (use_bl_rng) then
                      call UniformRNG(random,engine=rng_eng_init_mass%p)
                   else
                      call UniformRNG(random)
                   end if
                   c_loc = abs(smoothing_width)*random
                   do n=1,nspecies
                      c(i,j,n) = c_loc*(c_init(1,n)) + (1.d0-c_loc)*c_init(2,n)
                   end do
                end do
             end if

             if (j .eq. 3*n_cells(2)/4) then
                do i=lo(1),hi(1)
                   x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                   if (use_bl_rng) then
                      call UniformRNG(random,engine=rng_eng_init_mass%p)
                   else
                      call UniformRNG(random)
                   end if
                   c_loc = abs(smoothing_width)*random
                   do n=1,nspecies
                      c(i,j,n) = c_loc*(c_init(1,n)) + (1.d0-c_loc)*c_init(2,n)
                   end do
                end do
             end if

          end do

       else if (smoothing_width .gt. 0.d0) then

          ! smoothed version
          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2)*(dble(j)+0.5d0) - y1
             do i=lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+0.5d0)*dx(1) - x1
                do n=1,nspecies

                   ! tanh smoothing in y
                   coeff=0.5d0*(tanh(y/(smoothing_width*dx(2)))+1.d0)*0.5d0*(tanh((-y+y2-y1)/(smoothing_width*dx(2)))+1.d0)

                   ! Donev: prob_type = -15: a special case for doing ternary diffusion NaCl + KCl
                   ! Here the last two species have a tanh profile in both x and y (species are Na+,Cl-,K+,water)
                   ! Aadd a tanh smoothing for central 50% of domain in x for second-to-last species
                   if( (prob_type==-15) .and. (n==nspecies-1) ) then 
                      coeff=0.5d0*(tanh(x/(smoothing_width*dx(1)))+1.d0)*0.5d0*(tanh((-x+x2-x1)/(smoothing_width*dx(1)))+1.d0)*coeff
                   end if  
 
                   ! smooth between c_init(1,:) and c_init(2,:)
                   c_loc = c_init(2,n) + (c_init(1,n)-c_init(2,n))*coeff
                   c(i,j,n) = c_loc

                   ! for 4-species test, need to add Cl to central square to balance the K
                   if ( (prob_type==-15) .and. (nspecies .eq. 4) .and. (n .eq. nspecies-1) ) then
                      c(i,j,2) = c(i,j,2) - charge_per_mass(3)/charge_per_mass(2)*c_loc
                   end if

                end do
             end do
          end do

       else if (smoothing_width .eq. -2.d0) then

          ! discontinuous version with sinusoidal perturbation
          do j=lo(2),hi(2)
             y = prob_lo(2) + (j+0.5d0)*dx(2)
             if (y .lt. y2 .and. y .gt. y1) then
                do n=1,nspecies
                   c(lo(1):hi(1),j,n) = c_init(1,n)
                end do
             else
                do n=1,nspecies
                   c(lo(1):hi(1),j,n) = c_init(2,n)
                end do
             end if

             if (j .eq. 3*n_cells(2)/4) then
                do i=lo(1),hi(1)
                   x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                   c_loc = 0.5d0*(cos(4.d0*M_PI*x/L(1))+1.d0)
                   do n=1,nspecies
                      c(i,j,n) = c_loc*(c_init(1,n)) + (1.d0-c_loc)*c_init(2,n)
                   end do
                end do
             end if

             if (j .eq. n_cells(2)/4) then
                do i=lo(1),hi(1)
                   x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                   c_loc = 0.5d0*(cos(4.d0*M_PI*x/L(1))+1.d0)
                   do n=1,nspecies
                      c(i,j,n) = c_loc*(c_init(2,n)) + (1.d0-c_loc)*c_init(1,n)
                   end do
                end do
             end if

          end do

       else

          call bl_error("init_rho_and_umac_2d: smoothing_width not compatible with prob_type")

       end if


    case default

       call bl_error("Desired prob_type not supported in 2D")

    end select

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          ! set final c_i such that sum(c_i) = 1 to within roundoff
          sum = 0.d0
          do n=1,nspecies-1
             sum = sum + c(i,j,n)
          end do
          c(i,j,nspecies) = 1.d0 - sum

          ! calculate rho_total from eos
          if (algorithm_type .eq. 6) then
             rho_total = rho0
          else
             sum = 0.d0
             do n=1,nspecies
                ! sum represents rhoinv
                sum = sum + c(i,j,n)/rhobar(n)
             end do
             rho_total = 1.d0/sum
          end if

          ! add mass fluctuations
          if (abs(initial_variance_mass) .gt. 0.d0) then
             call add_mass_fluctuations(c(i,j,1:nspecies),dx,abs(initial_variance_mass),rho_total)

             ! calculate rho_total from eos including fluctuations
             if (algorithm_type .ne. 6) then
                sum = 0.d0
                do n=1,nspecies
                   ! sum represents rhoinv
                   sum = sum + c(i,j,n)/rhobar(n)
                end do
                rho_total = 1.d0/sum
             end if

          end if

          ! calculate rho_i
          rho(i,j,1:nspecies) = rho_total*c(i,j,1:nspecies)

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
    integer          :: i,j,k,n
    real(kind=dp_t)  :: x,y,z,rad,L(3),sum,c_loc,y1,r,r1,r2,m_e,gradToverT

    real(kind=dp_t) :: random

    real(kind=dp_t)  :: rho_total

    L(1:3) = prob_hi(1:3)-prob_lo(1:3) ! Domain length
    
    select case (abs(prob_type))

    case (1) 

       !=============================================================
       ! bubble with radius = 1/4 of domain in x
       ! c=c_init(1,:) inside, c=c_init(2,:) outside
       ! can be discontinous or smooth depending on smoothing_width
       !=============================================================

       u = 0.d0
       v = 0.d0
       w = 0.d0

       rad = L(1)/4.d0

       !$omp parallel do private(i,j,k,x,y,z,r)
       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+half)*dx(3) - half*(prob_lo(3)+prob_hi(3))
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+half)*dx(2) - half*(prob_lo(2)+prob_hi(2))
             do i=lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+half)*dx(1) - half*(prob_lo(1)+prob_hi(1))

                r = sqrt(x**2 + y**2 + z**2)

                if (smoothing_width .eq. 0) then

                   ! discontinuous interface
                   if (r .lt. rad) then
                      c(i,j,k,1:nspecies) = c_init(1,1:nspecies)
                   else
                      c(i,j,k,1:nspecies) = c_init(2,1:nspecies)
                   end if

                else

                   ! smooth interface
                   c(i,j,k,1:nspecies-1) = c_init(1,1:nspecies-1) + &
                        (c_init(2,1:nspecies-1) - c_init(1,1:nspecies-1))* &
                        0.5d0*(1.d0 + tanh((r-rad)/(smoothing_width*dx(1))))

                end if

             end do
          end do
       end do
       !$omp end parallel do

    case (2) 

       !=========================================================
       ! constant concentration gradient along y
       ! c=c_init(1,:) on bottom, c=c_init(2,:) on top
       !=========================================================

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

                c(i,j,k,1:nspecies) = c_init(1,1:nspecies) + &
                     (c_init(2,1:nspecies) - c_init(1,1:nspecies))*(y-prob_lo(2))/L(2)

             end do
          end do
       end do
       !$omp end parallel do

    case (3) 

       !=============================================================
       ! 1 fluid on top of another
       ! c = c_init(1,:) on bottom; c = c_init(2,:) on top
       ! smoothing_width > 0 is a tanh smoothed interface where smoothing width is approx the # of grid 
       !   cells and then c = rand*c_init(1,:) + (1-rand)*c_init(2,:)
       ! smoothing_width between 0 and -1 is random perturbation where rand = abs(smoothing_width)*rand()
       ! smoothing width of -2 is a sinusoidal perturbation
       ! x-vel = u_init(1) below centerline, u_init(2) above centerline
       !=============================================================

       u = 0.d0
       v = 0.d0
       w = 0.d0

       ! middle of domain
       y1 = (prob_lo(2)+prob_hi(2)) / 2.d0

       ! x-velocity = u_init(1) below centerline
       !              u_init(2) above centerline
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             y = prob_lo(2) + (j+0.5d0)*dx(2)
             if (y .lt. y1) then
                u(:,j,k) = u_init(1)
             else
                u(:,j,k) = u_init(2)
             end if
          end do
       end do

       if (smoothing_width .le. 0.d0 .and. smoothing_width .ge. -1.d0) then

          ! discontinuous version with random perturbation
          do j=lo(2),hi(2)
             y = prob_lo(2) + (j+0.5d0)*dx(2)
             if (y .lt. y1) then
                do n=1,nspecies
                   c(lo(1):hi(1),j,lo(3):hi(3),n) = c_init(1,n)
                end do
             else
                do n=1,nspecies
                   c(lo(1):hi(1),j,lo(3):hi(3),n) = c_init(2,n)
                end do
             end if

             if (j .eq. n_cells(2)/2) then
                do k=lo(3),hi(3)
                   do i=lo(1),hi(1)
                      if (use_bl_rng) then
                         call UniformRNG(random,engine=rng_eng_init_mass%p)
                      else
                         call UniformRNG(random)
                      end if
                      c_loc = abs(smoothing_width)*random
                      do n=1,nspecies
                         c(i,j,k,n) = c_loc*(c_init(1,n)) + (1.d0-c_loc)*c_init(2,n)
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
                c_loc = c_init(1,n) + (c_init(2,n)-c_init(1,n))*0.5d0*(tanh(y/(smoothing_width*dx(2)))+1.d0)
                c(lo(1):hi(1),j,lo(3):hi(3),n) = c_loc
             end do
          end do

       else if (smoothing_width .eq. -2.d0) then

          ! discontinuous version with sinusoidal perturbation
          do j=lo(2),hi(2)
             y = prob_lo(2) + (j+0.5d0)*dx(2)
             if (y .lt. y1) then
                do n=1,nspecies
                   c(lo(1):hi(1),j,lo(3):hi(3),n) = c_init(1,n)
                end do
             else
                do n=1,nspecies
                   c(lo(1):hi(1),j,lo(3):hi(3),n) = c_init(2,n)
                end do
             end if
             if (j .eq. n_cells(2)/2) then
                do k=lo(3),hi(3)
                   do i=lo(1),hi(1)
                      x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                      c_loc = 0.5d0*(cos(4.d0*M_PI*x/L(1))+1.d0)
                      do n=1,nspecies
                         c(i,j,k,n) = c_loc*(c_init(1,n)) + (1.d0-c_loc)*c_init(2,n)
                      end do
                   end do
                end do
             end if

          end do

       else

          call bl_error("init_rho_and_umac_3d: smoothing_width not compatible with prob_type")

       end if

    case (6)

       !=============================================================
       ! Two Gaussian bubbles with different centers
       ! c(1) peak is c_init(1,1)
       ! c(2) peak is c_init(2,2)
       !=============================================================

       u = 0.d0
       v = 0.d0
       w = 0.d0
       
       do k=lo(3),hi(3)
          z = prob_lo(3) + dx(3) * (dble(k)+0.5d0)
          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2) * (dble(j)+0.5d0)
             do i=lo(1),hi(1)
                x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)
                
                r1 = sqrt (   (x-(0.5d0*prob_lo(1)+0.5d0*prob_hi(1)))**2 &
                            + (y-(0.6d0*prob_lo(2)+0.4d0*prob_hi(2)))**2 &
                            + (z-(0.5d0*prob_lo(3)+0.5d0*prob_hi(3)))**2)
                r2 = sqrt (   (x-(0.5d0*prob_lo(1)+0.5d0*prob_hi(1)))**2 &
                            + (y-(0.4d0*prob_lo(2)+0.6d0*prob_hi(2)))**2 &
                            + (z-(0.5d0*prob_lo(3)+0.5d0*prob_hi(3)))**2)
                
                ! set c using Gaussian bump
                c(i,j,k,1) = c_init(1,1)*exp(-75.d0*r1**2)
                c(i,j,k,2) = c_init(2,2)*exp(-75.d0*r2**2)
                
             enddo
          enddo
       enddo

    case (7)

       !=============================================================
       ! 1D sin wave of c1
       !=============================================================

       u = 0.d0
       v = 0.d0

       ! note: c(:,:,:,3) will be computed below to enforce sum(c)=1
       c(:,:,:,2) = 0.2d0

       do k=lo(3),hi(3)
          z = prob_lo(3) + dx(3) * (dble(k)+0.5d0)
          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2) * (dble(j)+0.5d0)
             do i=lo(1),hi(1)
                x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)

                c(i,j,k,1) = 0.2d0 + 0.01d0*sin(2.d0*M_PI*x/L(1))
                
             enddo
          enddo
       enddo

    case (11)

       !=============================================================
       ! Discontinuous square in the central 25% of domain
       ! c=c_init(1,:) inside; c=c_init(2,:) outside
       !=============================================================

       u = 0.d0
       v = 0.d0
       w = 0.d0

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)

                ! initialize c to a square region
                if (i .ge. n_cells(1)/4 .and. i .le. 3*n_cells(1)/4-1 .and. &
                    j .ge. n_cells(2)/4 .and. j .le. 3*n_cells(2)/4-1 .and. &
                    k .ge. n_cells(3)/4 .and. k .le. 3*n_cells(3)/4-1) then
                   c(i,j,k,1:nspecies) = c_init(1,1:nspecies)
                else
                   c(i,j,k,1:nspecies) = c_init(2,1:nspecies)
                end if

             enddo
          enddo
       enddo

    case (12)

       !=============================================================
       ! Gaussian bubble centered in domain
       ! c=c_init(1,:) inside; c=c_init(2,:) outside
       ! lo- and hi-y walls move with prescribed velocity,
       ! see inhomogeneous_bc_val.f90
       !=============================================================

       u = 0.d0
       v = 0.d0
       w = 0.d0

       do k=lo(3),hi(3)
          z = prob_lo(3) + dx(3) * (dble(k)+0.5d0) - 0.5d0*(prob_lo(3)+prob_hi(3))
          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2) * (dble(j)+0.5d0) - 0.5d0*(prob_lo(2)+prob_hi(2))
             do i=lo(1),hi(1)
                x = prob_lo(1) + dx(1) * (dble(i)+0.5d0) - 0.5d0*(prob_lo(1)+prob_hi(1))

                r = sqrt (x**2 + y**2 + z**2)

                ! set c using Gaussian bump
                c(i,j,k,1:nspecies-1) = c_init(1,1:nspecies-1)*exp(-75.d0*r**2)

             enddo
          enddo
       enddo

    case (13)

       !=============================================================
       ! stratified multispecies due to barodiffusion
       ! approximate analytical steady solution
       ! assumes the final species is the dominant component
       !=============================================================

       u = 0.d0
       v = 0.d0
       w = 0.d0

       do n=1,nspecies-1
          m_e = (rhobar(nspecies)/rhobar(n) - 1.d0)*molmass(n)
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                y = prob_lo(2) + dx(2)*(dble(j)+0.5d0)
                do i=lo(1),hi(1)
                   c(i,j,k,n) = c_init(1,n)*exp(-m_e*grav(2)*y/(k_B*T_init(1)))
                enddo
             enddo
          enddo
       enddo

    case (14)

       !=============================================================
       ! stratified multispecies due to thermodiffusion
       ! approximate analytical steady solution
       ! assumes nspecies=3
       ! assumes the final species is the dominant component
       !=============================================================

       if (nspecies .ne. 3) then
          call bl_error("prob_type=14 requires nspecies=3")
       end if

       gradToverT = (T_init(2)-T_init(1))/(T_init(1)*(prob_hi(2)-prob_lo(2)))

       u = 0.d0
       v = 0.d0
       w = 0.d0

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2)*(dble(j)+0.5d0)
             do i=lo(1),hi(1)
             
                c(i,j,k,1) = c_init(1,1)*exp((Dtherm(1)-Dtherm(3))*gradToverT*y/Dbar(2))
                c(i,j,k,2) = c_init(1,2)*exp((Dtherm(2)-Dtherm(3))*gradToverT*y/Dbar(3))

             enddo
          enddo
       enddo

    case default

       call bl_error("Desired prob_type not supported in 3D")

    end select

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ! set final c_i such that sum(c_i) = 1 to within roundoff
             sum = 0.d0
             do n=1,nspecies-1
                sum = sum + c(i,j,k,n)
             end do
             c(i,j,k,nspecies) = 1.d0 - sum

             ! calculate rho_total from eos
             if (algorithm_type .eq. 6) then
                rho_total = rho0
             else
                sum = 0.d0
                do n=1,nspecies
                   ! sum represents rhoinv
                   sum = sum + c(i,j,k,n)/rhobar(n)
                end do
                rho_total = 1.d0/sum
             end if

             ! add mass fluctuations
             if (abs(initial_variance_mass) .gt. 0.d0) then
                call add_mass_fluctuations(c(i,j,k,1:nspecies),dx,abs(initial_variance_mass),rho_total)

                if (algorithm_type /= 6) then ! Recompute rho from EOS with fluctuations added
                   sum = 0.d0
                   do n=1,nspecies
                      ! sum represents rhoinv
                      sum = sum + c(i,j,k,n)/rhobar(n)
                   end do
                   rho_total = 1.d0/sum
                end if
             end if

             ! calculate rho_i
             rho(i,j,k,1:nspecies) = rho_total*c(i,j,k,1:nspecies)

          end do
       end do
    end do

  end subroutine init_rho_and_umac_3d

  ! Add local natural random fluctuations to the mass fractions assuming an ideal mixture
  ! This routine does not require singling out a solvent but only works for ideal mixtures
  subroutine add_mass_fluctuations(c,dx,variance,rho_tot)
    real(dp_t), intent(inout) :: c(nspecies) ! Mass fractions
    real(dp_t), intent(in)    :: dx(:)
    real(dp_t), intent(in)    :: rho_tot, variance

    ! local variables
    real(dp_t), dimension(nspecies) :: z, dc, c0
    real(dp_t) :: factor

    real(dp_t) :: tmp
    integer    :: i,j

    factor = sqrt(variance/product(dx(1:MAX_SPACEDIM)))

    c0 = c ! Store this  

    ! construct random vector z having nspecies N(0,1) random variables
    if (use_bl_rng) then
       call NormalRNGs(z, nspecies, engine=rng_eng_init_mass%p)
    else
       call NormalRNGs(z, nspecies)
    end if
    
    ! for ideal mixture, dw = 1/sqrt(rho)*(I-w*1^T)*P*sqrt(W)*sqrt(M)*N(0,1)
    ! Here P is a projection onto charge-neutrality, if this is necessary
    ! P=(I-W*M*z*z^T / (z^T*W*M*z))
    dc = factor*sqrt(c0*molmass(1:nspecies)/rho_tot)*z ! Unprojected mass fluctuations

    if(use_charged_fluid .and. electroneutral) then ! Project onto z^T*w=0 first
       factor = sum(molmass(1:nspecies)*c0*(charge_per_mass(1:nspecies)**2)) ! Related to Debye length
       ! Project onto z^T*w=0
       dc = dc - sum(dc*charge_per_mass(1:nspecies))/factor * molmass(1:nspecies)*c0*charge_per_mass(1:nspecies)
       ! Project onto 1^T*w=1
       dc = dc - sum(dc)*c0 ! Make it sum to zero
       c = c + dc ! add the fluctuations -- this can produce negative values
       ! We don't try to make them positive since charge neutrality has to be obeyed strictly
    else
       dc = dc - sum(dc)*c0 ! Make it sum to zero
       c = c + dc ! add the fluctuations        
       ! replace negative values by zero and normalize
       c = max(c,0.d0)
       c = c/sum(c)  
    end if

  end subroutine add_mass_fluctuations

  subroutine init_Temp(Temp,dx,time,the_bc_level)

    type(multifab) , intent(inout) :: Temp(:)            
    real(kind=dp_t), intent(in   ) :: dx(:,:)           
    real(kind=dp_t), intent(in   ) :: time 
    type(bc_level) , intent(in   ) :: the_bc_level(:)
 
    ! local variables
    integer                        :: lo(Temp(1)%dim), hi(Temp(1)%dim)
    integer                        :: dm, ng, i, n, nlevs
    real(kind=dp_t), pointer       :: dp1(:,:,:,:)  ! pointer for Temp 

    type(bl_prof_timer), save :: bpt

    call build(bpt,"init_Temp")

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
       ! FIXME: it is more accurate to write loops over the ghost cells below to fill in the functional values
       !        than average the interior+ghost -> ghost to store the Dirichlet value
       call multifab_coefbc(Temp(n),1,1,the_bc_level(n))

    end do

    call destroy(bpt)

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
    select case (temp_type)

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

       Temp = T_init(1)

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
    select case (temp_type)
    
    case(0) 
    !================================================================================
    ! Thermodynamic equilibrium
    !================================================================================
       Temp = T_init(1)
    
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

       Temp = T_init(1)

   end select
   
  end subroutine init_Temp_3d

end module init_lowmach_module
