module init_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use convert_stag_module
  use bc_module
  use bl_constants_module
  use probin_binarylm_module, only: rhobar, diff_coef, &
                                    smoothing_width, c_init, material_properties, &
                                    u_init
  use probin_common_module , only: prob_lo, prob_hi, prob_type, visc_type, visc_coef, &
                                   diff_type, n_cells

  implicit none

  private

  public :: init, compute_eta, compute_chi, compute_kappa

 ! prob_type controls the initial condition and selects specific problem to solve
 ! Positive means transport coefficients depend as rational functions on concentration
 ! e.g., eta=visc_coeff*(a+b*conc)/(1+c*conc)
 ! Negative means transport coefficients depend as polynomial functions on concentration
 ! e.g., eta=visc_coeff*(1+a*conc+b*conc^2+c*conc^3)
 ! KEY:
 ! 0 - linear gradient in c; c_init(1) at lo-y wall, c_init(2) at hi-y wall
 ! 1 = spherical bubble with c_init(1) in the interior, c_init(2) on the exterior
 ! 2 = bilayer interface (stripe)
   ! the lower third and upper third of the domain (in y) has c_init(1)
   ! the middle third of the domain has c_init(2)
 ! 3 = one fluid on top of another
    ! c_init(1), u_init(1) in lower half of domain (in y)
    ! c_init(2), u_init(2) in upper half
 ! 4 = Bell, Colella, Glaz 1989 jet in a doubly period geometry
 ! 5 = Kelvin-Helmholtz
 ! 6 = Exact solution for constant coefficient (Taylor vortex travelling wave)
 ! 7 = Steady state (diff_flux=0 i.e. rho(c(y))*chi(c(y))*dc/dy=const) -- NOT IMPLEMENTED, for now the same as prob_type=0
       ! See prob_type=13 in lowMach_explicit/exact_solutions.f90 for partial implementation

contains

  ! Important note: For periodic boundaries, the init routines should fill out
  ! *both* sides of the domain with values, even though this is duplicate information
  ! We ensure the two sides are bitwise identical, but low and/or high sides may win

  subroutine init(m,s,p,dx,mla,time)

    type(multifab) , intent(inout) :: m(:,:),s(:),p(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(ml_layout), intent(in   ) :: mla
    real(kind=dp_t), intent(in   ) :: time

    real(kind=dp_t), pointer :: mxp(:,:,:,:), myp(:,:,:,:), mzp(:,:,:,:)
    real(kind=dp_t), pointer :: sop(:,:,:,:), pp(:,:,:,:)
    integer :: lo(mla%dim),hi(mla%dim)
    integer :: nlevs,n,i,ng_m,ng_s,ng_p,dm

    nlevs = mla%nlevel
    dm = mla%dim
    
    ng_m = m(1,1)%ng
    ng_s = s(1)%ng
    ng_p = p(1)%ng

    do n=1,nlevs
       do i = 1, nfabs(s(n))
          mxp => dataptr(m(n,1),i)
          myp => dataptr(m(n,2),i)
          sop => dataptr(s(n),i)
          pp  => dataptr(p(n),i)
          lo =  lwb(get_box(s(n),i))
          hi =  upb(get_box(s(n),i))
          select case (dm)
          case (2)
             call init_2d(mxp(:,:,1,1), myp(:,:,1,1), &
                          sop(:,:,1,:), pp(:,:,1,1), &
                          lo, hi, ng_m, ng_s, ng_p, dx(n,:), time)
          case (3)
             mzp => dataptr(m(n,3),i)
             call init_3d(mxp(:,:,:,1), myp(:,:,:,1), mzp(:,:,:,1), &
                          sop(:,:,:,:), pp(:,:,:,1), &
                          lo, hi, ng_m, ng_s, ng_p, dx(n,:), time)
          end select
       end do

       ! For periodic boundaries, ensure the low and high side are consistent:
       ! Note: multifab_internal_sync compares the box number of the two boxes 
       ! with overlapping values and the data on the box with lower number wins. 
       do i=1,dm
          call multifab_internal_sync(m(n,i))
       end do

    enddo

  end subroutine init

  subroutine init_2d(mx,my,s,p,lo,hi,ng_m,ng_s,ng_p,dx,time)

    integer        , intent(in   ) :: lo(:), hi(:), ng_m, ng_s, ng_p
    real(kind=dp_t), intent(inout) ::  mx(lo(1)-ng_m:,lo(2)-ng_m:)
    real(kind=dp_t), intent(inout) ::  my(lo(1)-ng_m:,lo(2)-ng_m:)
    real(kind=dp_t), intent(inout) ::   s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(inout) ::   p(lo(1)-ng_p:,lo(2)-ng_p:)
    real(kind=dp_t), intent(in   ) :: dx(:),time

    ! local
    integer :: i,j,mid
    real(kind=dp_t) :: x,y,y1,y2,r,dy,c_loc,u_loc
    real(kind=dp_t) :: one_third_domain1,one_third_domain2
    real(kind=dp_t) :: cosxt,cosyt,freq,pfac,pfreq
    real(kind=dp_t) :: sinxt,sinyt,ucst,ufac,vcst,xm,xp,ym,yp,rand

    select case (abs(prob_type))
    case (0,7)

       ! linear gradient in c; c_init(1) at lo-y wall, c_init(2) at hi-y wall

       mx = 0.d0
       my = 0.d0

       p = 0.d0

       do j=lo(2),hi(2)

          ! compute distance from bottom of domain
          dy = dx(2)*(j+0.5d0)

          ! linear gradient in c
          c_loc = c_init(1) + (c_init(2)-c_init(1))*dy/(prob_hi(2)-prob_lo(2))

          ! compute rho with the eos
          s(lo(1):hi(1),j,1) = 1.0d0/(c_loc/rhobar(1)+(1.0d0-c_loc)/rhobar(2))

          ! compute rho*c
          s(lo(1):hi(1),j,2) = s(lo(1):hi(1),j,1)*c_loc

       end do

    case (1)

       ! spherical bubble with c_init(1) in the interior, c_init(2) on the exterior
       ! centered in domain where smoothing_width has units of GRID CELLS

       mx = 0.d0
       my = 0.d0

       p = 0.d0

       do j=lo(2),hi(2)
          y = prob_lo(2) + dx(2) * (dble(j)+0.5d0) - 0.5d0*(prob_lo(2)+prob_hi(2))
          do i=lo(1),hi(1)
             x = prob_lo(1) + dx(1) * (dble(i)+0.5d0) - 0.5d0*(prob_lo(1)+prob_hi(1))

             r = sqrt (x**2 + y**2)

             ! set c using tanh smoothing
             s(i,j,2) = c_init(1) + 0.5d0*(c_init(2)-c_init(1))* &
                  (1.d0 + tanh((r-2.5d0*smoothing_width*dx(1))/(smoothing_width*dx(1))))

             ! compute rho using eos
             s(i,j,1) = 1.0d0/(s(i,j,2)/rhobar(1)+(1.0d0-s(i,j,2))/rhobar(2))

             ! compute rho*c
             s(i,j,2) = s(i,j,1)*s(i,j,2)
          enddo
       enddo

    case (2)

       ! bilayer interface (stripe)
       ! the lower third and upper third of the domain (in y) has c_init(1)
       ! the middle third of the domain has c_init(2)
       ! smoothing_width has units of GRID CELLS
       ! if smoothing width is 0, use finite-volume averaging of sharp interface

       mx = 0.d0
       my = 0.d0

       p = 0.d0

       one_third_domain1=2.0d0/3.0d0*prob_lo(2)+1.0d0/3.0d0*prob_hi(2)
       one_third_domain2=1.0d0/3.0d0*prob_lo(2)+2.0d0/3.0d0*prob_hi(2)

       do j=lo(2),hi(2)
          y1 =(prob_lo(2) + dx(2)*(dble(j)+0.5d0) - one_third_domain1)
          y2 =(prob_lo(2) + dx(2)*(dble(j)+0.5d0) - one_third_domain2)
        
          ! tanh smoothing
          if(abs(smoothing_width)>epsilon(1.0d0)) then
             c_loc = c_init(1)+ 0.5d0*(c_init(2)-c_init(1))*&
                  (tanh(y1/(smoothing_width*dx(2))) - tanh(y2/(smoothing_width*dx(2))))
             s(lo(1):hi(1),j,1) = 1.0d0/(c_loc/rhobar(1)+(1.0d0-c_loc)/rhobar(2))
             s(lo(1):hi(1),j,2) = s(lo(1):hi(1),j,1)*c_loc
          else
             ! Try to initialize exactly as we do in the HDMD simulations,
             ! with finite-volume averaging of sharp interface
             if((y1<-0.5d0*dx(2)).or.(y2>0.5d0*dx(2))) then
                s(lo(1):hi(1),j,2) = 0
                s(lo(1):hi(1),j,1) = rhobar(2)
             else if((y1>0.5d0*dx(2)).and.(y2<-0.5d0*dx(2))) then
                s(lo(1):hi(1),j,2) = rhobar(1)
                s(lo(1):hi(1),j,1) = rhobar(1)
             else if(y1 <= 0.5d0*dx(2)) then
                s(lo(1):hi(1),j,2) = (max(0.0d0,min(0.5d0+y1/dx(2),1.0d0)))*rhobar(1)
                s(lo(1):hi(1),j,1) = s(lo(1):hi(1),j,2) &
                     + (1.0d0-max(0.0d0,min(0.5d0+y1/dx(2),1.0d0)))*rhobar(2)
             else 
                s(lo(1):hi(1),j,2) = (1.0d0-max(0.0d0,min(0.5d0+y2/dx(2),1.0d0)))*rhobar(1)
                s(lo(1):hi(1),j,1) = s(lo(1):hi(1),j,2) + (max(0.0d0,min(0.5d0+y2/dx(2),1.0d0)))*rhobar(2)
             end if
          end if
       enddo

    case (3)

       ! one fluid on top of another
       ! c_init(1), u_init(1) in lower half of domain (in y)
       ! c_init(2), u_init(2) in upper half

       my = 0.d0

       p = 0.d0

       ! middle of domain
       y1 = (prob_lo(2)+prob_hi(2)) / 2.d0

       if(abs(smoothing_width)>epsilon(1.d0)) then

          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2)*(dble(j)+0.5d0) - y1

             ! smoothed version
             ! c_init(1) in lower half of domain (in y)
             ! c_init(2) in upper half
             c_loc = c_init(1) + (c_init(2)-c_init(1))*0.5d0*(tanh(y/(smoothing_width*dx(2)))+1.d0)
             s(lo(1):hi(1),j,1) = 1.0d0/(c_loc/rhobar(1)+(1.0d0-c_loc)/rhobar(2))
             s(lo(1):hi(1),j,2) = s(lo(1):hi(1),j,1)*c_loc

          end do

       else

          ! c_init(1) in lower half of domain (in y)
          ! c_init(2) in upper half

          do j=lo(2),hi(2)
             y = prob_lo(2) + (j+0.5d0)*dx(2)

             if (y .lt. y1) then
                c_loc = c_init(1)
             else
                c_loc = c_init(2)
             end if

             s(lo(1):hi(1),j,1) = 1.0d0/(c_loc/rhobar(1)+(1.0d0-c_loc)/rhobar(2))
             s(lo(1):hi(1),j,2) = s(lo(1):hi(1),j,1)*c_loc
          
          end do

       end if

       do j=lo(2),hi(2)
          y = prob_lo(2) + (j+0.5d0)*dx(2)
          
          if (y .lt. y1) then
             u_loc = u_init(1)
          else
             u_loc = u_init(2)
          end if

          mx(lo(1):hi(1)+1,j) = u_loc*s(lo(1),j,1)
          
       end do

    case (4) 

       ! Bell, Colella, Glaz 1989
       ! jet in a doubly period geometry

       ! constant density
       s(:,:,1) = 1.d0/30.d0

       ! tracer
       do j=lo(2),hi(2)
          y = prob_lo(2) + (j+0.5d0)*dx(2)
          if (y .le. 0.5d0) then
             s(:,j,2) = 0.d0
          else
             s(:,j,2) = 1.d0/30.d0
          end if
       end do

       ! x velocity
       do j=lo(2),hi(2)
          y = prob_lo(2) + (j+0.5d0)*dx(2)
          if (y .le. 0.5d0) then
             mx(:,j) = tanh(30.d0*(y-0.25d0)) / 30.d0
          else
             mx(:,j) = tanh(30.d0*(0.75d0-y)) / 30.d0
          end if
       end do

       ! y-velocity
       do i=lo(1),hi(1)
          x = prob_lo(1) + (i+0.5d0)*dx(1)
          my(i,:) = (0.05d0/30.d0) * sin(2.d0*M_PI*x)
       end do

    case (5)

       ! Kelvin-Helmholtz

       ! density
       do j=lo(2),hi(2)
          y = prob_lo(2) + (j+0.5d0)*dx(2)
          if (y .le. 0.5d0*(prob_hi(2)+prob_lo(2))) then
             s(:,j,1) = 10.d0
          else
             s(:,j,1) = 1.d0
          end if
       end do

       ! density perturbation
       mid = n_cells(2)/2


       if (lo(2) .le. mid .and. hi(2) .ge. mid) then
          do i=lo(1),hi(1)
             call random_number(rand)
             s(i,mid,1) = s(i,mid,1)*(1.d0+0.01d0*rand)
          end do
!          s(0,mid,1) = 1.01d0*s(0,mid,1)
       end if

       ! tracer
       do j=lo(2),hi(2)
          y = prob_lo(2) + (j+0.5d0)*dx(2)
          if (y .le. 0.5d0*(prob_hi(2)+prob_lo(2))) then
             s(:,j,2) = 0.d0
          else
             s(:,j,2) = 1.d0
          end if
       end do

       ! x-momentum
       do j=lo(2),hi(2)
          y = prob_lo(2) + (j+0.5d0)*dx(2)
          if (y .le. 0.5d0*(prob_hi(2)+prob_lo(2))) then
             mx(:,j) = 0.d0
          else
             mx(:,j) = 1.d0
          end if
       end do

       ! y-momentum
       my = 0.d0

    case (6)

       ! traveling wave with exact solution

       s(:,:,1) = 1.d0

       freq  = 2.d0*M_PI
       pfreq = 4.d0*M_PI

       ucst = 0.75d0
       vcst = 0.75d0
       ufac = dexp(-2.0d0*freq*freq*time*visc_coef)
       pfac = dexp(-4.0d0*freq*freq*time*visc_coef)/64.0d0

       do j=lo(2),hi(2)
          y = dx(2) * (dble(j) + half) - vcst*time
          yp = y+0.5d0*dx(2)
          ym = y-0.5d0*dx(2)
          do i=lo(1),hi(1)
             x = dx(1) * (dble(i) + half) - ucst*time

             xp = x+0.5d0*dx(1)
             xm = x-0.5d0*dx(1)

             cosxt =  (sin(freq*xp)-sin(freq*xm))/(freq*dx(1))
             sinxt = -(cos(freq*xp)-cos(freq*xm))/(freq*dx(1))
             cosyt =  (sin(freq*yp)-sin(freq*ym))/(freq*dx(2))
             sinyt = -(cos(freq*yp)-cos(freq*ym))/(freq*dx(2))

             s(i,j,2) = ucst + ufac * 0.25d0*cosxt*sinyt
          enddo
       enddo

       do j=lo(2),hi(2)
          y = dx(2) * (dble(j) + half) - vcst*time
          yp = y+0.5d0*dx(2)
          ym = y-0.5d0*dx(2)
          do i=lo(1),hi(1)+1
             x = dx(1) * (dble(i) ) - ucst*time

             xp = x+0.5d0*dx(1)
             xm = x-0.5d0*dx(1)

             cosxt =  (sin(freq*xp)-sin(freq*xm))/(freq*dx(1))
             sinxt = -(cos(freq*xp)-cos(freq*xm))/(freq*dx(1))
             cosyt =  (sin(freq*yp)-sin(freq*ym))/(freq*dx(2))
             sinyt = -(cos(freq*yp)-cos(freq*ym))/(freq*dx(2))

             mx(i,j) = ucst + ufac * 0.25d0*cosxt*sinyt
          enddo
       enddo

       do j=lo(2),hi(2)+1
          y = dx(2) * (dble(j)) - vcst*time
          yp = y+0.5d0*dx(2)
          ym = y-0.5d0*dx(2)
          do i=lo(1),hi(1)
             x = dx(1) * (dble(i) + half) - ucst*time

             xp = x+0.5d0*dx(1)
             xm = x-0.5d0*dx(1)

             cosxt =  (sin(freq*xp)-sin(freq*xm))/(freq*dx(1))
             sinxt = -(cos(freq*xp)-cos(freq*xm))/(freq*dx(1))
             cosyt =  (sin(freq*yp)-sin(freq*ym))/(freq*dx(2))
             sinyt = -(cos(freq*yp)-cos(freq*ym))/(freq*dx(2))

             my(i,j) = vcst - ufac * 0.25d0*sinxt*cosyt
          enddo
       enddo

    case default

       call bl_error("init_2d: invalid prob_type")

    end select

  end subroutine init_2d

  subroutine init_3d(mx,my,mz,s,p,lo,hi,ng_m,ng_s,ng_p,dx,time)

    integer        , intent(in   ) :: lo(:), hi(:), ng_m, ng_s, ng_p
    real(kind=dp_t), intent(inout) ::  mx(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real(kind=dp_t), intent(inout) ::  my(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real(kind=dp_t), intent(inout) ::  mz(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real(kind=dp_t), intent(inout) ::   s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(inout) ::   p(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real(kind=dp_t), intent(in   ) :: dx(:),time

    ! local
    integer :: i,j,k
    real(kind=dp_t) :: x,y,y1,y2,z,r,dy,c_loc
    real(kind=dp_t) :: one_third_domain1,one_third_domain2

    select case (abs(prob_type))
    case (0,7)

       ! linear gradient in c; c_init(1) at lo-y wall, c_init(2) at hi-y wall

       mx = 0.d0
       my = 0.d0
       mz = 0.d0

       p = 0.d0

       do j=lo(2),hi(2)

          ! compute distance from bottom of domain
          dy = dx(2)*(j+0.5d0)

          ! linear gradient in c
          c_loc = c_init(1) + (c_init(2)-c_init(1))*dy/(prob_hi(2)-prob_lo(2))

          ! compute rho with the eos
          s(lo(1):hi(1),j,lo(3):hi(3),1) = 1.0d0/(c_loc/rhobar(1)+(1.0d0-c_loc)/rhobar(2))

          ! compute rho*c
          s(lo(1):hi(1),j,lo(3):hi(3),2) = s(lo(1):hi(1),j,lo(3):hi(3),1)*c_loc

       end do

    case (1)

       ! spherical bubble with c_init(1) in the interior, c_init(2) on the exterior
       ! centered in domain where smoothing_width has units of GRID CELLS

       mx = 0.d0
       my = 0.d0
       mz = 0.d0

       p = 0.d0

       do k=lo(3),hi(3)
          z = prob_lo(3) + dx(3) * (dble(k)+0.5d0) - 0.5d0*(prob_lo(3)+prob_hi(3))
          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2) * (dble(j)+0.5d0) - 0.5d0*(prob_lo(2)+prob_hi(2))
             do i=lo(1),hi(1)
                x = prob_lo(1) + dx(1) * (dble(i)+0.5d0) - 0.5d0*(prob_lo(1)+prob_hi(1))
                
                r = sqrt (x**2 + y**2 + z**2)

                ! set c using tanh smoothing
                s(i,j,k,2) = c_init(1) + 0.5d0*(c_init(2)-c_init(1))* &
                     (1.d0 + tanh((r-2.0d0*smoothing_width*dx(1))/(smoothing_width*dx(1))))
                
                ! compute rho using eos
                s(i,j,k,1) = 1.0d0/(s(i,j,k,2)/rhobar(1)+(1.0d0-s(i,j,k,2))/rhobar(2))

                ! compute rho*c
                s(i,j,k,2) = s(i,j,k,1)*s(i,j,k,2)
             enddo
          enddo
       end do

    case (2)

       ! bilayer interface (stripe)
       ! the lower third and upper third of the domain (in y) has c_init(1)
       ! the middle third of the domain has c_init(2)
       ! smoothing_width has units of GRID CELLS
       ! if smoothing width is 0, use finite-volume averaging of sharp interface

       mx = 0.d0
       my = 0.d0
       mz = 0.d0

       p = 0.d0

       one_third_domain1=2.0d0/3.0d0*prob_lo(2)+1.0d0/3.0d0*prob_hi(2)
       one_third_domain2=1.0d0/3.0d0*prob_lo(2)+2.0d0/3.0d0*prob_hi(2)

       do j=lo(2),hi(2)
          y1 =(prob_lo(2) + dx(2)*(dble(j)+0.5d0) - one_third_domain1)
          y2 =(prob_lo(2) + dx(2)*(dble(j)+0.5d0) - one_third_domain2)
        
          ! tanh smoothing
          if(abs(smoothing_width)>epsilon(1.0d0)) then
             c_loc = c_init(1)+ 0.5d0*(c_init(2)-c_init(1))*&
                  (tanh(y1/(smoothing_width*dx(2))) - tanh(y2/(smoothing_width*dx(2))))
             s(lo(1):hi(1),j,lo(3):hi(3),1) = 1.0d0/(c_loc/rhobar(1)+(1.0d0-c_loc)/rhobar(2))
             s(lo(1):hi(1),j,lo(3):hi(3),2) = s(lo(1):hi(1),j,lo(3):hi(3),1)*c_loc
          else
             ! Try to initialize exactly as we do in the HDMD simulations,
             ! with finite-volume averaging of sharp interface
             if((y1<-0.5d0*dx(2)).or.(y2>0.5d0*dx(2))) then
                s(lo(1):hi(1),j,lo(3):hi(3),2) = 0
                s(lo(1):hi(1),j,lo(3):hi(3),1) = rhobar(2)
             else if((y1>0.5d0*dx(2)).and.(y2<-0.5d0*dx(2))) then
                s(lo(1):hi(1),j,lo(3):hi(3),2) = rhobar(1)
                s(lo(1):hi(1),j,lo(3):hi(3),1) = rhobar(1)
             else if(y1 <= 0.5d0*dx(2)) then
                s(lo(1):hi(1),j,lo(3):hi(3),2) = (max(0.0d0,min(0.5d0+y1/dx(2),1.0d0)))*rhobar(1)
                s(lo(1):hi(1),j,lo(3):hi(3),1) = s(lo(1):hi(1),j,lo(3):hi(3),2) &
                     + (1.0d0-max(0.0d0,min(0.5d0+y1/dx(2),1.0d0)))*rhobar(2)
             else 
                s(lo(1):hi(1),j,lo(3):hi(3),2) = (1.0d0-max(0.0d0,min(0.5d0+y2/dx(2),1.0d0)))*rhobar(1)
                s(lo(1):hi(1),j,lo(3):hi(3),1) = s(lo(1):hi(1),j,lo(3):hi(3),2) &
                     + (max(0.0d0,min(0.5d0+y2/dx(2),1.0d0)))*rhobar(2)
             end if
          end if
       enddo

    case (3)

       ! one fluid on top of another
       ! c_init(1) in lower half of domain (in y)
       ! c_init(2) in upper half

       mx = 0.d0
       my = 0.d0
       mz = 0.d0

       p = 0.d0

       ! middle of domain
       y1 = (prob_lo(2)+prob_hi(2)) / 2.d0

       do j=lo(2),hi(2)
          y = prob_lo(2) + (j+0.5d0)*dx(2)

          if (y .lt. y1) then
             s(lo(1):hi(1),j,lo(3):hi(3),2) = c_init(1)
             s(lo(1):hi(1),j,lo(3):hi(3),1) = 1.0d0/(c_init(1)/rhobar(1)+(1.0d0-c_init(1))/rhobar(2))
          else
             s(lo(1):hi(1),j,lo(3):hi(3),2) = c_init(2)
             s(lo(1):hi(1),j,lo(3):hi(3),1) = 1.0d0/(c_init(2)/rhobar(1)+(1.0d0-c_init(2))/rhobar(2))   
          end if
          s(lo(1):hi(1),j,lo(3):hi(3),2) = s(lo(1):hi(1),j,lo(3):hi(3),1)*s(lo(1):hi(1),j,lo(3):hi(3),2)
          
       end do



    case default

       call bl_error("init_3d: invalid prob_type")
          
    end select

  end subroutine init_3d

!==============================================================
! Transport coefficients concentration dependence
!==============================================================

  ! Here we assume it to be of a simple rational or polynomial form depending on sign of prob_type
  subroutine compute_coeff_local(indx,conc,coeff)
    integer, intent(in) :: indx
    real(kind=dp_t), intent(in) ::  conc
    real(kind=dp_t), intent(out) :: coeff
     
    if(prob_type>=0) then ! Rational dependence
      coeff = (material_properties(1,indx) + material_properties(2,indx)*conc) / &
              (1.d0                        + material_properties(3,indx)*conc)
    else ! Polynomial dependence
      coeff = (1.0d0 + material_properties(1,indx)*conc    + &
                       material_properties(2,indx)*conc**2 + &
                       material_properties(3,indx)*conc**3)
    end if                                   
  
  end subroutine
  
  subroutine compute_chi(mla,chi,chi_fc,prim,dx,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: chi(:)
    type(multifab) , intent(inout) :: chi_fc(:,:)
    type(multifab) , intent(in   ) :: prim(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    integer :: nlevs,dm,i,n,ng_c,ng_p
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: cp(:,:,:,:)
    real(kind=dp_t), pointer :: pp(:,:,:,:)

    ng_c = chi(1)%ng
    ng_p = prim(1)%ng

    nlevs = mla%nlevel
    dm = mla%dim

    if(ng_p<ng_c) then
       call bl_error('ng_p must be no less than ng_c')
    end if

    do n=1,nlevs
       do i=1,nfabs(chi(n))
          cp => dataptr(chi(n), i)
          pp => dataptr(prim(n), i)
          lo = lwb(get_box(chi(n), i))
          hi = upb(get_box(chi(n), i))
          select case (dm)
          case (2)
             call compute_chi_2d(cp(:,:,1,1),ng_c,pp(:,:,1,:),ng_p,lo,hi,dx(n,:))
          case (3)
             call compute_chi_3d(cp(:,:,:,1),ng_c,pp(:,:,:,:),ng_p,lo,hi,dx(n,:))
          end select
       end do
    end do

    call average_cc_to_face(nlevs,chi,chi_fc,1,tran_bc_comp,1,the_bc_level)

  end subroutine compute_chi
  
  subroutine compute_chi_2d(chi,ng_c,prim,ng_p,lo,hi,dx)

    ! compute chi in valid AND ghost regions
    ! the ghost cells in prim have already been filled properly

    integer        , intent(in   ) :: lo(:), hi(:), ng_c, ng_p
    real(kind=dp_t), intent(inout) ::  chi(lo(1)-ng_c:,lo(2)-ng_c:)
    real(kind=dp_t), intent(in   ) :: prim(lo(1)-ng_p:,lo(2)-ng_p:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! local
    real(kind=dp_t) :: conc
    integer :: i,j

    if (diff_type > 0) then

       chi = diff_coef

    else

       do j=lo(2)-ng_c,hi(2)+ng_c
       do i=lo(1)-ng_c,hi(1)+ng_c
          conc = max(min(prim(i,j,2), 1.d0), 0.d0)
          call compute_coeff_local(1, conc, chi(i,j))
          chi(i,j) = diff_coef*chi(i,j)
          !write(101,*) j, chi(i,j)
       end do
       end do

    end if

  end subroutine compute_chi_2d

  subroutine compute_chi_3d(chi,ng_c,prim,ng_p,lo,hi,dx)

    ! compute chi in valid AND ghost regions
    ! the ghost cells in prim have already been filled properly

    integer        , intent(in   ) :: lo(:), hi(:), ng_c, ng_p
    real(kind=dp_t), intent(inout) ::  chi(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:)
    real(kind=dp_t), intent(in   ) :: prim(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! local
    real(kind=dp_t) :: conc
    integer :: i,j,k

    if (diff_type > 0) then

       chi = diff_coef

    else

       do k=lo(3)-ng_c,hi(3)+ng_c
       do j=lo(2)-ng_c,hi(2)+ng_c
       do i=lo(1)-ng_c,hi(1)+ng_c
          conc = max(min(prim(i,j,k,2), 1.d0), 0.d0)
          call compute_coeff_local(1, conc, chi(i,j,k))
          chi(i,j,k) = diff_coef*chi(i,j,k)
       end do
       end do
       end do

    end if

  end subroutine compute_chi_3d

  subroutine compute_eta(mla,eta,eta_ed,prim,dx,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: eta(:)
    type(multifab) , intent(inout) :: eta_ed(:,:) ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(in   ) :: prim(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    integer :: nlevs,dm,i,n,ng_e,ng_p
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: cp(:,:,:,:)
    real(kind=dp_t), pointer :: pp(:,:,:,:)

    ng_e = eta(1)%ng
    ng_p = prim(1)%ng

    nlevs = mla%nlevel
    dm = mla%dim

    if(ng_p<ng_e) then
       call bl_error('ng_p must be no less than ng_e')
    end if

    do n=1,nlevs
       do i=1,nfabs(eta(n))
          cp => dataptr(eta(n), i)
          pp => dataptr(prim(n), i)
          lo = lwb(get_box(eta(n), i))
          hi = upb(get_box(eta(n), i))
          select case (dm)
          case (2)
             call compute_eta_2d(cp(:,:,1,1),ng_e,pp(:,:,1,:),ng_p,lo,hi,dx(n,:))
          case (3)
             call compute_eta_3d(cp(:,:,:,1),ng_e,pp(:,:,:,:),ng_p,lo,hi,dx(n,:))
          end select
       end do
    end do

    if (dm .eq. 2) then
       call average_cc_to_node(nlevs,eta,eta_ed(:,1),1,tran_bc_comp,1,the_bc_level)
    else if (dm .eq. 3) then
       call average_cc_to_edge(nlevs,eta,eta_ed,1,tran_bc_comp,1,the_bc_level)
    end if

  end subroutine compute_eta

  subroutine compute_eta_2d(eta,ng_e,prim,ng_p,lo,hi,dx)

    ! compute eta in valid AND ghost regions
    ! the ghost cells in prim have already been filled properly

    integer        , intent(in   ) :: lo(:), hi(:), ng_e, ng_p
    real(kind=dp_t), intent(inout) ::  eta(lo(1)-ng_e:,lo(2)-ng_e:)
    real(kind=dp_t), intent(in   ) :: prim(lo(1)-ng_p:,lo(2)-ng_p:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! local
    real(kind=dp_t) :: conc
    integer :: i,j

    if (visc_type > 0) then

       eta = visc_coef

    else

       do j=lo(2)-ng_e,hi(2)+ng_e
       do i=lo(1)-ng_e,hi(1)+ng_e
          conc = max(min(prim(i,j,2), 1.d0), 0.d0)
          call compute_coeff_local(2, conc, eta(i,j))
          eta(i,j) = visc_coef*eta(i,j)
          !write(102,*) j, eta(i,j)
       end do
       end do

    end if

  end subroutine compute_eta_2d

  subroutine compute_eta_3d(eta,ng_e,prim,ng_p,lo,hi,dx)

    integer        , intent(in   ) :: lo(:), hi(:), ng_e, ng_p
    real(kind=dp_t), intent(inout) ::  eta(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(in   ) :: prim(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! local
    real(kind=dp_t) :: conc
    integer :: i,j,k

    if (visc_type > 0) then

       eta = visc_coef

    else

       do k=lo(3)-ng_e,hi(3)+ng_e
       do j=lo(2)-ng_e,hi(2)+ng_e
       do i=lo(1)-ng_e,hi(1)+ng_e
          conc = max(min(prim(i,j,k,2), 1.d0), 0.d0)
          call compute_coeff_local(2, conc, eta(i,j,k))
          eta(i,j,k) = visc_coef*eta(i,j,k)
       end do
       end do
       end do

    end if

  end subroutine compute_eta_3d

  subroutine compute_kappa(mla,kappa,prim,dx)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: kappa(:)
    type(multifab) , intent(in   ) :: prim(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)

    integer :: nlevs,dm,i,n,ng_k,ng_p
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: cp(:,:,:,:)
    real(kind=dp_t), pointer :: pp(:,:,:,:)

    ng_k = kappa(1)%ng
    ng_p = prim(1)%ng

    nlevs = mla%nlevel
    dm = mla%dim

    if(ng_p<ng_k) then
       call bl_error('ng_p must be no less than ng_k')
    end if

    do n=1,nlevs
       do i=1,nfabs(kappa(n))
          cp => dataptr(kappa(n), i)
          pp => dataptr(prim(n), i)
          lo = lwb(get_box(kappa(n), i))
          hi = upb(get_box(kappa(n), i))
          select case (dm)
          case (2)
             call compute_kappa_2d(cp(:,:,1,1),ng_k,pp(:,:,1,:),ng_p,lo,hi,dx(n,:))
          case (3)
             call compute_kappa_3d(cp(:,:,:,1),ng_k,pp(:,:,:,:),ng_p,lo,hi,dx(n,:))
          end select
       end do
    end do

  end subroutine compute_kappa

  subroutine compute_kappa_2d(kappa,ng_k,prim,ng_p,lo,hi,dx)

    ! compute kappa in valid AND ghost regions
    ! the ghost cells in prim have already been filled properly

    integer        , intent(in   ) :: lo(:), hi(:), ng_k, ng_p
    real(kind=dp_t), intent(inout) :: kappa(lo(1)-ng_k:,lo(2)-ng_k:)
    real(kind=dp_t), intent(in   ) ::  prim(lo(1)-ng_p:,lo(2)-ng_p:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    if(abs(visc_type)>=3) then
       call bl_error("Bulk viscosity not supported propertly yet")
    else
       kappa=0.0d0
    end if

  end subroutine compute_kappa_2d

  subroutine compute_kappa_3d(kappa,ng_k,prim,ng_p,lo,hi,dx)

    integer        , intent(in   ) :: lo(:), hi(:), ng_k, ng_p
    real(kind=dp_t), intent(inout) :: kappa(lo(1)-ng_k:,lo(2)-ng_k:,lo(3)-ng_k:)
    real(kind=dp_t), intent(in   ) ::  prim(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    if(abs(visc_type)>=3) then
       call bl_error("Bulk viscosity not supported propertly yet")
    else
       kappa=0.0d0
    end if

  end subroutine compute_kappa_3d

end module init_module
