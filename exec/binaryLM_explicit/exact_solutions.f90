module exact_solutions_module

  !   prob_type = 1  !  equilibrium (linear concentration profile)
  !   prob_type = 2  !  traveling wave
  !   prob_type = 3  !  advection diffusion on scalar
  !   prob_type = 4  !  variable density moving vortex
  !   prob_type = 5  !  constant density moving vortex
  !   prob_type = 6  !  mixing with hard walls (half-half)
  !   prob_type = 7  !  mixing for photo-bleaching (FREP) experiments (cyllinder)
  !   prob_type = 8  !  mixing with periodic BCs (stripe third-third-third)
  !   prob_type = 9  !  sawtooth concentration profile
  !   prob_type = 10 !  forced flow with BCs 
  !   prob_type = 11 !  incompressible modes with BCs
  !   prob_type = 12 !  lo-density spherical bubble
  !   prob_type = 13 !  constant c gradient
  
  ! Temporary stuff:
  !   prob_type = 100 ! For internal testing of fluctuation spectrum
  !   prob_type = -21  ! rho*chi=const to get strictly linear concentration profile

  use bl_types
  use bl_constants_module
  use bc_module
  use define_bc_module
  use multifab_module
  use multifab_fill_ghost_module
  use ml_layout_module
  use bl_error_module
  use multifab_physbc_module

  use probin_module, only: prob_type, prob_dir, nscal, ABC_coefs, &
       mode_coefs, diff_coef, visc_coef, rhobar, vel_init, c_init, &
       triangle_coeff, prob_lo, prob_hi, smoothing_width, material_properties, &
       n_cells
       
  implicit none

  private
  public :: exact_2d, exact_3d, exact_s_force_2d, exact_s_force_3d, &
               exact_m_force_2d, exact_m_force_3d, compute_chi, compute_eta, compute_kappa

contains

  ! 
  ! ---------------------------------------------------------------------------
  ! 2D forcing
  ! ---------------------------------------------------------------------------
  ! 

  subroutine exact_2d(mx,my,s,lo,hi,ng_m,ng_s,dx,time)

    integer, intent(in) :: lo(:), hi(:), ng_m, ng_s
    real (kind = dp_t), intent(out) :: mx(lo(1)-ng_m:,lo(2)-ng_m:)
    real (kind = dp_t), intent(out) :: my(lo(1)-ng_m:,lo(2)-ng_m:)
    real (kind = dp_t), intent(out) ::  s(lo(1)-ng_s:,lo(2)-ng_s:,:)  
    real (kind = dp_t), intent(in ) :: dx(:),time

    ! local
    real(kind=dp_t) :: c_face, rho_face
    real(kind=dp_t) :: x,y,freq,pfreq,hx,hy,xp,yp,xm,ym,xx,yy
    real(kind=dp_t) :: sinxt,cosxt,sinyt,cosyt,pressure,ww,dw
    real(kind=dp_t) :: ut,vt,pt,rhot,yc,xc,vcst,ucst,r,theta
    real(kind=dp_t) :: ufac,pfac,sfac,qc(1:3),dc(1:3),qt,length
    real(kind=dp_t) :: temp1,temp2,F,D
    integer         :: i,j,ii,jj
    
    ! mcai local
    real(kind=dp_t) :: one_third_domain1, one_third_domain2, y1, y2

    ! arrays to do Simpson's rule
    qc(1)=1.0d0/6.0d0
    qc(2)=4.0d0/6.0d0
    qc(3)=1.0d0/6.0d0
    dc(1)=-0.5d0
    dc(2)=0.0d0
    dc(3)=0.5d0

    hx = dx(1)
    hy = dx(2)

    ! Provide default values for momentum (zero velocity):  
    mx = 0
    my = 0
      
    select case ( abs(prob_type) )
    case (1, 100, 21)

       ! Equilibrium setup with linear concentration profile along y

       do j=lo(2),hi(2)
          y = prob_lo(2) + dx(2) * (dble(j) + half)
          do i=lo(1),hi(1)
             x = prob_lo(1) + dx(1) * (dble(i) + half)             
             call linear_profile(s(i,j,2), s(i,j,1), y)
             s(i,j,2) = s(i,j,1)*s(i,j,2) ! rho1=rho*c
          enddo
       enddo

       do j=lo(2),hi(2)
          y = dx(2) * (dble(j) + half)
          do i=lo(1),hi(1)+1
             x = dx(1) * (dble(i) )               
             call linear_profile(c_face, rho_face, y)
             mx(i,j) = rho_face*vel_init(1)
          enddo
       enddo

       do j=lo(2),hi(2)+1
          y = dx(2) * (dble(j))
          do i=lo(1),hi(1)
             x = dx(1) * (dble(i) + half)
             call linear_profile(c_face, rho_face, y)
             my(i,j) = rho_face*vel_init(2)
          enddo
       enddo

    case (2)

       ! traveling wave with exact solution

       freq  = 2.d0*M_PI
       pfreq = 4.d0*M_PI

       ucst = 0.75d0
       vcst = 0.75d0
       ufac = dexp(-2.0d0*freq*freq*time*visc_coef)
       pfac = dexp(-4.0d0*freq*freq*time*visc_coef)/64.0d0

       do j=lo(2),hi(2)
          y = dx(2) * (dble(j) + half) - vcst*time
          yp = y+0.5d0*hy
          ym = y-0.5d0*hy
          do i=lo(1),hi(1)
             x = dx(1) * (dble(i) + half) - ucst*time

             xp = x+0.5d0*hx
             xm = x-0.5d0*hx

             cosxt =  (sin(freq*xp)-sin(freq*xm))/(freq*hx)
             sinxt = -(cos(freq*xp)-cos(freq*xm))/(freq*hx)
             cosyt =  (sin(freq*yp)-sin(freq*ym))/(freq*hy)
             sinyt = -(cos(freq*yp)-cos(freq*ym))/(freq*hy)

             s(i,j,2) = ucst + ufac * 0.25d0*cosxt*sinyt
          enddo
       enddo

       do j=lo(2),hi(2)
          y = dx(2) * (dble(j) + half) - vcst*time
          yp = y+0.5d0*hy
          ym = y-0.5d0*hy
          do i=lo(1),hi(1)+1
             x = dx(1) * (dble(i) ) - ucst*time

             xp = x+0.5d0*hx
             xm = x-0.5d0*hx

             cosxt =  (sin(freq*xp)-sin(freq*xm))/(freq*hx)
             sinxt = -(cos(freq*xp)-cos(freq*xm))/(freq*hx)
             cosyt =  (sin(freq*yp)-sin(freq*ym))/(freq*hy)
             sinyt = -(cos(freq*yp)-cos(freq*ym))/(freq*hy)

             mx(i,j) = ucst + ufac * 0.25d0*cosxt*sinyt
          enddo
       enddo

       do j=lo(2),hi(2)+1
          y = dx(2) * (dble(j)) - vcst*time
          yp = y+0.5d0*hy
          ym = y-0.5d0*hy
          do i=lo(1),hi(1)
             x = dx(1) * (dble(i) + half) - ucst*time

             xp = x+0.5d0*hx
             xm = x-0.5d0*hx

             cosxt =  (sin(freq*xp)-sin(freq*xm))/(freq*hx)
             sinxt = -(cos(freq*xp)-cos(freq*xm))/(freq*hx)
             cosyt =  (sin(freq*yp)-sin(freq*ym))/(freq*hy)
             sinyt = -(cos(freq*yp)-cos(freq*ym))/(freq*hy)

             my(i,j) = vcst - ufac * 0.25d0*sinxt*cosyt
          enddo
       enddo

       ! We no longer maintain pressure, but here it is for the record:
       if(.false.) then
       do j=lo(2),hi(2)
          y = dx(2) * (dble(j) + half) - vcst*time
          yp = y+0.5d0*hy
          ym = y-0.5d0*hy
          do i=lo(1),hi(1)
             x = dx(1) * (dble(i) + half) - ucst*time
             xp = x+0.5d0*hx
             xm = x-0.5d0*hx
             pressure =-pfac*((sin(pfreq*xp)-sin(pfreq*xm))/(pfreq*hx) &
                       +(sin(pfreq*yp)-sin(pfreq*ym))/(pfreq*hy))
          enddo
       enddo
       end if

       s(:,:,      1) = 1.d0

    case (3)         

       ! Advection diffusion on scalar

       mx = 1.d0
       my = 1.d0
       s(:,:,1) = 1.d0

       sfac = exp(-two*two*M_PI*M_PI*time*diff_coef)
       ucst = mx(lo(1),lo(2))
       vcst = my(lo(1),lo(2))

       ! We define only on the valid region here, and let the periodic filling routines
       !   take care of the ghost cells
       do j=lo(2),hi(2)
          y = dx(2) * (dble(j) + half) - vcst*time
          do i=lo(1),hi(1)
             x = dx(1) * (dble(i) + half) - ucst*time
             s(i,j,2) = sfac * (cos(two*M_PI*(x+dx(1)*0.5))-cos(two*M_PI*(x-dx(1)*0.5)))/dx(1) &
                      + sfac * (cos(two*M_PI*(y+dx(2)*0.5))-cos(two*M_PI*(y-dx(2)*0.5)))/dx(2)
          enddo
       enddo

    case (4)

       ! Variable density traveling vortex

       ucst = 1.0d0
       vcst = 1.0d0

       xc = 0.5d0+ucst*time
       yc = 0.5d0+vcst*time

       s(:,:,1) = 0.d0
       do j = lo(2),hi(2)
          y = dx(2) * (dble(j) + half) - yc
          do i = lo(1),hi(1)+1
             x = dx(1) * (dble(i)) - xc
             ! Simpson's rule
             do jj = 1,3
                yy = y + dc(jj)*hy
                do ii = 1,3
                   xx = x + dc(ii)*hx
                   r=dsqrt(xx**2+yy**2)
                   theta = datan2(yy,xx)
                   call Gresho_vortex(ut,vt,pt,rhot,r,theta)             
                   qt = qc(ii)*qc(jj)
                   mx(i,j) = mx(i,j)+qt*ut
                enddo
             enddo
          enddo
       enddo
       do j = lo(2),hi(2)+1
          y = dx(2) * (dble(j)) - yc
          do i = lo(1),hi(1)
             x = dx(1) * (dble(i) + half) - xc
             ! Simpson's rule
             do jj = 1,3
                yy = y + dc(jj)*hy
                do ii = 1,3
                   xx = x + dc(ii)*hx
                   r=dsqrt(xx**2+yy**2)
                   theta = datan2(yy,xx)
                   call Gresho_vortex(ut,vt,pt,rhot,r,theta)             
                   qt = qc(ii)*qc(jj)
                   my(i,j) = my(i,j)+qt*vt
                enddo
             enddo
          enddo
       enddo
       do j = lo(2),hi(2)
          y = dx(2) * (dble(j) + half) - yc
          do i = lo(1),hi(1)
             x = dx(1) * (dble(i) + half) - xc
             ! Simpson's rule
             do jj = 1,3
                yy = y + dc(jj)*hy
                do ii = 1,3
                   xx = x + dc(ii)*hx
                   r=dsqrt(xx**2+yy**2)
                   theta = datan2(yy,xx)
                   call Gresho_vortex(ut,vt,pt,rhot,r,theta)             
                   qt = qc(ii)*qc(jj)
                   s(i,j,1) = s(i,j,1)+qt*rhot
                enddo
             enddo
          enddo
       enddo

       s(:,:,2:nscal) = 0.d0

    case (5)

       ! Constant density traveling vortex

       ucst = 1.0d0
       vcst = 0.0d0

       xc = 0.5d0+ucst*time
       yc = 0.5d0+vcst*time

       s(:,:,1) = 1.d0
       s(:,:,2:nscal) = 0.d0

       mx = ucst
       my = vcst
       do j = lo(2),hi(2)
          y = dx(2) * (dble(j) + half) - yc
          do i = lo(1),hi(1)+1
             x = dx(1) * (dble(i) ) - xc
             ! Simpson's rule
             do jj = 1,3
                yy = y + dc(jj)*hy
                do ii = 1,3
                   xx = x + dc(ii)*hx
                   r=dsqrt(xx**2+yy**2)
                   theta = datan2(yy,xx)
                   call Gresho_constant_vortex(ut,vt,pt,r,theta)             
                   qt = qc(ii)*qc(jj)
                   mx(i,j) = mx(i,j)+qt*ut
                enddo
             enddo
          enddo
       enddo
       do j = lo(2),hi(2)+1
          y = dx(2) * (dble(j) ) - yc
          do i = lo(1),hi(1)
             x = dx(1) * (dble(i) + half) - xc
             ! Simpson's rule
             do jj = 1,3
                yy = y + dc(jj)*hy
                do ii = 1,3
                   xx = x + dc(ii)*hx
                   r=dsqrt(xx**2+yy**2)
                   theta = datan2(yy,xx)
                   call Gresho_constant_vortex(ut,vt,pt,r,theta)             
                   qt = qc(ii)*qc(jj)
                   my(i,j) = my(i,j)+qt*vt
                enddo
             enddo
          enddo
       enddo

    case (6)

       ! mixing (starting with half-half interface)
       do j=lo(2),hi(2)
          y = prob_lo(2) + dx(2) * (dble(j) + half) - 0.5d0*(prob_lo(2)+prob_hi(2))
          do i=lo(1),hi(1)

             ! tanh smoothing
             s(i,j,2) = 0.5d0*(c_init(2)-c_init(1))*(1.d0 + tanh(y/(smoothing_width*dx(2)))) + c_init(1)

             s(i,j,1) = 1.0d0/(s(i,j,2)/rhobar(1)+(1.0d0-s(i,j,2))/rhobar(2))
             s(i,j,2) = s(i,j,1)*s(i,j,2)
          enddo
       enddo

    case (7)

       ! cylinder
       do j=lo(2),hi(2)
          y = prob_lo(2) + dx(2)*(dble(j) + half)
          do i=lo(1),hi(1)
             x = prob_lo(1) + dx(1)*(dble(i) + half)

             r = sqrt( (x - 0.5*(prob_lo(1)+prob_hi(1)) )**2 + (y - 0.5*(prob_lo(2)+prob_hi(2)) )**2 ) - &
                  0.25d0*(prob_lo(1)+prob_hi(1))

             ! tanh smoothing
             s(i,j,2) = 0.5d0*(c_init(2)-c_init(1))*(1.d0 + tanh(r/(smoothing_width*dx(2)))) + c_init(1)

             s(i,j,1) = 1.0d0/(s(i,j,2)/rhobar(1)+(1.0d0-s(i,j,2))/rhobar(2))
             s(i,j,2) = s(i,j,1)*s(i,j,2)
          enddo
       enddo

    case (8)   

       !mcai: bilayer interface (stripe)
       one_third_domain1=2.0d0/3.0d0*prob_lo(2)+1.0d0/3.0d0*prob_hi(2)
       one_third_domain2=1.0d0/3.0d0*prob_lo(2)+2.0d0/3.0d0*prob_hi(2)

       do j=lo(2),hi(2)
          y1 =(prob_lo(2) + dx(2) * (dble(j) + half) - one_third_domain1)
          y2 =(prob_lo(2) + dx(2) * (dble(j) + half) - one_third_domain2)
        
          do i=lo(1),hi(1)
             ! tanh smoothing
             if(abs(smoothing_width)>epsilon(1.0d0)) then
                s(i,j,2) = c_init(1)+ 0.5d0*(c_init(2)-c_init(1))*&
                   (tanh(y1/(smoothing_width*dx(2))) - tanh(y2/(smoothing_width*dx(2))))
                s(i,j,1) = 1.0d0/(s(i,j,2)/rhobar(1)+(1.0d0-s(i,j,2))/rhobar(2))
                s(i,j,2) = s(i,j,1)*s(i,j,2)
             else ! Try to initialize exactly as we do in the HDMD simulations, with finite-volume averaging of sharp interface
                if((y1<-0.5d0*dx(2)).or.(y2>0.5d0*dx(2))) then
                   s(i,j,2) = 0
                   s(i,j,1) = rhobar(2)
                else if((y1>0.5d0*dx(2)).and.(y2<-0.5d0*dx(2))) then
                   s(i,j,2) = rhobar(1)
                   s(i,j,1) = rhobar(1)
                else if(y1 <= 0.5d0*dx(2)) then
                   s(i,j,2) = (max(0.0d0,min(0.5d0+y1/dx(2),1.0d0)))*rhobar(1)
                   s(i,j,1) = s(i,j,2) + (1.0d0-max(0.0d0,min(0.5d0+y1/dx(2),1.0d0)))*rhobar(2)
                else 
                   s(i,j,2) = (1.0d0-max(0.0d0,min(0.5d0+y2/dx(2),1.0d0)))*rhobar(1)
                   s(i,j,1) = s(i,j,2) + (max(0.0d0,min(0.5d0+y2/dx(2),1.0d0)))*rhobar(2)
                end if   
             end if  
          enddo
       enddo         

    case (9)

       ! sawtooth initialization 
       length = prob_hi(2)-prob_lo(2)

       do j=lo(2),hi(2)
          y = dx(2) * (dble(j) + half)
          do i=lo(1),hi(1)

             if (y .le. length/2.0d0) then
                s(i,j,2) = c_init(1)+(c_init(2)-c_init(1))*2.0d0/length*y
             else
                s(i,j,2) = c_init(2)-(c_init(2)-c_init(1))*2.0d0/length*(y-length/2.0d0) 
             end if

             s(i,j,1) = 1.0d0/(s(i,j,2)/rhobar(1)+(1.0d0-s(i,j,2))/rhobar(2))
             s(i,j,2) = s(i,j,1)*s(i,j,2)
          enddo
       enddo

    case (10)

       ! forced flow with BCs 

       freq  = 2.d0*M_PI
       pfreq = 4.d0*M_PI
       ww = 1.0d0+sin(freq*time*time)
       dw = pfreq*time*cos(freq*time*time) 

       do j=lo(2),hi(2)
          y = dx(2) * (dble(j) + half)
          do i=lo(1),hi(1)+1
             x = dx(1) * (dble(i) )
               
             cosxt = cos(freq*(x-ww)) 
             yp = 3.0d0*y*y-2.0d0*y

             mx(i,j) = cosxt*yp 
          enddo
       enddo

       do j=lo(2),hi(2)+1
          y = dx(2) * (dble(j))
          do i=lo(1),hi(1)
             x = dx(1) * (dble(i) + half)

             sinxt = sin(freq*(x-ww)) 
             ym = y*y*(y-1.0d0)

             my(i,j) = freq*sinxt*ym 
          enddo
       enddo

       ! We no longer maintain pressure, but here it is for the record:
       if(.false.) then
       do j=lo(2),hi(2)
          y = dx(2) * (dble(j) + half)
          do i=lo(1),hi(1)
             x = dx(1) * (dble(i) + half)

             cosxt = cos(freq*(x-ww)) 
             yp = 3.0d0*y*y-2.0d0*y
             sinxt = sin(freq*(x-ww)) 
             ym = y*y*(y-1.0d0)

             pressure = -dw/freq*sinxt*(sin(freq*y)-freq*y+M_PI) &
                       -visc_coef*cosxt*(-2.0d0*sin(freq*y)+freq*y-M_PI) 
          enddo
       enddo
       end if

       s(:,:,1) = 1.d0

    case (11)

       ! incompressible modes with BCs       

       freq  = mode_coefs(1)*2.0*M_PI ! frequency in the periodic, "mod(prob_dir+1,DIM)" direction
       
       select case(mode_coefs(1))
       case(1)
          select case(mode_coefs(2))
          ! frequency in non-periodic, "prob_dir" direction
          case(1)
             pfreq  = 1.8583684196603974620d0 
          case(2)
             pfreq  = 5.4245857755826005234d0
          case(3)
             pfreq  = 8.8049901634565292604d0
          case default
             call bl_error('ERROR: Only mode_coefs(2)<=3 supported for prob_type=11 in exact_2d')
          end select
       case(2)
          select case(mode_coefs(2))
          ! frequency in non-periodic, "prob_dir" direction
          case(1)
             pfreq  = 2.1789521056066303662d0
          case(2)
             pfreq  = 5.7874274827134140679d0
          case(3)
             pfreq  = 9.0932821164949254163d0
          case default
             call bl_error('ERROR: Only mode_coefs(2)<=3 supported for prob_type=11 in exact_2d')
          end select   
       case default
          call bl_error('ERROR: Only mode_coefs(1)<=2 supported for prob_type=11 in exact_2d')
       end select            
       
       ucst = visc_coef*(freq*freq+pfreq*pfreq) 
       ! since this is an exact solution for the time-dependant Stokes equations, 
       ! we multiply by a small prefactor so that the nonlinearity is negligable
       ufac = (1E-8)*dexp(-ucst*time)

       if (prob_dir .eq. 1) then

          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2) * (dble(j) + half)
             do i=lo(1),hi(1)+1
                x = prob_lo(1) + dx(1) * (dble(i) )
                
                mx(i,j) = ufac*sin(freq*y)/cos(pfreq)*(cos(pfreq*x)*cosh(freq) &
                     -cosh(freq*x)*cos(pfreq)) 
             enddo
          enddo

          do j=lo(2),hi(2)+1
             y = prob_lo(2) + dx(2) * (dble(j))
             do i=lo(1),hi(1)
                x = prob_lo(1) + dx(1) * (dble(i) + half)
                
                my(i,j) = ufac*cos(freq*y)/sin(pfreq)*(sin(pfreq*x)*sinh(freq) &
                     -sinh(freq*x)*sin(pfreq)) 
             enddo
          enddo

       else if (prob_dir .eq. 2) then

          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2) * (dble(j) + half)
             do i=lo(1),hi(1)+1
                x = prob_lo(1) + dx(1) * (dble(i) )
                
                mx(i,j) = ufac*cos(freq*x)/sin(pfreq)*(sin(pfreq*y)*sinh(freq) &
                     -sinh(freq*y)*sin(pfreq)) 

             enddo
          enddo

          do j=lo(2),hi(2)+1
             y = prob_lo(2) + dx(2) * (dble(j))
             do i=lo(1),hi(1)
                x = prob_lo(1) + dx(1) * (dble(i) + half)
                
                my(i,j) = ufac*sin(freq*x)/cos(pfreq)*(cos(pfreq*y)*cosh(freq) &
                     -cosh(freq*y)*cos(pfreq)) 
             enddo
          enddo

       else
          call bl_error('exact_solutions.f90: invalid prob_dir for prob_type=11')
       end if

       ! We no longer maintain pressure, but here it is for the record:
       if(.false.) then
       do j=lo(2),hi(2)
          y = dx(2) * (dble(j) + half)
          do i=lo(1),hi(1)
             x = dx(1) * (dble(i) + half)
             ! this may be wrong--check pfac!?!
             pressure = -ufac*sinh(freq*x)*sin(freq*y)/freq
          enddo
       enddo
       end if

       s(:,:,      1) = 1.d0



    case (12)

       ! lo density spherical bubble
       do j=lo(2),hi(2)
          y = prob_lo(2) + dx(2) * (dble(j) + half) - 0.5d0*(prob_lo(2)+prob_hi(2))
          do i=lo(1),hi(1)
             x = prob_lo(1) + dx(1) * (dble(i) + half) - 0.5d0*(prob_lo(1)+prob_hi(1))

             r = sqrt (x**2 + y**2)

             ! tanh smoothing
             s(i,j,2) = c_init(1) + 0.5d0*(c_init(2)-c_init(1))* &
                  (1.d0 + tanh((r-2.5d0*smoothing_width*dx(1))/(smoothing_width*dx(1))))

             s(i,j,1) = 1.0d0/(s(i,j,2)/rhobar(1)+(1.0d0-s(i,j,2))/rhobar(2))
             s(i,j,2) = s(i,j,1)*s(i,j,2)
          enddo
       enddo

    case (13)

       ! constant c gradient, convert to rho and rho*c later
       do j=lo(2),hi(2)
          y = prob_lo(2) + dx(2) * (dble(j) + half)
          do i=lo(1),hi(1)
             x = prob_lo(1) + dx(1) * (dble(i) + half)
             s(i,j,2) = c_init(1) + &
                  (y-prob_lo(2))/(prob_hi(2)-prob_lo(2))*(c_init(2)-c_init(1))
          enddo
       enddo

       ! we need a profile such that rho grad c = constant so the deterministic problem will 'sit still'
       ! rho grad c = CONST
       ! rho = (c/rho1bar + (1-c)/rho2bar)^-1
       ! dc/dz = CONST
       ! integrating, we have CONST = int_{c_lo}^{c_hi} rho dc / (prob_hi(2)-prob_lo(2))
       ! Then, int_{c_lo}^{c(z)} = CONST*z
       ! gives c = (rho1bar - exp((1/rho1bar - 1/rho2bar)*(CONST*z-I(c_lo))) ) / (rho1bar - rho2bar)

!       temp1 = rhobar(1)*rhobar(2)*log( rhobar(1)*(1.d0-c_init(1)) + rhobar(2)*c_init(1) ) / (rhobar(2)-rhobar(1))
!       temp2 = rhobar(1)*rhobar(2)*log( rhobar(1)*(1.d0-c_init(2)) + rhobar(2)*c_init(2) ) / (rhobar(2)-rhobar(1))
!       F = (temp2 - temp1) / (prob_hi(2)-prob_lo(2))
!
!       ! compute c, will multiply by rho later
!       do j=lo(2),hi(2)
!          y = prob_lo(2) + dx(2) * (dble(j) + half)
!          do i=lo(1),hi(1)
!             x = prob_lo(1) + dx(1) * (dble(i) + half)
!
!             D = F*(y-prob_lo(2)) + temp1
!             s(i,j,2) = ( rhobar(1) - exp((1.d0/rhobar(1)-1.d0/rhobar(2))*D) ) / (rhobar(1)-rhobar(2))
!          enddo
!       enddo
             
       ! convert c to rho and rho*c
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             s(i,j,1) = 1.0d0/(s(i,j,2)/rhobar(1)+(1.0d0-s(i,j,2))/rhobar(2))
             s(i,j,2) = s(i,j,1)*s(i,j,2)
          enddo
       enddo
       
    case default
       call bl_error('ERROR: bad or unimplemented choice for prob_type in exact_2d')
    end select

  end subroutine exact_2d

  ! Note: Thi routine *increments* s_update: Do not reset it to zero! 
  subroutine exact_s_force_2d(s_update,s,ng_u,ng_s,lo,hi,dx,time)
 
    integer        , intent(in   ) :: lo(:),hi(:),ng_u,ng_s
    real(kind=dp_t), intent(inout) :: s_update(lo(1)-ng_u:,lo(2)-ng_u:,:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),time
    
    ! local
    integer :: i,j

    real(kind=dp_t) :: tm, t, tp, y, x
    real(kind=dp_t) :: pfac, ucst, vcst, freq, pfreq, length 

    select case ( abs(prob_type) )
    case (2) ! Taylor (traveling wave) vortices

       ucst = 0.75d0
       vcst = 0.75d0
       freq  = 2.d0*M_PI
       pfreq  = 4.d0*M_PI
       pfac = dexp(-4.0d0*freq*freq*time*diff_coef)/64.0d0

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             x = dx(1) * (dble(i) + half) - ucst*time
             s_update(i,j,2) = s_update(i,j,2) - &
                       pfac*pfreq*sin(pfreq*x)
          enddo
       enddo
    
    case (9) ! Sawtooth (triangle)   

       length = prob_hi(2)-prob_lo(2)

       do j = lo(2),hi(2)
          y = dx(2) * (dble(j) + half)
          do i = lo(1),hi(1)
             call triang(y,t,c_init,dx,length)
             call triang(y-dx(2),tm,c_init,dx,length)
             call triang(y+dx(2),tp,c_init,dx,length)
             s_update(i,j,2) = s_update(i,j,2) - 0.5*triangle_coeff/(dx(2)*dx(2)) &
                  *((s(i,j,1)+s(i,j+1,1))*(tp-t)&     
                  -(s(i,j-1,1)+s(i,j,1))*(t-tm))  
          end do
       end do

    end select

  end subroutine exact_s_force_2d

  ! Note: This routine *increments* m_update: Do not reset it to zero! 
  subroutine exact_m_force_2d(m_updatex,m_updatey,ng_u,lo,hi,dx,time)
 
    integer        , intent(in   ) :: lo(:),hi(:),ng_u
    real(kind=dp_t), intent(inout) :: m_updatex(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t), intent(inout) :: m_updatey(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t), intent(in   ) :: dx(:),time
    
    ! local
    integer :: i,j

    real(kind=dp_t) :: y, x
    real(kind=dp_t) :: ww, dw, cosxt, sinxt
    real(kind=dp_t) :: yp, ym, ut, vt, freq, pfreq 

    select case ( abs(prob_type) )

    case (10) ! forced flow with BCs 

       freq  = 2.d0*M_PI
       pfreq = 4.d0*M_PI
       ww = 1.0d0+sin(freq*time*time)
       dw = pfreq*time*cos(freq*time*time) 

       do j=lo(2),hi(2)
          y = dx(2) * (dble(j) + half)
          do i=lo(1),hi(1)+1
             x = dx(1) * (dble(i) )
               
             cosxt = cos(freq*(x-ww)) 
             yp = 3.0d0*y*y-2.0d0*y
             sinxt = sin(freq*(x-ww)) 
             ym = y*y*(y-1.0d0)

             ut = cosxt*yp 
             vt = freq*sinxt*ym 

             m_updatex(i,j) = m_updatex(i,j) &
                        +freq*dw*sinxt*yp &
                        +ut*(-freq*sinxt*yp) &
                        +vt*cosxt*(6.0d0*y-2.0d0) &
                        -dw*cosxt*(sin(freq*y)-freq*y+M_PI) & 
                        +visc_coef*freq*sinxt*(-2.0d0*sin(freq*y)+freq*y-M_PI) &
                        -visc_coef*(-pfreq*M_PI*cosxt*yp &
                        +6.0d0*cosxt) 

          enddo
       enddo

       do j=lo(2),hi(2)+1
          y = dx(2) * (dble(j))
          do i=lo(1),hi(1)
             x = dx(1) * (dble(i) + half)
               
             cosxt = cos(freq*(x-ww)) 
             yp = 3.0d0*y*y-2.0d0*y
             sinxt = sin(freq*(x-ww)) 
             ym = y*y*(y-1.0d0)

             ut = cosxt*yp 
             vt = freq*sinxt*ym 

             m_updatey(i,j) = m_updatey(i,j) &
                        -pfreq*M_PI*dw*cosxt*ym &
                        +ut*(pfreq*M_PI*cosxt*ym) &
                        +vt*freq*sinxt*y*(3.0d0*y-2.0d0) &
                        -dw*sinxt*(cos(freq*y)-1.0d0) & 
                        -visc_coef*freq*cosxt*(-2.0d0*cos(freq*y)+1.0d0) &
                        -visc_coef*(-pfreq*freq*M_PI*sinxt*ym &
                        +freq*sinxt*(6.0d0*y-2.0d0)) 

          enddo
       enddo

    end select

  end subroutine exact_m_force_2d

  subroutine triang(a,b,c_init,dx,length)

    real (kind = dp_t), intent(in   ) :: dx(:),a,c_init(2),length
    real (kind = dp_t), intent(inout) :: b

    ! Local variables
    real(kind=dp_t) :: x,l

    l = length 
    x = modulo(a,l)

    !if (x .le. l/4.0d0) then
    !   b = 0.5d0-2.0d0/l*x
    !else if (x .gt. 0.75d0*l) then
    !   b = 1.0d0-2.0d0/l*(x-0.75d0*l)
    !else
    !   b = 2.0d0/l*(x-0.25d0*l)
    !end if

    if (x .le. l/2.0d0) then
       b = c_init(1)+(c_init(2)-c_init(1))*2.0d0/l*x
    else
       b = c_init(2)-(c_init(2)-c_init(1))*2.0d0/l*(x-l/2.0d0)  
    end if

  end subroutine triang

  ! 
  ! ---------------------------------------------------------------------------
  ! 

  subroutine Gresho_vortex(u,v,p,rho,r,theta)
    real(kind=dp_t), intent(out) :: u,v,p,rho
    real(kind=dp_t), intent(in) :: r,theta  ! radius  and theta
    !  Return the starting value of 0-Mach vortex problem

    real(kind=dp_t) :: R1,R2,R6,Rval  !  Scaled radius squared
    real(kind=dp_t) :: nm(25),dm(25),coeff(25),pcst
    integer ::  k

    data nm/ 1,  6,15, 74,57,174,269,450,153,1564,510,204,1473,1014,1053,558,783,54,38,222,609,184,9, 12, 1/
    data dm/72,-35,17,-33,32, 31,-15, 29,  8, -27, 13,  5,- 16,  23,  22, -7, 20,19,-9,-17, 32,-15,2,-13,12/

    !  Compute the constant in pressure polynomial ,i.e. p(R)
    pcst = 0.0d0
    do k = 1,25
       coeff(k) = nm(k) / dm(k)
       pcst = pcst + coeff(k)
    end do

    R1=r/(0.4d0)
    R2=R1*R1
    R6=R2*R2*R2

    if (R1 > 1.0d0) then

       rho = 0.5d0
       u = 1.0d0*rho
       v = 0
       p = 0

    else

       rho = 0.5d0 + 0.5d0*(1.d0-R2)**6
       Rval = ((1.0d0-R1)**6)*R6

       u = -1024.0d0*Rval*sin(theta)
       v =  1024.0d0*Rval*cos(theta)

       !  Return the pressure for the variable density Gresho vortex  problem
       !  (See Eq. 96 in KKM)

       !  Horners rule

       p = coeff(1) * R1 
       do k = 2,24
          p = (p + coeff(k)) * R1
       end do
       p = p + coeff(25)

       p = p*R6*R6
       p = p-pcst
       p = p*(1024.0d0)*(1024.0d0)

       u=(1.0d0+u)*rho
       v=v*rho
    end if

  end subroutine Gresho_vortex

  ! 
  ! ---------------------------------------------------------------------------
  ! 

  subroutine Gresho_constant_vortex(u,v,p,r,theta)
    real(kind=dp_t), intent(out) :: u,v,p
    real(kind=dp_t), intent(in) :: r,theta  ! radius  and theta
    !  Return the starting value of 0-Mach vortex problem

    real(kind=dp_t) :: R1,R2,R6,Rval  !  Scaled radius squared
    real(kind=dp_t) :: nm(13),dm(13),coeff(13),pcst
    integer ::  k

    data nm/1.d0, -12.D0, 3.d0,-220.d0, 99.d0,-792.d0, 154.d0, -792.d0, 495.d0, -44.d0, &
         33.d0, -12.d0, 1.d0/ 
    data dm/24.d0, 23.d0, 1.d0,  21.d0,  4.d0,  19.d0, 3.d0,     17.d0,  16.d0,   3.d0, &
         7.d0, 13.d0, 12.d0/ 

    !  Compute the constant in pressure polynomial ,i.e. p(R)
    pcst = 0.0d0
    do k = 1,13
       coeff(k) = nm(k) / dm(k)
       pcst = pcst + coeff(k)
    end do

    R1=r/(0.4d0)
    R2=R1*R1
    R6=R2*R2*R2

    if (R1 > 1.0d0) then

       u = 0
       v = 0

    else

       Rval = ((1.0d0-R1)**6)*R6

       u = -1024.0d0*Rval*sin(theta)
       v =  1024.0d0*Rval*cos(theta)

       !  Return the pressure for the constant density Gresho vortex problem

       !  Horners rule
       p = coeff(1) * R1 
       do k = 2,12
          p = (p + coeff(k)) * R1
       end do
       p = p + coeff(13)

       p = p*R6*R6
       p = p-pcst
       p = p*(1024.0d0)*(1024.0d0)

    end if

  end subroutine Gresho_constant_vortex

   subroutine linear_profile(concentration, density, y)
      real(kind=dp_t), intent(out) :: concentration, density
      real(kind=dp_t), intent(in) :: y

      ! Concentration c is linear along y going from c_init(1) to c_init(2) from boundary-to-boundary
      concentration = (c_init(2)-c_init(1))*((y-prob_lo(2))/(prob_hi(2)-prob_lo(2))) + c_init(1)
      ! Density follows from the EOS:      
      density = 1.0d0/(concentration/rhobar(1)+(1.0d0-concentration)/rhobar(2))
    end subroutine  

  ! 
  ! ---------------------------------------------------------------------------
  ! 3D forcing
  ! ---------------------------------------------------------------------------
  ! 

  subroutine exact_3d(mx,my,mz,s,lo,hi,ng_m,ng_s,dx,time)

    integer, intent(in) :: lo(:), hi(:), ng_m, ng_s
    real (kind = dp_t), intent(out) :: mx(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real (kind = dp_t), intent(out) :: my(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real (kind = dp_t), intent(out) :: mz(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real (kind = dp_t), intent(out) ::  s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)  
    real (kind = dp_t), intent(in ) :: dx(:),time

    ! local
    real(kind=dp_t) :: c_face, rho_face
    real(kind=dp_t) :: x,y,z,ufac,vfac,pfac,hx,hy,hz
    real(kind=dp_t) :: freq,pfreq,freqx,freqy
    real(kind=dp_t) :: xp,yp,zp,xm,ym,zm,xg,yg,zg
    real(kind=dp_t) :: ucst,vcst,wcst
    real(kind=dp_t) :: sinxt,cosxt,sinyt,cosyt,sinzt,coszt
    real(kind=dp_t) :: sfac,r
    integer         :: i,j,k,ii,jj,kk

    ! mcai local
    real(kind=dp_t) :: one_third_domain1, one_third_domain2, y1, y2

    hx = dx(1)
    hy = dx(2)
    hz = dx(3)
    freq = two *M_PI

    ! Provide default values of momentum (zero velocity)
    mx = 0
    my = 0
    mz = 0

    select case ( abs(prob_type) )
    case (1, 100, 21)

       ! Equilibrium setup with constant velocity and linear concentration profile along y
       
       do k=lo(3),hi(3)
          z = prob_lo(3) + dx(3) * (dble(k) + half)
          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2) * (dble(j) + half) 
             do i=lo(1),hi(1)
                x = prob_lo(1) + dx(1) * (dble(i) + half)                
                ! Concentration c is linear along y going from c_init(1) to c_init(2) from boundary-to-boundary
                call linear_profile(s(i,j,k,2), s(i,j,k,1), y)
                s(i,j,k,2) = s(i,j,k,1)*s(i,j,k,2)
             enddo
          enddo
       end do

       ! constant velocity
       do k=lo(3),hi(3)
          z = dx(3) * (dble(k) + half)
          do j=lo(2),hi(2)
             y = dx(2) * (dble(j) + half)
             do i=lo(1),hi(1)+1
                x = dx(1) * (dble(i))
                call linear_profile(c_face, rho_face, y)
                mx(i,j,k) = rho_face*vel_init(1)
             enddo
          enddo
       enddo

       do k=lo(3),hi(3)
          z = dx(3) * (dble(k) + half)
          do j=lo(2),hi(2)+1
             y = dx(2) * (dble(j))
             do i=lo(1),hi(1)
                x = dx(1) * (dble(i) + half)
                call linear_profile(c_face, rho_face, y)
                my(i,j,k) = rho_face*vel_init(2)
             enddo
          enddo
       enddo

       do k=lo(3),hi(3)+1
          z = dx(3) * (dble(k))
          do j=lo(2),hi(2)
             y = dx(2) * (dble(j) + half)
             do i=lo(1),hi(1)
                x = dx(1) * (dble(i) + half)
                call linear_profile(c_face, rho_face, y)
                mz(i,j,k) = rho_face*vel_init(3)
             enddo
          enddo
       enddo

    case (2)

       ! ABC flow 

       freq  = 2.d0*M_PI

       ufac = dexp(-freq*freq*visc_coef*time)
       pfac = dexp(-2.0d0*freq*freq*visc_coef*time)

       hx = ABC_coefs(1)
       hy = ABC_coefs(2)
       hz = ABC_coefs(3)

       do k=lo(3),hi(3)
          z = dx(3) * (dble(k) + half)
          do j=lo(2),hi(2)
             y = dx(2) * (dble(j) + half)
             do i=lo(1),hi(1)
                x = dx(1) * (dble(i) + half)

                s(i,j,k,2) = 1.0d0 + ufac*(hz*cos(freq*(y-time))+hx*sin(freq*(z-time)))

             enddo
          enddo
       enddo

       do k=lo(3),hi(3)
          z = dx(3) * (dble(k) + half)
          do j=lo(2),hi(2)
             y = dx(2) * (dble(j) + half)
             do i=lo(1),hi(1)+1
                x = dx(1) * (dble(i))

                mx(i,j,k) = 1.0d0 + ufac*(hz*cos(freq*(y-time))+hx*sin(freq*(z-time)))

             enddo
          enddo
       enddo

       do k=lo(3),hi(3)
          z = dx(3) * (dble(k) + half)
          do j=lo(2),hi(2)+1
             y = dx(2) * (dble(j))
             do i=lo(1),hi(1)
                x = dx(1) * (dble(i) + half)

                my(i,j,k) = 1.0d0 + ufac*(hy*sin(freq*(x-time))+hx*cos(freq*(z-time)))

             enddo
          enddo
       enddo

       do k=lo(3),hi(3)+1
          z = dx(3) * (dble(k))
          do j=lo(2),hi(2)
             y = dx(2) * (dble(j) + half)
             do i=lo(1),hi(1)
                x = dx(1) * (dble(i) + half)

                mz(i,j,k) = 1.0d0 + ufac*(hy*cos(freq*(x-time))+hz*sin(freq*(y-time)))

             enddo
          enddo
       enddo

       s(:,:,:,      1) = 1.d0

    case (3)         

       ! Advection diffusion on scalar

       mx = 1.d0
       my = 1.d0
       mz = 1.d0
       s(:,:,:,1) = 1.d0

       sfac = exp(-two*two*M_PI*M_PI*time*diff_coef)
       ucst = mx(lo(1),lo(2),lo(3))
       vcst = my(lo(1),lo(2),lo(3))
       wcst = mz(lo(1),lo(2),lo(3))

       ! We define only on the valid region here, and let the periodic filling routines
       !   take care of the ghost cells
       do k=lo(3),hi(3)
          z = dx(3) * (dble(k) + half) - wcst*time
          do j=lo(2),hi(2)
             y = dx(2) * (dble(j) + half) - vcst*time 
             do i=lo(1),hi(1)
                x = dx(1) * (dble(i) + half) - ucst*time
                s(i,j,k,2) = (cos(two*M_PI*(x+dx(1)*0.5))-cos(two*M_PI*(x-dx(1)*0.5)))/dx(1) &
                            +(cos(two*M_PI*(y+dx(2)*0.5))-cos(two*M_PI*(y-dx(2)*0.5)))/dx(2) &
                            +(cos(two*M_PI*(z+dx(3)*0.5))-cos(two*M_PI*(z-dx(3)*0.5)))/dx(3)
                s(i,j,k,2) = sfac*s(i,j,k,2)
             enddo
          enddo
       enddo

    case (6)

       ! mixing variant
       do k=lo(3),hi(3)
          z = prob_lo(3) + dx(3) * (dble(k) + half) - 0.5d0*(prob_lo(3)+prob_hi(3))
          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2) * (dble(j) + half) - 0.5d0*(prob_lo(2)+prob_hi(2))
             do i=lo(1),hi(1)

                ! tanh smoothing
                !s(i,j,k,2) = 0.5d0*(c_init(2)-c_init(1))*(1.d0 + tanh(z/(smoothing_width*dx(3)))) + c_init(1)
                s(i,j,k,2) = 0.5d0*(c_init(2)-c_init(1))*(1.d0 + tanh(y/(smoothing_width*dx(2)))) + c_init(1)

                s(i,j,k,1) = 1.0d0/(s(i,j,k,2)/rhobar(1)+(1.0d0-s(i,j,k,2))/rhobar(2))
                s(i,j,k,2) = s(i,j,k,1)*s(i,j,k,2)
             enddo
          enddo
       end do

    case (7)

       ! cylinder
       do j=lo(2),hi(2)
          y = prob_lo(2) + dx(2)*(dble(j) + half)
          do i=lo(1),hi(1)
             x = prob_lo(1) + dx(1)*(dble(i) + half)

             r = sqrt( (x - 0.5*(prob_lo(1)+prob_hi(1)) )**2 + (y - 0.5*(prob_lo(2)+prob_hi(2)) )**2 ) - &
                  0.25d0*(prob_lo(1)+prob_hi(1))

             ! tanh smoothing
             s(i,j,:,2) = 0.5d0*(c_init(2)-c_init(1))*(1.d0 + tanh(r/(smoothing_width*dx(2)))) + c_init(1)

             s(i,j,:,1) = 1.0d0/(s(i,j,:,2)/rhobar(1)+(1.0d0-s(i,j,:,2))/rhobar(2))
             s(i,j,:,2) = s(i,j,:,1)*s(i,j,:,2)
          enddo
       enddo

    case (8)

       ! bilayer interface (stripe)
       one_third_domain1=2.0d0/3.0d0*prob_lo(2)+1.0d0/3.0d0*prob_hi(2)
       one_third_domain2=1.0d0/3.0d0*prob_lo(2)+2.0d0/3.0d0*prob_hi(2)

       do k=lo(3),hi(3)
          z = prob_lo(3) + dx(3) * (dble(k) + half) - 0.5d0*(prob_lo(3)+prob_hi(3))
          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2) * (dble(j) + half) - 0.5d0*(prob_lo(2)+prob_hi(2))
             y1 =(prob_lo(2) + dx(2) * (dble(j) + half) - one_third_domain1)
             y2 =(prob_lo(2) + dx(2) * (dble(j) + half) - one_third_domain2)
             do i=lo(1),hi(1)

                ! tanh smoothing
                s(i,j,k,2) = 0.5d0*(c_init(2)-c_init(1))*(tanh(y1/(smoothing_width*dx(2))) - &
                                                          tanh(y2/(smoothing_width*dx(2)))) + c_init(1)

                s(i,j,k,1) = 1.0d0/(s(i,j,k,2)/rhobar(1)+(1.0d0-s(i,j,k,2))/rhobar(2))
                s(i,j,k,2) = s(i,j,k,1)*s(i,j,k,2)
             enddo
          enddo
       end do
       
    case (11)

       ! incompressible modes with BCs
       ! Only quasi-2D flow supported for now

       freq  = mode_coefs(1)*2.0*M_PI ! frequency in the periodic, "mod(prob_dir+1,DIM)" direction

       select case(mode_coefs(1))
       case(1)
          select case(mode_coefs(2))
          ! frequency in non-periodic, "prob_dir" direction
          case(1)
             pfreq  = 1.8583684196603974620d0 
          case(2)
             pfreq  = 5.4245857755826005234d0
          case(3)
             pfreq  = 8.8049901634565292604d0
          case default
             call bl_error('ERROR: Only mode_coefs(2)<=3 supported for prob_type=11 in exact_2d')
          end select
       case(2)
          select case(mode_coefs(2))
          ! frequency in non-periodic, "prob_dir" direction
          case(1)
             pfreq  = 2.1789521056066303662d0
          case(2)
             pfreq  = 5.7874274827134140679d0
          case(3)
             pfreq  = 9.0932821164949254163d0
          case default
             call bl_error('ERROR: Only mode_coefs(2)<=3 supported for prob_type=11 in exact_2d')
          end select   
       case default
          call bl_error('ERROR: Only mode_coefs(1)<=2 supported for prob_type=11 in exact_2d')
       end select

       ucst = visc_coef * (freq*freq + pfreq*pfreq)
       ! since this is an exact solution for the time-dependant Stokes equations, 
       ! we multiply by a small prefactor so that the nonlinearity is negligable
       ufac = (1E-8)*dexp(-ucst*time)

       if (prob_dir .eq. 1) then

          do k=lo(3),hi(3)
             z = prob_lo(3) + dx(3) * (dble(k) + half)
             do j=lo(2),hi(2)
                y = prob_lo(2) + dx(2) * (dble(j) + half)
                do i=lo(1),hi(1)+1
                   x = prob_lo(1) + dx(1) * (dble(i) )
                
                   mx(i,j,k) = ufac*sin(freq*y)/cos(pfreq)*(cos(pfreq*x)*cosh(freq) &
                        -cosh(freq*x)*cos(pfreq)) 
                enddo
             enddo
          enddo

          do k=lo(3),hi(3)
             z = prob_lo(3) + dx(3) * (dble(k) + half)
             do j=lo(2),hi(2)+1
                y = prob_lo(2) + dx(2) * (dble(j))
                do i=lo(1),hi(1)
                   x = prob_lo(1) + dx(1) * (dble(i) + half)
                
                   my(i,j,k) = ufac*cos(freq*y)/sin(pfreq)*(sin(pfreq*x)*sinh(freq) &
                        -sinh(freq*x)*sin(pfreq)) 
                enddo
             enddo
          enddo

          do k=lo(3),hi(3)+1
             z = prob_lo(3) + dx(3) * (dble(k))
             do j=lo(2),hi(2)
                y = prob_lo(2) + dx(2) * (dble(j) + half)
                do i=lo(1),hi(1)
                   x = prob_lo(1) + dx(1) * (dble(i) + half)
                
                   mz(i,j,k) = 0.d0

                enddo
             enddo
          enddo

       else if (prob_dir .eq. 2) then

          do k=lo(3),hi(3)
             z = prob_lo(3) + dx(3) * (dble(k) + half)
             do j=lo(2),hi(2)
                y = prob_lo(2) + dx(2) * (dble(j) + half)
                do i=lo(1),hi(1)+1
                   x = prob_lo(1) + dx(1) * (dble(i) )
                
                   mx(i,j,k) = 0.d0

                enddo
             enddo
          enddo

          do k=lo(3),hi(3)
             z = prob_lo(3) + dx(3) * (dble(k) + half)
             do j=lo(2),hi(2)+1
                y = prob_lo(2) + dx(2) * (dble(j))
                do i=lo(1),hi(1)
                   x = prob_lo(1) + dx(1) * (dble(i) + half)

                   my(i,j,k) = ufac*sin(freq*z)/cos(pfreq)*(cos(pfreq*y)*cosh(freq) &
                        -cosh(freq*y)*cos(pfreq))

                enddo
             enddo
          enddo

          do k=lo(3),hi(3)+1
             z = prob_lo(3) + dx(3) * (dble(k))
             do j=lo(2),hi(2)
                y = prob_lo(2) + dx(2) * (dble(j) + half)
                do i=lo(1),hi(1)
                   x = prob_lo(1) + dx(1) * (dble(i) + half)
                                
                   mz(i,j,k) = ufac*cos(freq*z)/sin(pfreq)*(sin(pfreq*y)*sinh(freq) &
                        -sinh(freq*y)*sin(pfreq)) 

                enddo
             enddo
          enddo

       else if (prob_dir .eq. 3) then

          do k=lo(3),hi(3)
             z = prob_lo(3) + dx(3) * (dble(k) + half)
             do j=lo(2),hi(2)
                y = prob_lo(2) + dx(2) * (dble(j) + half)
                do i=lo(1),hi(1)+1
                   x = prob_lo(1) + dx(1) * (dble(i) )
                
                   mx(i,j,k) = ufac*cos(freq*x)/sin(pfreq)*(sin(pfreq*z)*sinh(freq) &
                        -sinh(freq*z)*sin(pfreq)) 

                enddo
             enddo
          enddo

          do k=lo(3),hi(3)
             z = prob_lo(3) + dx(3) * (dble(k) + half)
             do j=lo(2),hi(2)+1
                y = prob_lo(2) + dx(2) * (dble(j))
                do i=lo(1),hi(1)
                   x = prob_lo(1) + dx(1) * (dble(i) + half)

                   my(i,j,k) = 0.d0

                enddo
             enddo
          enddo

          do k=lo(3),hi(3)+1
             z = prob_lo(3) + dx(3) * (dble(k))
             do j=lo(2),hi(2)
                y = prob_lo(2) + dx(2) * (dble(j) + half)
                do i=lo(1),hi(1)
                   x = prob_lo(1) + dx(1) * (dble(i) + half)
                                
                   mz(i,j,k) = ufac*sin(freq*x)/cos(pfreq)*(cos(pfreq*z)*cosh(freq) &
                        -cosh(freq*z)*cos(pfreq))

                enddo
             enddo
          enddo

       else
          call bl_error('exact_solutions.f90: invalid prob_dir for prob_type=11')
       end if

       s(:,:,:,1) = 1.d0

    case (12)

       ! lo density spherical bubble
       do k=lo(3),hi(3)
          z = prob_lo(3) + dx(3) * (dble(k) + half) - 0.5d0*(prob_lo(3)+prob_hi(3))
          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2) * (dble(j) + half) - 0.5d0*(prob_lo(2)+prob_hi(2))
             do i=lo(1),hi(1)
                x = prob_lo(1) + dx(1) * (dble(i) + half) - 0.5d0*(prob_lo(1)+prob_hi(1))

                r = sqrt (x**2 + y**2 + z**2)

                ! tanh smoothing
                s(i,j,k,2) = c_init(1) + 0.5d0*(c_init(2)-c_init(1))* &
                     (1.d0 + tanh((r-2.5d0*smoothing_width*dx(1))/(smoothing_width*dx(1))))

                s(i,j,k,1) = 1.0d0/(s(i,j,k,2)/rhobar(1)+(1.0d0-s(i,j,k,2))/rhobar(2))
                s(i,j,k,2) = s(i,j,k,1)*s(i,j,k,2)
             enddo
          enddo
       end do

    case default
       call bl_error('ERROR: bad or unimplemented choice for prob_type in exact_3d')
    end select

  end subroutine exact_3d

  ! Note: This routine *increments* s_update: Do not reset it to zero! 
  subroutine exact_s_force_3d(s_update,s,ng_u,ng_s,lo,hi,dx,time)
 
    integer        , intent(in   ) :: lo(:),hi(:),ng_u,ng_s
    real(kind=dp_t), intent(inout) :: s_update(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:), time
    
    ! local
    integer :: i,j,k

    real(kind=dp_t) :: x, y, z, hx, hy, hz
    real(kind=dp_t) :: pfac, freq, pfreq 

    select case ( abs(prob_type) )
    case (2) ! Taylor (traveling wave) vortices

       freq  = 2.d0*M_PI
       pfac = dexp(-2.0d0*freq*freq*diff_coef*time)

       hx = ABC_coefs(1)
       hy = ABC_coefs(2)
       hz = ABC_coefs(3)

       do k=lo(3),hi(3)
          z = dx(3) * (dble(k) + half)
          do j=lo(2),hi(2)
             y = dx(2) * (dble(j) + half)
             do i=lo(1),hi(1)
                x = dx(1) * (dble(i) + half)

                s_update(i,j,k,2) = s_update(i,j,k,2) + &
                          pfac*(hx*hy*freq*cos(freq*(x-time))*cos(freq*(z-time)) - &
                          hy*hz*freq*sin(freq*(x-time))*sin(freq*(y-time)))

             enddo
          enddo
       enddo

       
    case (9) ! Sawtooth (triangle)   
       call bl_error('exact_s_force_3d: prob_type=9 not supported yet')
    end select

  end subroutine exact_s_force_3d

  ! Note: This routine *increments* m_update: Do not reset it to zero! 
  subroutine exact_m_force_3d(m_updatex,m_updatey,m_updatez,ng_u,lo,hi,dx,time)
 
    integer        , intent(in   ) :: lo(:),hi(:),ng_u
    real(kind=dp_t), intent(inout) :: m_updatex(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(inout) :: m_updatey(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(inout) :: m_updatez(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(in   ) :: dx(:), time
    
    select case ( abs(prob_type) )
    !not yet implemented
    end select

  end subroutine exact_m_force_3d

  ! 
  ! ---------------------------------------------------------------------------
  ! Variable-coefficients: chi
  ! ---------------------------------------------------------------------------
  ! 

  subroutine compute_chi(mla,chi,prim,dx)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: chi(:)
    type(multifab) , intent(in   ) :: prim(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)

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
          
  end subroutine compute_chi

  subroutine compute_chi_2d(chi,ng_c,prim,ng_p,lo,hi,dx)
    
    ! compute chi in valid AND ghost regions
    ! the ghost cells in prim have already been filled properly

    integer        , intent(in   ) :: lo(:), hi(:), ng_c, ng_p
    real(kind=dp_t), intent(inout) ::  chi(lo(1)-ng_c:,lo(2)-ng_c:)
    real(kind=dp_t), intent(in   ) :: prim(lo(1)-ng_p:,lo(2)-ng_p:,:)
    real (kind = dp_t), intent(in ) :: dx(:)

    ! Local
    integer :: i,j
    real(dp_t) :: x,y,stats(3),conc
    
    select case ( prob_type )
    case (-21) ! Keep rho*chi=const so as to get linear concentration profile
    
       do j=lo(2)-ng_c,hi(2)+ng_c
          do i=lo(1)-ng_c,hi(1)+ng_c
             chi(i,j) = abs(diff_coef) / prim(i,j,1)
          end do
       end do
      
    case (-100) ! Test the case of spatially-dependent diffusivity
      
       ! Spatially-varying transport coefficient:
       stats=(/0_dp_t, huge(0_dp_t), -huge(0_dp_t)/)
       do j=lo(2)-ng_c,hi(2)+ng_c
          y = dx(2) * (dble(j) + half)
          do i=lo(1)-ng_c,hi(1)+ng_c
             x = dx(1) * (dble(i) + half)
             chi(i,j) = abs(diff_coef) &
                * ( 0.75d0 + 0.25d0*sin(2*M_PI*(x-prob_lo(1))/(prob_hi(1)-prob_lo(1))) ) &
                * ( 0.75d0 + 0.25d0*sin(2*M_PI*(y-prob_lo(2))/(prob_hi(2)-prob_lo(2))) )
             stats = (/ stats(1)+chi(i,j), min(stats(2),chi(i,j)), max(stats(3),chi(i,j)) /) 
          end do
       end do
       stats(1)=stats(1)/size(chi)
       !write(*,*) "Mean/min/min chi = ", real(stats) ! TESTing only

    case (-1,-6,-7,-8,-9,-12)

       ! simple test function
       ! chi = chi0*(1 + a*c + b*c^2)
       do j=lo(2)-ng_c,hi(2)+ng_c
          do i=lo(1)-ng_c,hi(1)+ng_c
             conc = max(min(prim(i,j,2), 1.0_dp_t), 0.0_dp_t)
             chi(i,j) = abs(diff_coef)*(1 + material_properties(1,1)*conc + &
                                            material_properties(2,1)*conc**2)
          end do
       end do
      
    case default
      ! No need to do anything as chi has valid values already
    end select

  end subroutine compute_chi_2d

  subroutine compute_chi_3d(chi,ng_c,prim,ng_p,lo,hi,dx)
    
    integer        , intent(in   ) :: lo(:), hi(:), ng_c, ng_p
    real(kind=dp_t), intent(inout) ::  chi(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:)
    real(kind=dp_t), intent(in   ) :: prim(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:,:)
    real (kind = dp_t), intent(in ) :: dx(:)

    ! Local
    integer :: i,j,k
    real(dp_t) :: x,y,z,stats(3),conc
    
    select case ( prob_type )
    case (-21) ! Keep rho*chi=const so as to get linear concentration profile
    
       ! simple test function
       do k=lo(3)-ng_c,hi(3)+ng_c
          do j=lo(2)-ng_c,hi(2)+ng_c
             do i=lo(1)-ng_c,hi(1)+ng_c
               chi(i,j,k) = abs(diff_coef) / prim(i,j,k,1)
             end do  
          end do
       end do
      
    case (-100) ! Test the case of spatially-dependent diffusivity
      
       ! Spatially-varying transport coefficient:
       stats=(/0_dp_t, huge(0_dp_t), -huge(0_dp_t)/)
       
       do k=lo(3)-ng_c,hi(3)+ng_c
          z = dx(3) * (dble(k) + half)
       do j=lo(2)-ng_c,hi(2)+ng_c
          y = dx(2) * (dble(j) + half)
       do i=lo(1)-ng_c,hi(1)+ng_c
          x = dx(1) * (dble(i) + half)
          chi(i,j,k) = abs(diff_coef) &
             * ( 0.75d0 + 0.25d0*sin(2*M_PI*(x-prob_lo(1))/(prob_hi(1)-prob_lo(1))) ) &
             * ( 0.75d0 + 0.25d0*sin(2*M_PI*(y-prob_lo(2))/(prob_hi(2)-prob_lo(2))) ) &
             * ( 0.75d0 + 0.25d0*sin(2*M_PI*(z-prob_lo(3))/(prob_hi(3)-prob_lo(3))) )
          stats = (/ stats(1)+chi(i,j,k), min(stats(2),chi(i,j,k)), max(stats(3),chi(i,j,k)) /) 
       end do
       end do
       end do
       stats(1)=stats(1)/size(chi)
       !write(*,*) "Mean/min/min chi = ", real(stats) ! TESTing only
      
    case (-1,-6,-7,-8,-9,-12)

       ! simple test function
       do k=lo(3)-ng_c,hi(3)+ng_c
          do j=lo(2)-ng_c,hi(2)+ng_c
             do i=lo(1)-ng_c,hi(1)+ng_c
                conc = max(min(prim(i,j,k,2), 1.0_dp_t), 0.0_dp_t)
                chi(i,j,k) = abs(diff_coef)*(1 + material_properties(1,1)*conc + &
                                                 material_properties(2,1)*conc**2)
             end do
          end do
       end do

    case default
      ! No need to do anything as chi has valid values already
    end select

  end subroutine compute_chi_3d

  ! 
  ! ---------------------------------------------------------------------------
  ! Variable-coefficients: eta
  ! ---------------------------------------------------------------------------
  ! 

  subroutine compute_eta(mla,eta,prim,dx)

    ! compute eta in valid AND ghost regions
    ! the ghost cells in prim have already been filled properly

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: eta(:)
    type(multifab) , intent(in   ) :: prim(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)

    integer :: nlevs,dm,i,n,ng_e,ng_p
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: ep(:,:,:,:)
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
          ep => dataptr(eta(n), i)
          pp => dataptr(prim(n), i)
          lo = lwb(get_box(eta(n), i))
          hi = upb(get_box(eta(n), i))
          select case (dm)
          case (2)
             call compute_eta_2d(ep(:,:,1,1),ng_e,pp(:,:,1,:),ng_p,lo,hi,dx(n,:))
          case (3)
             call compute_eta_3d(ep(:,:,:,1),ng_e,pp(:,:,:,:),ng_p,lo,hi,dx(n,:))
          end select
       end do
    end do
          
  end subroutine compute_eta

  subroutine compute_eta_2d(eta,ng_e,prim,ng_p,lo,hi,dx)
    
    integer        , intent(in   ) :: lo(:), hi(:), ng_e, ng_p
    real(kind=dp_t), intent(inout) ::  eta(lo(1)-ng_e:,lo(2)-ng_e:)
    real(kind=dp_t), intent(in   ) :: prim(lo(1)-ng_p:,lo(2)-ng_p:,:)
    real (kind = dp_t), intent(in ) :: dx(:)

    ! Local
    integer :: i,j
    real(dp_t) :: x,y,stats(3),conc
    
    select case ( prob_type )
    case (-100) ! Test the case of spatially-dependent diffusivity
      
       ! Spatially-varying transport coefficient:
       stats=(/0_dp_t, huge(0_dp_t), -huge(0_dp_t)/)
       do j=lo(2)-ng_e,hi(2)+ng_e
          y = dx(2) * (dble(j) + half)
          do i=lo(1)-ng_e,hi(1)+ng_e
             x = dx(1) * (dble(i) + half)
             eta(i,j) = abs(visc_coef) &
                * ( 0.75d0 + 0.25d0*sin(2*M_PI*(x-prob_lo(1))/(prob_hi(1)-prob_lo(1))) ) &
                * ( 0.75d0 + 0.25d0*sin(2*M_PI*(y-prob_lo(2))/(prob_hi(2)-prob_lo(2))) )
             stats = (/ stats(1)+eta(i,j), min(stats(2),eta(i,j)), max(stats(3),eta(i,j)) /) 
          end do
       end do
       stats(1)=stats(1)/size(eta)
       !write(*,*) "Mean/min/min eta = ", real(stats) ! TESTing only

    case (-1,-6,-7,-8,-9,-12,-21)

       ! simple test function
       ! eta = eta0*(1+a*c+b*c^2)
       ! or
       ! eta = eta0*(1+a*c) / (1+b*c)
       stats=(/0_dp_t, huge(0_dp_t), -huge(0_dp_t)/)
       !stats=0
       do j=lo(2)-ng_e,hi(2)+ng_e
          do i=lo(1)-ng_e,hi(1)+ng_e
          conc = max(min(prim(i,j,2), 1.0_dp_t), 0.0_dp_t)
          if((prob_type==-1).or.(prob_type==-21)) then ! Quadratic function in c
             eta(i,j) = abs(visc_coef)*(1 + material_properties(1,2)*conc + &
                                            material_properties(2,2)*conc**2)
          else ! Rational function in c
             eta(i,j) = abs(visc_coef)*(1 + material_properties(1,2)*conc) / &
                                       (1 + material_properties(2,2)*conc)
             !stats = (/ stats(1)+eta(i,j), min(stats(2),eta(i,j)), max(stats(3),eta(i,j)) /) 
             stats = (/ stats(1)+eta(i,j)/prim(i,j,1), min(stats(2),eta(i,j)/prim(i,j,1)), max(stats(3),eta(i,j)/prim(i,j,1)) /)
             !stats = (/ stats(1)+eta(i,j), stats(2)+prim(i,j,1), stats(3)+eta(i,j)/prim(i,j,1) /) 
          end if
          end do
       end do
       stats(1)=stats(1)/size(eta)
       !stats=stats/size(eta)
       !write(*,*) "Mean/min/min nu = ", real(stats) ! TESTing only
      
    case default
      ! No need to do anything as eta has valid values already
    end select

  end subroutine compute_eta_2d

  subroutine compute_eta_3d(eta,ng_e,prim,ng_p,lo,hi,dx)
    
    integer        , intent(in   ) :: lo(:), hi(:), ng_e, ng_p
    real(kind=dp_t), intent(inout) ::  eta(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(in   ) :: prim(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:,:)
    real (kind = dp_t), intent(in ) :: dx(:)

    ! Local
    integer :: i,j,k
    real(dp_t) :: x,y,z,stats(3),conc
    
    select case ( prob_type )
    case (-100) ! Test the case of spatially-dependent diffusivity
      
       ! Spatially-varying transport coefficient:
       stats=(/0_dp_t, huge(0_dp_t), -huge(0_dp_t)/)
       
       do k=lo(3)-ng_e,hi(3)+ng_e
          z = dx(3) * (dble(k) + half)
       do j=lo(2)-ng_e,hi(2)+ng_e
          y = dx(2) * (dble(j) + half)
       do i=lo(1)-ng_e,hi(1)+ng_e
          x = dx(1) * (dble(i) + half)
          eta(i,j,k) = abs(visc_coef) &
             * ( 0.75d0 + 0.25d0*sin(2*M_PI*(x-prob_lo(1))/(prob_hi(1)-prob_lo(1))) ) &
             * ( 0.75d0 + 0.25d0*sin(2*M_PI*(y-prob_lo(2))/(prob_hi(2)-prob_lo(2))) ) &
             * ( 0.75d0 + 0.25d0*sin(2*M_PI*(z-prob_lo(3))/(prob_hi(3)-prob_lo(3))) )
          stats = (/ stats(1)+eta(i,j,k), min(stats(2),eta(i,j,k)), max(stats(3),eta(i,j,k)) /) 
       end do
       end do
       end do
       stats(1)=stats(1)/size(eta)
       !write(*,*) "Mean/min/min eta = ", real(stats) ! TESTing only
      
    case (-1,-6,-7,-8,-9,-12,-21)

       ! simple test function
       do k=lo(3)-ng_e,hi(3)+ng_e
       do j=lo(2)-ng_e,hi(2)+ng_e
       do i=lo(1)-ng_e,hi(1)+ng_e
         conc = max(min(prim(i,j,k,2), 1.0_dp_t), 0.0_dp_t)
         if((prob_type==-1).or.(prob_type==-21)) then ! Quadratic function in c
            eta(i,j,k) = abs(visc_coef)*(1 + material_properties(1,2)*conc + &
                                             material_properties(2,2)*conc**2)
         else ! Rational function in c
            eta(i,j,k) = abs(visc_coef)*(1 + material_properties(1,2)*conc) / &
                                        (1 + material_properties(2,2)*conc)
         end if
       end do
       end do
       end do
      
    case default
      ! No need to do anything as eta has valid values already
    end select

  end subroutine compute_eta_3d

  subroutine compute_kappa(mla,kappa,prim,dx)

    ! compute kappa in valid AND ghost regions
    ! the ghost cells in prim have already been filled properly

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: kappa(:)
    type(multifab) , intent(in   ) :: prim(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)

    integer :: nlevs,dm,i,n,ng_k,ng_p
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: kp(:,:,:,:)
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
          kp => dataptr(kappa(n), i)
          pp => dataptr(prim(n), i)
          lo = lwb(get_box(kappa(n), i))
          hi = upb(get_box(kappa(n), i))
          select case (dm)
          case (2)
             call compute_kappa_2d(kp(:,:,1,1),ng_k,pp(:,:,1,:),ng_p,lo,hi,dx(n,:))
          case (3)
             call compute_kappa_3d(kp(:,:,:,1),ng_k,pp(:,:,:,:),ng_p,lo,hi,dx(n,:))
          end select
       end do
    end do
          
  end subroutine compute_kappa

  subroutine compute_kappa_2d(kappa,ng_k,prim,ng_p,lo,hi,dx)
    
    integer        , intent(in   ) :: lo(:), hi(:), ng_k, ng_p
    real(kind=dp_t), intent(inout) :: kappa(lo(1)-ng_k:,lo(2)-ng_k:)
    real(kind=dp_t), intent(in   ) :: prim(lo(1)-ng_p:,lo(2)-ng_p:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

  end subroutine compute_kappa_2d

  subroutine compute_kappa_3d(kappa,ng_k,prim,ng_p,lo,hi,dx)
    
    integer        , intent(in   ) :: lo(:), hi(:), ng_k, ng_p
    real(kind=dp_t), intent(inout) :: kappa(lo(1)-ng_k:,lo(2)-ng_k:,lo(3)-ng_k:)
    real(kind=dp_t), intent(in   ) ::  prim(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

  end subroutine compute_kappa_3d

end module exact_solutions_module
