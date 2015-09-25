module SSA_Gene
use iso_c_binding
use gsl_interface
use NonUniformRNGs
use ContinuousPoisson
implicit none
public

...

! Use negative value of method for second-order (-1 or -2)
integer :: method=3 ! 0=SSA, 1=tau leaping, 2=Chemical Langevin Equation (CLE)
logical :: use_gsl=.false.
integer, parameter :: wp=kind(0.0d0)
integer, parameter :: nreactions=4, nspecies=2 ! In the real code we want these to vary at runtime!

integer*8 :: n_bad_samples=0, n_rejections=0

real(wp), parameter :: theta=0.5_wp ! Parameter for Mattingly predictor-correoctor method
real(wp) :: alpha1, alpha2

contains

subroutine SSA_Gene_update( x, k, dV, dt, n_steps )
   real(wp), intent(inout)  :: x(nspecies) ! State
   real(wp), intent(in   )  :: k(nreactions), dV
   real(wp), intent(inout)  :: dt
      ! On input: Time interval to simulate, negative if only one step desired
      ! On output: Last sampled time in SSA
   integer, intent(out), optional :: n_steps ! How many SSA steps were taken
   
   real(wp) :: t, kRates(nreactions), r(nreactions), rTotal, rSum, tau, alpha1, alpha2
   real(wp) :: mean, dice, n_react(nreactions), rates_p(nreactions)
   real(wp) :: n_X(nspecies)  ! Number of molecules, if desired
   integer  :: i, j, iReaction
   real :: rr
   
   ! Table of stochiometric coefficients: to be read from input namelist
   real(wp), dimension(nspecies,nreactions), parameter :: stochiometric = &
      reshape((/ 1, 0, &
                 0, 1, &
                -1, 0, &
                 0, -1  /), (/nspecies,nreactions/));

   t = 0; ! Local time
   n_X = x*dV ! Convert to number of molecules instead of number densities (not required)
               
   if(present(n_steps)) n_steps=0

if(method==0) then ! SSA
   
   EventLoop: do
       call compute_rates()
       
       !write(*,*) "propensity = ", real(r(1)), real(r(3:8)); stop
       rTotal = sum(r);
          
       ! Generate exp. dist. random number with mean = 1/rSum
       call random_uniform(rr)
       tau = -log(1-rr)/rTotal;
       t = t + tau;

       ! Exit routine if t > dt (Normal exit)
       if( (dt>0.0_wp) .and. (t > dt) ) exit EventLoop

       ! Select the next reaction according to relative rates
       call random_uniform(rr)
       rr = rr*rTotal;
       rSum = 0;
       FindReaction: do i=1,size(r)
           rSum = rSum + r(i);
           iReaction = i;
           if( rSum >= rr ) then
               exit FindReaction
           end if
       end do FindReaction

       ! Change particle numbers accor
       n_X = n_X + stochiometric(:,iReaction)
       x = n_X / dV
              
       if(present(n_steps)) n_steps = n_steps + 1
       !write(*,*) t, iReaction, n_U, n_V, n_W

       if( dt<0.0_wp ) then ! Done
           dt = t;
           exit EventLoop
       end if

    end do EventLoop
     
else ! Some form of tau leaping: execute all events in the time interval at once

    call compute_rates()
    rates_p = r ! Save these for predictor stage
    
    ! Decide how many reactions will happen for each reaction by sampling a Poisson number
    do iReaction=1, nreactions
    if(r(iReaction)>0.0_wp) then
   
    ! Predictor step (Euler-Maruyama so first-order if only predictor):
        
       mean=r(iReaction)*dt ! Mean number of events
       !write(*,*) "predictor=", mean
       if(mean<0) then
          write(*,*) "Negative mean in tau leaping rates=", r, " state=", n_X
          mean=0.0_wp
       end if
       
       if(method<0) mean=tau*mean ! Only go up to fraction tau of time step
    
       select case (abs(method))
       case(1) ! Traditional tau leaping 

          if(use_gsl) then
             n_react(iReaction)=gsl_ran_poisson (base_rng, real(mean,c_double)) ! Use GSL
          else   
             n_react(iReaction)=random_Poisson(real(mean), first=.true.) ! Faster but sometimes (rarely) segfaults?
          end if   
              
       case(2) ! CLE
       
          dice=random_normal()
          n_react(iReaction)=mean+sqrt(mean)*dice
          
       case default
          stop "Unsupported method selection"   
       end select   
       !write(*,*) iReaction, n_react(iReaction), kRates(iReaction)

       ! Change particle numbers accordingly
       n_X = n_X + n_react(iReaction)*stochiometric(:,iReaction);
       x = n_X / dV
    end if   
    end do
    !write(*,*) "end=", n_x, n_react
    n_steps = sum(n_react)

if(method<0) then ! Do corrector stage of second-order Anderson/Mattingly method

   alpha1=1.0_wp/(2*theta*(1-theta));
   alpha2=alpha1-1.0_wp;
      
   call compute_rates()
      
   ! Decide how many reactions will happen for each reaction by sampling a Poisson number
   do iReaction=1, nreactions
   if(r(iReaction)>0.0_wp) then
   
      mean = (alpha1*rates_p(iReaction)-alpha2*r(iReaction))*(1.0_wp-theta)*dt      
      !write(*,*) "corrector=", mean
      if(mean<0) then
         n_rejections = n_rejections + 1
         mean=0         
      end if

      select case (method)
      case(-1) ! Traditional tau leaping 

          if(use_gsl) then
             n_react(iReaction)=gsl_ran_poisson (base_rng, real(mean,c_double)) ! Use GSL
          else   
             n_react(iReaction)=random_Poisson(real(mean), first=.true.) ! Faster but sometimes (rarely) segfaults?
          end if   
         
      case(-2) ! CLE  

         dice=random_normal()
         n_react(iReaction)=mean+sqrt(mean)*dice

      case default
         stop "Unsupported method selection in corrector stage"   
      end select   
      !write(*,*) iReaction, n_react(iReaction), kRates(iReaction)

      ! Change particle numbers accordingly
      n_X = n_X + n_react(iReaction)*stochiometric(:,iReaction);
      x = n_X / dV
   end if
   end do
   !write(*,*) "end=", n_x, n_react
   n_steps = n_steps + sum(n_react)
     
end if

end if

contains

   subroutine compute_rates() ! Computes reaction rates given number densities or numbers of molecules
      ! This needs to compute stochiometric factors -- I need to write code to do this in a general setting
   end subroutine
      
end subroutine

end module
