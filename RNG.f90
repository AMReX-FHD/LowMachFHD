module BoxLibRNGs
   ! These are written in C and part of BoxLib
   use Random_Numbers ! Used to generate unpredictable seeds or initial seeds for Marsenne twister
   use NonUniformRNGs ! Generates samples from several nonuniform distributions (we use our own for Guassian)
   implicit none
   public
   
   integer, parameter, private :: dp=kind(0.0d0), sp=kind(0.0)
   
   ! These are written in C and part of BoxLib
   interface SeedRNG_C
      subroutine srandgen(seed) bind(c)
         ! Seed the generator (should be called at the beginning)
         ! Use the shell script $RANDOM to generate a truly random seed
         integer, value :: seed
      end subroutine
   end interface

! This is not used at present and since it requires a C++ compiler I leave it out here   
!   interface SeedParallelRNG_C
!      subroutine SeedParallelRNG(seed) bind(c, name="SeedParallelRNG")
!         ! Seed the generator on each processor (should be called at the beginning)
!         ! This requires MPI to be initialized and running and is compiled with C++
!         integer, value :: seed
!      end subroutine
!   end interface   

   ! Donev: We have observed some difficulties calling this from Fortran
   ! It either has to do with unsigned versus signed integers or passing by value
   ! Appears to be a compiler bug in gfortran as it works in ifort.
   ! For now, best avoid this routine and use the floating-point RNG to compute integers
   interface UniformInteger
      ! void genrandint (unsigned long int *r, unsigned long int n)
      subroutine genrandint(number, range) bind(c)
         ! Returns an integer uniformly distributed in the range [1,range]
         import
         integer, intent(out) :: number
         integer, intent(in), value :: range
         !integer, intent(in) :: range ! Donev: Temporary change to pass by address here
      end subroutine   
   end interface

   interface UniformRNG
      subroutine genrand(number) bind(c)
         ! Returns pseudorandom number in interval [0,1).
         import
         real(dp), intent(out) :: number
      end subroutine   
      
      module procedure genrand_sp
   end interface
      
   interface NormalRNG  ! The Fortran version of this is below and may be faster
      subroutine genrandn(number) bind(c)
         ! Returns a normally-distributed number with mean 0 and variance 1
         import
         real(dp), intent(out) :: number
      end subroutine 
      
      module procedure genrandn_sp  
   end interface

   interface UniformRNGVec  ! The Fortran version of this is below and may be faster
      module procedure UniformRNGs_sp
      module procedure UniformRNGs
   end interface

   interface NormalRNGVec  ! The Fortran version of this is below and may be faster
      module procedure NormalRNGs_sp
      module procedure NormalRNGs
   end interface
   
   interface PoissonRNG ! This is only scalar for now
      module procedure PoissonRNG_sp
      module procedure PoissonRNG_dp
   end interface   

   interface BinomialRNG ! This is only scalar for now
      module procedure BinomialRNG_sp
      module procedure BinomialRNG_dp
   end interface   

contains ! It is likely that vectorized versions will do better here

   ! This is in principle callable by C directly, but we go through the wrapper here for safety
   ! void SeedRNG(int *seed); // On output, the actual seed used
   !
   subroutine SeedRNG(seed) BIND(C,NAME="SeedRNG")
      integer, intent(inout) :: seed ! If zero, the clock will be used to find an unpredictable seed

      if(seed==0) then
         call UnpredictableSeeds(seed)
         write(*,*) "SeedRNG @ BoxLibRNGs: Generated unpredictable SEED=", seed
      end if
      
      call SeedRNG_C(seed)

   end subroutine

   subroutine genrand_sp(number)
      ! Returns pseudorandom number in interval [0,1).
      real(sp), intent(out) :: number
      
      real(dp) :: number_dp
      
      call genrand(number_dp)
      number=number_dp
      ! In single precision, we may get 1.0 here so we need to do some hack      
      if(number>=1.0) number=number-epsilon(number)
      
   end subroutine   

   subroutine genrandn_sp(number)
      ! Returns pseudorandom number in interval [0,1).
      real(sp), intent(out) :: number
      
      real(dp) :: number_dp
      call genrandn(number_dp)
      number=number_dp
   end subroutine   

  subroutine UniformRNGs(numbers, n_numbers)
    integer, intent(in) :: n_numbers
    real(dp), intent(out) :: numbers(n_numbers)

    integer :: i

    do i=1, n_numbers
       call UniformRNG(numbers(i)) ! Marsenne-Twister in C
    end do   

  end subroutine

  subroutine UniformRNGs_sp(numbers, n_numbers)
    integer, intent(in) :: n_numbers
    real(sp), intent(out) :: numbers(n_numbers)

    integer :: i

    do i=1, n_numbers
       call UniformRNG(numbers(i)) ! Marsenne-Twister in C
    end do   

  end subroutine

  subroutine NormalRNGs(numbers, n_numbers)
    integer, intent(in) :: n_numbers
    real(dp), intent(out) :: numbers(n_numbers)

    integer :: i

    do i=1, n_numbers
       call NormalRNG(numbers(i))
    end do   

  end subroutine

  subroutine NormalRNGs_sp(numbers, n_numbers)
    integer, intent(in) :: n_numbers
    real(sp), intent(out) :: numbers(n_numbers)

    integer :: i

    do i=1, n_numbers
       call NormalRNG(numbers(i))
    end do   

  end subroutine 

  subroutine NormalRNG_Fortran(invnormdist)
      ! This is the Fortran equivalent of the C blinvnormdist, just for the record
      real(dp), intent(inout) :: invnormdist

      real(dp)     :: p

      real(dp) ::  q,r
      real(dp), parameter :: a1=-39.6968302866538d0
      real(dp), parameter :: a2=220.946098424521d0
      real(dp), parameter :: a3=-275.928510446969d0
      real(dp), parameter :: a4=138.357751867269d0
      real(dp), parameter :: a5=-30.6647980661472d0
      real(dp), parameter :: a6=2.50662827745924d0
      real(dp), parameter :: b1=-54.4760987982241d0
      real(dp), parameter :: b2=161.585836858041d0
      real(dp), parameter :: b3=-155.698979859887d0
      real(dp), parameter :: b4=66.8013118877197d0
      real(dp), parameter :: b5=-13.2806815528857d0
      real(dp), parameter :: c1=-0.00778489400243029d0
      real(dp), parameter :: c2=-0.322396458041136d0
      real(dp), parameter :: c3=-2.40075827716184d0
      real(dp), parameter :: c4=-2.54973253934373d0
      real(dp), parameter :: c5=4.37466414146497d0
      real(dp), parameter :: c6=2.93816398269878d0
      real(dp), parameter :: d1=0.00778469570904146d0
      real(dp), parameter :: d2=0.32246712907004d0
      real(dp), parameter :: d3=2.445134137143d0
      real(dp), parameter :: d4=3.75440866190742d0
      real(dp), parameter :: p_low=0.02425d0
      real(dp), parameter :: p_high=0.9575d0
      
      call UniformRNG(p)
      
      if(p.lt.p_low) then
         q=dsqrt(-2.d0*dlog(p))
         invnormdist = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/   &
             ((((d1*q+d2)*q+d3)*q+d4)*q+1.d0)
      elseif (p.le.p_high)then
         q=p-0.5d0
         r=q*q
         invnormdist = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/  &
                   (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.d0)
      else
         q=dsqrt(-2.d0*dlog(1.d0-p))
         invnormdist = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/   &
               ((((d1*q+d2)*q+d3)*q+d4)*q+1.d0)
      endif

  end subroutine

  subroutine NormalRNGFast(p)
  real(dp)     :: u,p
  real(dp), parameter :: f = 3.46410161514d0
  
  call UniformRNG(u)
  p = f*(u-0.5_dp)

  end subroutine 

 SUBROUTINE PoissonRNG_dp(number,mean)
    INTEGER, INTENT(OUT) :: number
    REAL(dp), INTENT(IN) :: mean

    number=random_Poisson(mu=real(mean), first=.true.)

 END SUBROUTINE

 SUBROUTINE PoissonRNG_sp(number,mean)
    INTEGER, INTENT(OUT) :: number
    REAL(sp), INTENT(IN) :: mean

    number=random_Poisson(mu=mean, first=.true.)

 END SUBROUTINE

 SUBROUTINE BinomialRNG_dp(number,n_trials,success_prob)
    INTEGER, INTENT(OUT) :: number
    INTEGER, INTENT(IN) :: n_trials ! Number of trials
    REAL(dp), INTENT(IN) :: success_prob ! Probability of successful trial

    number=random_binomial(n=n_trials, pp=real(success_prob), first=.true.)

 END SUBROUTINE

 SUBROUTINE BinomialRNG_sp(number,n_trials,success_prob)
    INTEGER, INTENT(OUT) :: number
    INTEGER, INTENT(IN) :: n_trials ! Number of trials
    REAL(sp), INTENT(IN) :: success_prob ! Probability of successful trial

    number=random_binomial(n=n_trials, pp=success_prob, first=.true.)
    
 END SUBROUTINE
  
end module
