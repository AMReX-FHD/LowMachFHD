module BoxLibRNGs
   use RNGEngine ! Interfaces for uniform/integer RNG generators written in C and part of BoxLib
   use Random_Numbers ! Used to generate unpredictable seeds or initial seeds for Marsenne twister
   use NonUniformRNGs ! Generates samples from several nonuniform distributions (we use our own for Guassian)
   use iso_c_binding
   implicit none
   public ! Export all names of public entities so this is the one unified Fortran RNG interface for the outside world
   
   integer, parameter, private :: dp=kind(0.0d0), sp=kind(0.0)

   ! Returns pseudorandom number in interval [0,1).
   ! Donev: The single precision generator should be tested at least once!
   interface UniformRNG
      module procedure UniformRNG_dp
      module procedure UniformRNG_sp
   end interface UniformRNG
      
   interface NormalRNG  ! The Fortran version of this is below and may be faster
      module procedure NormalRNG_dp
      module procedure NormalRNG_sp  
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

   subroutine UniformRNG_dp(number, engine)
      ! Returns pseudorandom number in interval [0,1).
      real(dp), intent(out) :: number
      type(hg_rng_engine), intent(inout), optional :: engine
      if (present(engine)) then
         call UniformRNG_C(number, engine)
      else
         call UniformRNG_C(number)
      end if
    end subroutine 

   subroutine UniformRNG_sp(number, engine)
      ! Returns pseudorandom number in interval [0,1).
      real(sp), intent(out) :: number
      type(hg_rng_engine), intent(inout), optional :: engine
      
      real(dp) :: number_dp
      
      if (present(engine)) then
         ! Donev: Consider adding support for single-precision reals here
         call UniformRNG_C(number_dp, engine)
         number=number_dp
         ! In single precision, we may get 1.0 here so we need to do some hack      
         if(number>=1.0) number=number-epsilon(number)
      else
         call UniformRNG_C(number)
      end if
      
   end subroutine   

   ! Donev: Consider replacing this with NormalRNG_Fortran for this and removing NormalRNG_C from the code  
   subroutine NormalRNG_dp(number, engine)
      ! Returns pseudorandom number in interval [0,1).
      real(dp), intent(out) :: number
      type(hg_rng_engine), intent(inout), optional :: engine
      if (present(engine)) then
         call NormalRNG_C(number, engine)
      else
         call NormalRNG_C(number)
      end if
   end subroutine   

   subroutine NormalRNG_sp(number, engine)
      ! Returns pseudorandom number in interval [0,1).
      real(sp), intent(out) :: number
      type(hg_rng_engine), intent(inout), optional :: engine
      real(dp) :: number_dp
      if (present(engine)) then
         call NormalRNG_C(number_dp, engine)
      else
         call NormalRNG_C(number_dp)
      end if
      number=number_dp
   end subroutine   

  subroutine UniformRNGs(numbers, n_numbers, engine)
    integer, intent(in) :: n_numbers
    real(dp), intent(out) :: numbers(n_numbers)
    type(hg_rng_engine), intent(inout), optional :: engine

    integer :: i

    do i=1, n_numbers
       call UniformRNG(numbers(i), engine) ! Marsenne-Twister in C
    end do   

  end subroutine

  subroutine UniformRNGs_sp(numbers, n_numbers, engine)
    integer, intent(in) :: n_numbers
    real(sp), intent(out) :: numbers(n_numbers)
    type(hg_rng_engine), intent(inout), optional :: engine

    integer :: i

    do i=1, n_numbers
       call UniformRNG(numbers(i), engine) ! Marsenne-Twister in C
    end do   

  end subroutine

  subroutine NormalRNGs(numbers, n_numbers, engine)
    integer, intent(in) :: n_numbers
    real(dp), intent(out) :: numbers(n_numbers)
    type(hg_rng_engine), intent(inout), optional :: engine

    integer :: i

    do i=1, n_numbers
       call NormalRNG(numbers(i), engine)
    end do   

  end subroutine

  subroutine NormalRNGs_sp(numbers, n_numbers, engine)
    integer, intent(in) :: n_numbers
    real(sp), intent(out) :: numbers(n_numbers)
    type(hg_rng_engine), intent(inout), optional :: engine

    integer :: i

    do i=1, n_numbers
       call NormalRNG(numbers(i), engine)
    end do   

  end subroutine 

  ! Donev: Consider renaming this NormalRNG_dp so it is used instead of the C code
  ! Donev: It should be tested for histogram at least once or checked that given the same random integer
  !   it produces the same result as the C code
  subroutine NormalRNG_Fortran(invnormdist, engine)
      ! This is the Fortran equivalent of the C blinvnormdist, just for the record
      real(dp), intent(inout) :: invnormdist
      type(hg_rng_engine), intent(inout), optional :: engine

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
      
      call UniformRNG(p, engine)
      
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

  subroutine NormalRNGFast(p, engine) ! This has the correct first and second moment only
    ! It is not an actual normal random generator!
    type(hg_rng_engine), intent(inout), optional :: engine
    real(dp)     :: u,p
    real(dp), parameter :: f = 3.46410161514d0

    call UniformRNG(u, engine)
    p = f*(u-0.5_dp)

  end subroutine 

 SUBROUTINE PoissonRNG_dp(number,mean,engine)
    INTEGER, INTENT(OUT) :: number
    REAL(dp), INTENT(IN) :: mean
    type(hg_rng_engine), intent(inout), optional :: engine

    number=random_Poisson(mu=real(mean), engine=engine)

 END SUBROUTINE

 SUBROUTINE PoissonRNG_sp(number,mean,engine)
    INTEGER, INTENT(OUT) :: number
    REAL(sp), INTENT(IN) :: mean
    type(hg_rng_engine), intent(inout), optional :: engine

    number=random_Poisson(mu=mean, engine=engine)

 END SUBROUTINE

 SUBROUTINE BinomialRNG_dp(number,n_trials,success_prob, engine)
    INTEGER, INTENT(OUT) :: number
    INTEGER, INTENT(IN) :: n_trials ! Number of trials
    REAL(dp), INTENT(IN) :: success_prob ! Probability of successful trial
    type(hg_rng_engine), intent(inout), optional :: engine

    number=random_binomial(n=n_trials, pp=real(success_prob), engine=engine)

 END SUBROUTINE

 SUBROUTINE BinomialRNG_sp(number,n_trials,success_prob, engine)
    INTEGER, INTENT(OUT) :: number
    INTEGER, INTENT(IN) :: n_trials ! Number of trials
    REAL(sp), INTENT(IN) :: success_prob ! Probability of successful trial
    type(hg_rng_engine), intent(inout), optional :: engine

    number=random_binomial(n=n_trials, pp=success_prob, engine=engine)
    
 END SUBROUTINE

 ! This samples from a multinomial distribution
 ! The last sample is not sampled explicitly since it is just N-sum(samples)
 subroutine MultinomialRNG(samples, n_samples, N, p, engine)
    integer, intent(in) :: n_samples, N
    integer, intent(out) :: samples(n_samples)
    real(dp), intent(in) :: p(n_samples)
    type(hg_rng_engine), intent(inout), optional :: engine

    real(dp) :: sum_p
    integer :: sample, sum_n

    if(sum(p)>1.0_dp) stop "Sum of probabilities must be less than 1"

    sum_p=0
    sum_n=0
    do sample=1, n_samples
       call BinomialRNG(number=samples(sample), n_trials=N-sum_n, &
               success_prob=p(sample)/(1.0_dp-sum_p), engine=engine)
       sum_n = sum_n + samples(sample)
       sum_p = sum_p + p(sample)
    end do      

 end subroutine
  
end module
