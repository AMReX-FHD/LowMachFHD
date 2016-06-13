module bl_rng_module

  use bl_types
  use bl_random_module
  use parallel
  use bl_error_module
  use probin_common_module, only: restart
  use probin_reactdiff_module, only: temporal_integrator, diffusion_type, reaction_type, &
                                     use_Poisson_rng, seed_diffusion, seed_reaction, &
                                     seed_init, integer_populations

  implicit none

  private

  public :: rng_init, rng_destroy, &
            rng_binomial_diffusion, &
            rng_normal_diffusion, &
            rng_poisson_reaction, &
            rng_normal_reaction, &
            rng_uniform_real_reaction, &
            rng_poisson_init, &
            rng_normal_init, &
            bl_MultinomialRNG

  ! diffusion
  type(bl_rng_binomial)    , save :: rng_binomial_diffusion     ! multinomial
  type(bl_rng_normal)      , save :: rng_normal_diffusion       ! fluctuating hydro

  ! reactions
  type(bl_rng_poisson)     , save :: rng_poisson_reaction       ! tau-leaping
  type(bl_rng_normal)      , save :: rng_normal_reaction        ! CLE
  type(bl_rng_uniform_real), save :: rng_uniform_real_reaction  ! SSA

  ! initialization (limited support for now)
  type(bl_rng_poisson)      , save :: rng_poisson_init
  type(bl_rng_normal)       , save :: rng_normal_init

contains

  subroutine rng_init()

    if (seed_diffusion .eq. -1 .and. restart .lt. 0) then
       call bl_error("seed_diffusion = -1 requires restart")
    end if

    if (seed_reaction .eq. -1 .and. restart .lt. 0) then
       call bl_error("seed_reaction = -1 requires restart")
    end if

    !!!!!!!!!!!!!!!!!!
    ! diffusion
    !!!!!!!!!!!!!!!!!!

    if (seed_diffusion .eq. 0) then
       seed_diffusion = bl_rng_random_uint_c()
       call parallel_bcast(seed_diffusion)
       if (parallel_IOProcessor()) then
          print*,'seed_diffusion = 0 --> picking a random root seed'
       end if
    end if
    if (parallel_IOProcessor()) then
       print*,'seed_diffusion =',seed_diffusion
    end if

    if (seed_diffusion .ne. -1) then

       ! multinomial diffusion - calls a sequence of binomial random numbers
       ! (initialize the trials to 1 and probability to 0.5; this will be overridden)
       call bl_rng_build(rng_binomial_diffusion,seed_diffusion,1,0.5d0)

       ! fluctuating hydro (mean 0, standard deviation 1)
       call bl_rng_build(rng_normal_diffusion,seed_diffusion,0.d0,1.d0)

    end if

    !!!!!!!!!!!!!!!!!!
    ! reactions
    !!!!!!!!!!!!!!!!!!

    if (seed_reaction .eq. 0) then
       seed_reaction = bl_rng_random_uint_c()
       call parallel_bcast(seed_reaction)
       if (parallel_IOProcessor()) then
          print*,'seed_reaction = 0 --> picking a random root seed'
       end if
    end if
    if (parallel_IOProcessor()) then
       print*,'seed_reaction =',seed_reaction
    end if

    if (seed_diffusion .ne. -1) then

       ! tau-leaping (initialize mean to 1; this will be overridden)
       call bl_rng_build(rng_poisson_reaction,seed_reaction,1.d0)

       ! CLE (mean 0, standard deviation 1)
       call bl_rng_build(rng_normal_reaction,seed_reaction,0.d0,1.d0)

       ! SSA (in interval [0,1))
       call bl_rng_build(rng_uniform_real_reaction,seed_reaction,0.d0,1.d0)

    end if

    !!!!!!!!!!!!!!!!!!
    ! initilization
    !!!!!!!!!!!!!!!!!!

    if (seed_init .eq. 0) then
       seed_init = bl_rng_random_uint_c()
       call parallel_bcast(seed_init)
       if (parallel_IOProcessor()) then
          print*,'seed_init = 0 --> picking a random root seed'
       end if
    end if
    if (parallel_IOProcessor()) then
       print*,'seed_init =',seed_init
    end if
    
    ! random integer population (initialize mean to 1; this will be overridden)
    call bl_rng_build(rng_poisson_init,seed_init,1.d0)

    ! Gaussian noise (mean 0, standard deviation 1)
    call bl_rng_build(rng_normal_init,seed_init,0.d0,1.d0)

  end subroutine rng_init

  subroutine rng_destroy()

    !!!!!!!!!!!!!!!!!!
    ! diffusion
    !!!!!!!!!!!!!!!!!!

    ! multinomial diffusion
    call bl_rng_destroy(rng_binomial_diffusion)

    ! fluctuating hydro
    call bl_rng_destroy(rng_normal_diffusion)

    !!!!!!!!!!!!!!!!!!
    ! reactions
    !!!!!!!!!!!!!!!!!!

    ! tau-leaping
    call bl_rng_destroy(rng_poisson_reaction)

    ! CLE
    call bl_rng_destroy(rng_normal_reaction)

    ! SSA
    call bl_rng_destroy(rng_uniform_real_reaction)

    !!!!!!!!!!!!!!!!!!
    ! initialization
    !!!!!!!!!!!!!!!!!!

    ! random integer population
    call bl_rng_destroy(rng_poisson_init)

    ! Gaussian noise
    call bl_rng_destroy(rng_normal_init)

  end subroutine rng_destroy

  ! This samples from a multinomial distribution
  ! The last sample is not sampled explicitly since it is just N-sum(samples)
  subroutine bl_MultinomialRNG(samples, n_samples, N, p)

    integer   , intent(in)  :: n_samples, N
    integer   , intent(out) :: samples(n_samples)
    real(dp_t), intent(in)  :: p(n_samples)

    real(dp_t) :: sum_p
    integer :: sample, sum_n

    if(sum(p)>1.d0) stop "Sum of probabilities must be less than 1"

    sum_p=0
    sum_n=0
    do sample=1, n_samples

       call bl_rng_change_distribution(rng_binomial_diffusion, &
                                       N-sum_n,p(sample)/(1.d0-sum_p))
       samples(sample) = bl_rng_get(rng_binomial_diffusion)
       sum_n = sum_n + samples(sample)
       sum_p = sum_p + p(sample)
    end do      

 end subroutine

end module bl_rng_module
