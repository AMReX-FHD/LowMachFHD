module bl_rng_module

  use bl_types
  use bl_random_module
  use probin_reactdiff_module, only: temporal_integrator, diffusion_type, reaction_type, &
                                     use_Poisson_rng, seed_diffusion, seed_reaction, &
                                     seed_init

  implicit none

  private

  public :: rng_init, rng_destroy, &
            rng_binomial_diffusion, &
            rng_normal_diffusion, &
            rng_poisson_reaction, &
            rng_normal_reaction, &
            rng_uniform_real_reaction, &
            bl_MultinomialRNG

  ! diffusion
  type(bl_rng_binomial)    , save :: rng_binomial_diffusion     ! multinomial
  type(bl_rng_normal)      , save :: rng_normal_diffusion       ! fluctuating hydro

  ! reactions
  type(bl_rng_poisson)     , save :: rng_poisson_reaction       ! tau-leaping
  type(bl_rng_normal)      , save :: rng_normal_reaction        ! CLE
  type(bl_rng_uniform_real), save :: rng_uniform_real_reaction  ! SSA

  ! initialization
  type(bl_rng_normal)      , save :: rng_normal_init

contains

  subroutine rng_init()

    ! diffusion
    if (temporal_integrator .eq. -3 .or. diffusion_type .eq. 3) then
       ! multinomial diffusion - calls a sequence of binomial random numbers
       ! (initialize the trials to 1 and probability to 0.5; this will be overridden)
       call bl_rng_build(rng_binomial_diffusion,seed_diffusion,1,0.5d0)
    else
       ! fluctuating hydro (mean 0, standard deviation 1)
       call bl_rng_build(rng_normal_diffusion,seed_diffusion,0.d0,1.d0)
    end if

    ! reactions
    if (reaction_type .eq. 0 .or. reaction_type .eq. 1) then
       if (use_Poisson_rng .eq. 1) then
          ! tau-leaping (initialize mean to 1; this will be overridden)
          call bl_rng_build(rng_poisson_reaction,seed_reaction,1.d0)
       else if (use_Poisson_rng .eq. 0) then
          ! CLE (mean 0, standard deviation 1)
          call bl_rng_build(rng_normal_reaction,seed_reaction,0.d0,1.d0)
       end if
    else if (reaction_type .eq. 2) then
       ! SSA (in interval [0,1))
       call bl_rng_build(rng_uniform_real_reaction,seed_reaction,0.d0,1.d0)
    end if


  end subroutine rng_init

  subroutine rng_destroy()

    ! diffusion
    if (temporal_integrator .eq. -3 .or. diffusion_type .eq. 3) then
       ! multinomial diffusion
       call bl_rng_destroy(rng_binomial_diffusion)
    else
       ! fluctuating hydro
       call bl_rng_destroy(rng_normal_diffusion)
    end if

    ! reactions
    if (reaction_type .eq. 0 .or. reaction_type .eq. 1) then
       if (use_Poisson_rng .eq. 1) then
          ! tau-leaping
          call bl_rng_destroy(rng_poisson_reaction)
       else if (use_Poisson_rng .eq. 0) then
          ! CLE
          call bl_rng_destroy(rng_normal_reaction)
       end if
    else if (reaction_type .eq. 1) then
       ! SSA
       call bl_rng_destroy(rng_uniform_real_reaction)
    end if

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
