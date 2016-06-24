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
            rng_eng_diffusion, &
            rng_eng_reaction, &
            rng_eng_init, &
            rng_dist_binomial_diffusion, &
            rng_dist_normal_diffusion, &
            rng_dist_poisson_reaction, &
            rng_dist_normal_reaction, &
            rng_dist_uniform_real_reaction, &
            rng_dist_poisson_init, &
            rng_dist_normal_init, &
            bl_MultinomialRNG

  ! randon number engines
  type(bl_rng_engine), save :: rng_eng_diffusion
  type(bl_rng_engine), save :: rng_eng_reaction
  type(bl_rng_engine), save :: rng_eng_init

  ! distributions
  type(bl_rng_binomial)    , save :: rng_dist_binomial_diffusion
  type(bl_rng_normal)      , save :: rng_dist_normal_diffusion
  type(bl_rng_poisson)     , save :: rng_dist_poisson_reaction
  type(bl_rng_normal)      , save :: rng_dist_normal_reaction
  type(bl_rng_uniform_real), save :: rng_dist_uniform_real_reaction
  type(bl_rng_poisson)     , save :: rng_dist_poisson_init
  type(bl_rng_normal)      , save :: rng_dist_normal_init

contains

  subroutine rng_init()

    if (seed_diffusion .eq. -1 .and. restart .lt. 0) then
       call bl_error("seed_diffusion = -1 requires restart")
    end if

    if (seed_reaction .eq. -1 .and. restart .lt. 0) then
       call bl_error("seed_reaction = -1 requires restart")
    end if

    !!!!!!!!!!!!!!!!!!!!!!
    ! build engines
    !!!!!!!!!!!!!!!!!!!!!!

    if (seed_diffusion .ne. -1) then
       call bl_rng_build_engine(rng_eng_diffusion, seed_diffusion)
    end if
    
    if (seed_reaction .ne. -1) then
       call bl_rng_build_engine(rng_eng_reaction, seed_reaction)
    end if

    call bl_rng_build_engine(rng_eng_init, seed_init)

    !!!!!!!!!!!!!!!!!!!!!!
    ! build distributions
    !!!!!!!!!!!!!!!!!!!!!!

    ! binomial; t=1 (trials), p=0.5; this will be overridden
    call bl_rng_build_distro(rng_dist_binomial_diffusion, 1, 0.5d0)

    ! normal; mean=0, std=1
    call bl_rng_build_distro(rng_dist_normal_diffusion, 0.d0, 1.d0)

    ! poisson; mean=1; this will be overridden
    call bl_rng_build_distro(rng_dist_poisson_reaction, 1.d0)

    ! normal; mean=0, std=1
    call bl_rng_build_distro(rng_dist_normal_reaction, 0.d0, 1.d0)

    ! uniform real distribution: [0.d0, 1.d0)
    call bl_rng_build_distro(rng_dist_uniform_real_reaction, 0.d0, 1.d0)

    ! poisson; mean=1; this will be overridden
    call bl_rng_build_distro(rng_dist_poisson_init, 1.d0)

    ! normal; mean=0, std=1
    call bl_rng_build_distro(rng_dist_normal_init, 0.d0, 1.d0)

  end subroutine rng_init

  subroutine rng_destroy()

    call bl_rng_destroy_engine(rng_eng_diffusion)
    call bl_rng_destroy_engine(rng_eng_reaction)
    call bl_rng_destroy_engine(rng_eng_init)

    call bl_rng_destroy_distro(rng_dist_binomial_diffusion)
    call bl_rng_destroy_distro(rng_dist_normal_diffusion)
    call bl_rng_destroy_distro(rng_dist_poisson_reaction)
    call bl_rng_destroy_distro(rng_dist_normal_reaction)
    call bl_rng_destroy_distro(rng_dist_uniform_real_reaction)
    call bl_rng_destroy_distro(rng_dist_poisson_init)
    call bl_rng_destroy_distro(rng_dist_normal_init)

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

       call bl_rng_destroy_distro(rng_dist_binomial_diffusion)
       call bl_rng_build_distro(rng_dist_binomial_diffusion,N-sum_n,p(sample)/(1.d0-sum_p))
       samples(sample) = bl_rng_get(rng_dist_binomial_diffusion,rng_eng_diffusion)
       sum_n = sum_n + samples(sample)
       sum_p = sum_p + p(sample)
    end do      

 end subroutine

end module bl_rng_module
