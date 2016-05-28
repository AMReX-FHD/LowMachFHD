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
            rng_uniform_real_reaction

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
       ! multinomial diffusion

    else
       ! fluctuating hydro

    end if

    ! reactions
    if (reaction_type .eq. 0 .or. reaction_type .eq. 1) then
       if (use_Poisson_rng .eq. 1) then
          ! tau-leaping (initialize mean to 1; this will be overridden)
          call bl_rng_build(rng_poisson_reaction,seed_reaction,1.d0)
       else if (use_Poisson_rng .eq. 0) then
          ! CLE

       end if
    else if (reaction_type .eq. 1) then
       ! SSA

    end if


  end subroutine rng_init

  subroutine rng_destroy()

    ! diffusion
    if (temporal_integrator .eq. -3 .or. diffusion_type .eq. 3) then
       ! multinomial diffusion

    else
       ! fluctuating hydro

    end if

    ! reactions
    if (reaction_type .eq. 0 .or. reaction_type .eq. 1) then
       if (use_Poisson_rng .eq. 1) then
          ! tau-leaping
          call bl_rng_destroy(rng_poisson_reaction)
       else if (use_Poisson_rng .eq. 0) then
          ! CLE

       end if
    else if (reaction_type .eq. 1) then
       ! SSA

    end if

  end subroutine rng_destroy

end module bl_rng_module
