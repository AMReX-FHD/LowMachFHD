module bl_rng_module

  use bl_types
  use bl_random_module
  use parallel
  use bl_error_module
  use probin_common_module, only: restart
  use probin_reactdiff_module, only: temporal_integrator, diffusion_type, &
                                     use_Poisson_rng, seed_diffusion, seed_reaction, &
                                     seed_init, integer_populations
  use hg_rng_engine_module, only : hg_rng_engine

  implicit none

  private

  public :: rng_init, rng_destroy, &
            rng_eng_diffusion, &
            rng_eng_diffusion_e, &
            rng_eng_diffusion_d, &
            rng_eng_reaction, &
            rng_eng_reaction_e, &
            rng_eng_reaction_d, &
            rng_eng_init, &
            rng_dist_poisson_init, &
            rng_dist_normal_init

  ! randon number engines
  type(hg_rng_engine)      , save :: rng_eng_diffusion
  type(bl_rng_engine)      , save :: rng_eng_diffusion_e
  type(bl_rng_uniform_real), save :: rng_eng_diffusion_d
  type(hg_rng_engine)      , save :: rng_eng_reaction
  type(bl_rng_engine)      , save :: rng_eng_reaction_e
  type(bl_rng_uniform_real), save :: rng_eng_reaction_d

  !
  type(bl_rng_engine)      , save :: rng_eng_init
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
       call bl_rng_build_engine(rng_eng_diffusion_e, seed_diffusion)
       ! uniform real distribution: [0.d0, 1.d0)
       call bl_rng_build_distro(rng_eng_diffusion_d, 0.d0, 1.d0)
       rng_eng_diffusion%eng = rng_eng_diffusion_e%p
       rng_eng_diffusion%dis = rng_eng_diffusion_d%p
    end if
    
    if (seed_reaction .ne. -1) then
       call bl_rng_build_engine(rng_eng_reaction_e, seed_reaction)
       ! uniform real distribution: [0.d0, 1.d0)
       call bl_rng_build_distro(rng_eng_reaction_d, 0.d0, 1.d0)
       rng_eng_reaction%eng = rng_eng_reaction_e%p
       rng_eng_reaction%dis = rng_eng_reaction_d%p
    end if

    call bl_rng_build_engine(rng_eng_init, seed_init)
    ! poisson; mean=1; this will be overridden
    call bl_rng_build_distro(rng_dist_poisson_init, 1.d0)
    ! normal; mean=0, std=1
    call bl_rng_build_distro(rng_dist_normal_init, 0.d0, 1.d0)

  end subroutine rng_init

  subroutine rng_destroy()

    call bl_rng_destroy_engine(rng_eng_diffusion_e)
    call bl_rng_destroy_distro(rng_eng_diffusion_d)

    call bl_rng_destroy_engine(rng_eng_reaction_e)
    call bl_rng_destroy_distro(rng_eng_reaction_d)

    call bl_rng_destroy_engine(rng_eng_init)
    call bl_rng_destroy_distro(rng_dist_poisson_init)
    call bl_rng_destroy_distro(rng_dist_normal_init)

  end subroutine rng_destroy

end module bl_rng_module
