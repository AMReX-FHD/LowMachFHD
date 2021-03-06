module bl_rng_module

  use bl_types
  use bl_random_module
  use parallel
  use bl_error_module
  use probin_common_module, only: restart, seed_momentum, seed_diffusion, &
                                  seed_reaction, seed_init_mass, seed_init_momentum

  implicit none

  private

  public :: rng_init, rng_destroy, &
            rng_eng_momentum, &
            rng_eng_diffusion, &
            rng_eng_diffusion_chk, &
            rng_eng_reaction, &
            rng_eng_reaction_chk, &
            rng_eng_init_mass, &
            rng_eng_init_momentum

  ! randon number engines
  type(bl_rng_engine)      , save :: rng_eng_momentum
  type(bl_rng_engine)      , save :: rng_eng_diffusion
  type(bl_rng_engine)      , save :: rng_eng_diffusion_chk
  type(bl_rng_engine)      , save :: rng_eng_reaction
  type(bl_rng_engine)      , save :: rng_eng_reaction_chk
  type(bl_rng_engine)      , save :: rng_eng_init_mass
  type(bl_rng_engine)      , save :: rng_eng_init_momentum

contains

  subroutine rng_init()

    if (seed_momentum .eq. -1 .and. restart .lt. 0) then
       call bl_error("seed_momentum = -1 requires restart")
    end if

    if (seed_diffusion .eq. -1 .and. restart .lt. 0) then
       call bl_error("seed_diffusion = -1 requires restart")
    end if

    if (seed_reaction .eq. -1 .and. restart .lt. 0) then
       call bl_error("seed_reaction = -1 requires restart")
    end if

    !!!!!!!!!!!!!!!!!!!!!!
    ! build engines
    !!!!!!!!!!!!!!!!!!!!!!

    ! momentum
    if (seed_momentum .ne. -1) then
       call bl_rng_build_engine(rng_eng_momentum, seed_momentum)
    end if

    ! mass diffusion
    if (seed_diffusion .ne. -1) then
       call bl_rng_build_engine(rng_eng_diffusion, seed_diffusion)
    end if
    ! build this - seed doesn't matter since this engine is overwritten
    call bl_rng_build_engine(rng_eng_diffusion_chk, 1)
    
    ! reactions
    if (seed_reaction .ne. -1) then
       call bl_rng_build_engine(rng_eng_reaction, seed_reaction)
    end if
    ! build this - seed doesn't matter since this engine is overwritten
    call bl_rng_build_engine(rng_eng_reaction_chk, 1)

    ! initialization
    call bl_rng_build_engine(rng_eng_init_mass,     seed_init_mass)
    call bl_rng_build_engine(rng_eng_init_momentum, seed_init_momentum)

  end subroutine rng_init

  subroutine rng_destroy()

    call bl_rng_destroy_engine(rng_eng_momentum)
    call bl_rng_destroy_engine(rng_eng_diffusion)
    call bl_rng_destroy_engine(rng_eng_diffusion_chk)
    call bl_rng_destroy_engine(rng_eng_reaction)
    call bl_rng_destroy_engine(rng_eng_reaction_chk)
    call bl_rng_destroy_engine(rng_eng_init_mass)
    call bl_rng_destroy_engine(rng_eng_init_momentum)

  end subroutine rng_destroy

end module bl_rng_module
