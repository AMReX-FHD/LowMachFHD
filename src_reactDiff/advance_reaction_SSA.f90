module advance_reaction_SSA_module

  use ml_layout_module
  use bc_module
  use define_bc_module
  use multifab_physbc_module
  use BoxLibRNGs
  use bl_rng_module
  use bl_random_module
  use compute_reaction_rates_module
  use probin_reactdiff_module, only: nspecies, nreactions, stoichiometric_factors, &
                                     use_bl_rng

  implicit none

  private

  public :: advance_reaction_SSA_cell
  
contains

  ! advance_reaction_SSA solves dn/dt = f(n) 
  !  where f(n) are the chemical production rates (deterministic or stochastic)

  subroutine advance_reaction_SSA_cell(n_old,n_new,dv,dt)

    real(kind=dp_t), intent(in   ) :: n_old(:)
    real(kind=dp_t), intent(inout) :: n_new(:)
    real(kind=dp_t), intent(in   ) :: dv, dt

    real(kind=dp_t) :: avg_react_rate(1:nreactions)
    real(kind=dp_t) :: rTotal, rr, rSum, tau, t_local

    integer :: spec,reaction, which_reaction
    integer :: n_steps_SSA
    
    ! copy old state into new
    n_new = n_old

    t_local = 0.d0
    n_steps_SSA = 0

    EventLoop: do

       ! compute reaction rates in units (# reactions) / (unit time) / (unit volume)
       call compute_reaction_rates(n_new(1:nspecies),avg_react_rate,dv)

       ! compute reaction rates in units (# reactions) / (unit time)
       avg_react_rate = max(0.0d0,avg_react_rate*dv)

       ! sum the reaction rates
       rTotal = sum(avg_react_rate(1:nreactions))

       ! generate pseudorandom number in interval [0,1).
       if (use_bl_rng) then
          call UniformRNG(rr, rng_eng_reaction)
       else
          call UniformRNG(rr)
       end if
       ! tau is how long until the next reaction occurs
       tau = -log(1-rr)/rTotal
       t_local = t_local + tau;

       if (t_local .gt. dt) exit EventLoop

       ! Select the next reaction according to relative rates
       if (use_bl_rng) then
          call UniformRNG(rr, rng_eng_reaction)
       else
          call UniformRNG(rr)
       end if
       rr = rr*rTotal
       rSum = 0
       FindReaction: do reaction=1,nreactions
          rSum = rSum + avg_react_rate(reaction)
          which_reaction = reaction
          if( rSum >= rr ) exit FindReaction
       end do FindReaction

       ! update number densities for this reaction
       do spec=1,nspecies
          n_new(spec) = n_new(spec) + &
               (stoichiometric_factors(spec,2,which_reaction)-stoichiometric_factors(spec,1,which_reaction)) / dv
       end do
       
       n_steps_SSA = n_steps_SSA+1

    end do EventLoop

  end subroutine advance_reaction_SSA_cell

end module advance_reaction_SSA_module
