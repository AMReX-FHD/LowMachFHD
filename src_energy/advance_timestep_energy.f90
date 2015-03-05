module advance_timestep_energy_module

  implicit none

  private

  public :: advance_timestep_energy

contains

  ! this routine performs "Step 1" and "Step 2" of the algorithm:
  ! (Corrector Step)
  ! -Update P_0 and compute v^{*,n+1} using a GMRES solve
  ! -Advance rho_i and (rho h)
  ! -Compute volume discrepancy correction and return to beginning of this step
  !
  ! (Predictor Step)
  ! -Update P_0 and compute v^n with a GMRES solve
  ! -Advance rho_i and (rho h).
  ! -Compute volume discrepancy correction and return to beginning of this step
  subroutine advance_timestep_energy()

  end subroutine advance_timestep_energy

end module advance_timestep_energy_module
