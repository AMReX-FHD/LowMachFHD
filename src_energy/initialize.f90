module initialize_module

  implicit none

  private

  public :: initialize

contains

  ! this routine performs "Step 0" of the algorithm:
  ! -Update P_0 and compute v^n with a projection
  ! -Advance rho_i and (rho h).
  ! -Compute volume discrepancy correction and return to beginning of this step
  subroutine initialize()

  end subroutine initialize

end module initialize_module
