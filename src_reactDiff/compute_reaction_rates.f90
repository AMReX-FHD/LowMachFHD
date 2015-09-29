module compute_reaction_rates_module

  use bl_types

  implicit none

  private

  public :: compute_reaction_rates

contains

  subroutine compute_reaction_rates(n_in,reaction_rates,dv)
    
    real(kind=dp_t), intent(in   ) :: n_in(:),dv
    real(kind=dp_t), intent(inout) :: reaction_rates(:)

    ! reaction 1: n1 + n2 -> n3
    reaction_rates(1) = 0.0001d0*n_in(1)*n_in(2)

    ! reaction 2: n3 -> n1 + n2
    reaction_rates(2) = 0.d0
    
  end subroutine compute_reaction_rates

end module compute_reaction_rates_module
