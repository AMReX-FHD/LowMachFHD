module energy_eos_wrapper_module

  use ml_layout_module
  use multifab_physbc_module
  use bc_module
  use define_bc_module
  use energy_eos_module
  use convert_variables_module

  implicit none

  private

  public :: convert_c_to_x

contains

  subroutine convert_c_to_x()

  end subroutine convert_c_to_x

end module energy_eos_wrapper_module
