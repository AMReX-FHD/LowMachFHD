module inhomogeneous_bc_val_module

  use bl_types

  implicit none

  private

  public :: inhomogeneous_bc_val_2d, inhomogeneous_bc_val_3d

contains
 
 function inhomogeneous_bc_val_2d(comp,x,y) result(val)

   integer        , intent(in   ) :: comp
   real(kind=dp_t), intent(in   ) :: x,y
   real(kind=dp_t)                :: val

   call bl_error("inhomogeneous_bc_val: write a custom routine and put in local source directory"

 end function inhomogeneous_bc_val_2d

 function inhomogeneous_bc_val_3d(comp,x,y,z) result(val)

   integer        , intent(in   ) :: comp
   real(kind=dp_t), intent(in   ) :: x,y,z
   real(kind=dp_t)                :: val

   call bl_error("inhomogeneous_bc_val: write a custom routine and put in local source directory"

 end function inhomogeneous_bc_val_3d

end module inhomogeneous_bc_val_module
