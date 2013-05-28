module inhomogeneous_bc_val_module

  use bl_types
  use probin_common_module,  only: wallspeed_hi, prob_hi
  use probin_lowmach_module, only: prob_type

  implicit none

  private

  public :: inhomogeneous_bc_val_2d, inhomogeneous_bc_val_3d

contains
 
 function inhomogeneous_bc_val_2d(comp,x,y) result(val)

   integer        , intent(in   ) :: comp
   real(kind=dp_t), intent(in   ) :: x,y
   real(kind=dp_t)                :: val

   select case (prob_type)
   case (1)
      if (y .eq. prob_hi(2) .and. comp .eq. 1) then
         val = wallspeed_hi(1,2)
      else
         val = 0.d0
      end if
   case default
      val = 0.d0
   end select

 end function inhomogeneous_bc_val_2d

 function inhomogeneous_bc_val_3d(comp,x,y,z) result(val)

   integer        , intent(in   ) :: comp
   real(kind=dp_t), intent(in   ) :: x,y,z
   real(kind=dp_t)                :: val

   select case (prob_type)
   case (1)
      if (z .eq. prob_hi(3)) then
         if (comp .eq. 1) then
            val = wallspeed_hi(1,3)
         else if (comp .eq. 2) then
            val = wallspeed_hi(2,3)
         else
            val = 0.d0
         end if
      else
         val = 0.d0
      end if
   case default
      val = 0.d0
   end select

 end function inhomogeneous_bc_val_3d

end module inhomogeneous_bc_val_module
