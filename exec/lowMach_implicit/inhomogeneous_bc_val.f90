module inhomogeneous_bc_val_module

  use bl_types
  use probin_common_module,  only: wallspeed_hi, prob_lo, prob_hi
  use probin_lowmach_module, only: prob_type, rhobar, c_bc

  implicit none

  private

  public :: inhomogeneous_bc_val_2d, inhomogeneous_bc_val_3d

contains
 
 function inhomogeneous_bc_val_2d(comp,x,y) result(val)

   integer        , intent(in   ) :: comp
   real(kind=dp_t), intent(in   ) :: x,y
   real(kind=dp_t)                :: val

   select case (prob_type)
   case (0)

      if (comp .eq. 4) then
         ! density
         if (x .eq. prob_lo(1)) then
            val = 1.d0/(c_bc(1,1)/rhobar(1) + (1.d0-c_bc(1,1))/rhobar(2))
         else if (x .eq. prob_hi(1)) then
            val = 1.d0/(c_bc(1,2)/rhobar(1) + (1.d0-c_bc(1,2))/rhobar(2))
         else if (y .eq. prob_lo(2)) then
            val = 1.d0/(c_bc(2,1)/rhobar(1) + (1.d0-c_bc(2,1))/rhobar(2))
         else if (y .eq. prob_hi(2)) then
            val = 1.d0/(c_bc(2,2)/rhobar(1) + (1.d0-c_bc(2,2))/rhobar(2))
         end if
      else if (comp .eq. 5) then
         ! concentration
         if (x .eq. prob_lo(1)) then
            val = c_bc(1,1)
         else if (x .eq. prob_hi(1)) then
            val = c_bc(1,2)
         else if (y .eq. prob_lo(2)) then
            val = c_bc(2,1)
         else if (y .eq. prob_hi(2)) then
            val = c_bc(2,2)
         end if
      else
         val = 0.d0
      end if

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
   case (0)

      if (comp .eq. 5) then
         ! density
         if (x .eq. prob_lo(1)) then
            val = 1.d0/(c_bc(1,1)/rhobar(1) + (1.d0-c_bc(1,1))/rhobar(2))
         else if (x .eq. prob_hi(1)) then
            val = 1.d0/(c_bc(1,2)/rhobar(1) + (1.d0-c_bc(1,2))/rhobar(2))
         else if (y .eq. prob_lo(2)) then
            val = 1.d0/(c_bc(2,1)/rhobar(1) + (1.d0-c_bc(2,1))/rhobar(2))
         else if (y .eq. prob_hi(2)) then
            val = 1.d0/(c_bc(2,2)/rhobar(1) + (1.d0-c_bc(2,2))/rhobar(2))
         else if (z .eq. prob_lo(3)) then
            val = 1.d0/(c_bc(3,1)/rhobar(1) + (1.d0-c_bc(3,1))/rhobar(2))
         else if (z .eq. prob_hi(3)) then
            val = 1.d0/(c_bc(3,2)/rhobar(1) + (1.d0-c_bc(3,2))/rhobar(2))
         end if
      else if (comp .eq. 6) then
         ! concentration
         if (x .eq. prob_lo(1)) then
            val = c_bc(1,1)
         else if (x .eq. prob_hi(1)) then
            val = c_bc(1,2)
         else if (y .eq. prob_lo(2)) then
            val = c_bc(2,1)
         else if (y .eq. prob_hi(2)) then
            val = c_bc(2,2)
         else if (z .eq. prob_lo(3)) then
            val = c_bc(3,1)
         else if (z .eq. prob_hi(3)) then
            val = c_bc(3,2)
         end if
      else
         val = 0.d0
      end if

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
