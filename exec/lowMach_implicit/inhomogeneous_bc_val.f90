module inhomogeneous_bc_val_module

  use bl_types
  use bc_module
  use probin_common_module,  only: wallspeed_lo, wallspeed_hi, prob_lo, prob_hi
  use probin_lowmach_module, only: rhobar, c_bc

  implicit none

  private

  public :: inhomogeneous_bc_val_2d, inhomogeneous_bc_val_3d

contains
 
 function inhomogeneous_bc_val_2d(comp,x,y) result(val)

   integer        , intent(in   ) :: comp
   real(kind=dp_t), intent(in   ) :: x,y
   real(kind=dp_t)                :: val

   if (comp .eq. vel_bc_comp) then
      ! x-vel

      if (y .eq. prob_lo(2)) then
         val = wallspeed_lo(1,2)
      else if (y .eq. prob_hi(2)) then
         val = wallspeed_hi(1,2)
      else
         val = 0.d0
      end if
   
   else if (comp .eq. vel_bc_comp+1) then
      ! y-vel

      if (x .eq. prob_lo(1)) then
         val = wallspeed_lo(1,1)
      else if (x .eq. prob_hi(1)) then
         val = wallspeed_hi(1,1)
      else
         val = 0.d0
      end if

   else if (comp .eq. scal_bc_comp) then
      ! density

      if (x .eq. prob_lo(1)) then
         val = 1.d0/(c_bc(1,1)/rhobar(1) + (1.d0-c_bc(1,1))/rhobar(2))
      else if (x .eq. prob_hi(1)) then
         val = 1.d0/(c_bc(1,2)/rhobar(1) + (1.d0-c_bc(1,2))/rhobar(2))
      else if (y .eq. prob_lo(2)) then
         val = 1.d0/(c_bc(2,1)/rhobar(1) + (1.d0-c_bc(2,1))/rhobar(2))
      else if (y .eq. prob_hi(2)) then
         val = 1.d0/(c_bc(2,2)/rhobar(1) + (1.d0-c_bc(2,2))/rhobar(2))
      else
         val = 0.d0
      end if
      
   else if (comp .eq. scal_bc_comp+1) then
      ! concentration

      if (x .eq. prob_lo(1)) then
         val = c_bc(1,1)
      else if (x .eq. prob_hi(1)) then
         val = c_bc(1,2)
      else if (y .eq. prob_lo(2)) then
         val = c_bc(2,1)
      else if (y .eq. prob_hi(2)) then
         val = c_bc(2,2)
      else
         val = 0.d0
      end if

   else
      val = 0.d0
   end if

 end function inhomogeneous_bc_val_2d

 function inhomogeneous_bc_val_3d(comp,x,y,z) result(val)

   integer        , intent(in   ) :: comp
   real(kind=dp_t), intent(in   ) :: x,y,z
   real(kind=dp_t)                :: val

   if (comp .eq. vel_bc_comp) then
      ! x-vel

      if (y .eq. prob_lo(2)) then
         val = wallspeed_lo(1,2)
      else if (y .eq. prob_hi(2)) then
         val = wallspeed_hi(1,2)
      else if (z .eq. prob_lo(3)) then
         val = wallspeed_lo(1,3)
      else if (z .eq. prob_hi(3)) then
         val = wallspeed_hi(1,3)
      else
         val = 0.d0
      end if

   else if (comp .eq. vel_bc_comp+1) then
      ! y-vel

      if (x .eq. prob_lo(1)) then
         val = wallspeed_lo(1,1)
      else if (x .eq. prob_hi(1)) then
         val = wallspeed_hi(1,1)
      else if (z .eq. prob_lo(3)) then
         val = wallspeed_lo(2,3)
      else if (z .eq. prob_hi(3)) then
         val = wallspeed_hi(2,3)
      else
         val = 0.d0
      end if

   else if (comp .eq. vel_bc_comp+2) then
      ! z-vel

      if (x .eq. prob_lo(1)) then
         val = wallspeed_lo(2,1)
      else if (x .eq. prob_hi(1)) then
         val = wallspeed_hi(2,1)
      else if (y .eq. prob_lo(2)) then
         val = wallspeed_lo(2,2)
      else if (y .eq. prob_hi(2)) then
         val = wallspeed_hi(2,2)
      else
         val = 0.d0
      end if

   else if (comp .eq. scal_bc_comp) then
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
      else
         val = 0.d0
      end if

   else if (comp .eq. scal_bc_comp+1) then
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
      else
         val = 0.d0
      end if

   else
      val = 0.d0
   end if

 end function inhomogeneous_bc_val_3d

end module inhomogeneous_bc_val_module
