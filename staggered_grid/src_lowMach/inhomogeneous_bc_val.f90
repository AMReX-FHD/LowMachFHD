module inhomogeneous_bc_val_module

  use bl_types
  use bc_module
  use bl_error_module
  use bl_constants_module
  use probin_multispecies_module, only: c_bc, rhotot_bc
  use probin_common_module, only: prob_lo, prob_hi, wallspeed_lo, wallspeed_hi, prob_type, &
                                  nspecies, algorithm_type, rho0, n_cells
  use probin_charged_module, only: Epot_wall, Epot_wall_bc_type, zero_charge_on_wall_type, bc_function_type, &
                                   L_pos, L_trans, L_zero, induced_charge_eo, E_ext_value, ac_iceo

  implicit none

  private

  public :: scalar_bc, transport_bc, inhomogeneous_bc_val_2d, inhomogeneous_bc_val_3d, alternating_current_efield

contains

  subroutine scalar_bc(phys_bc, bc_code)

    ! variable ordering for low Mach code (num_scal_bc = 2*nspecies+3):
    ! 1                       = total density
    ! 2:nspecies+1            = concentrations
    ! nspecies+2:2*nspecies+1 = molfrac or massfrac (dimensionless fractions
    ! 2*nspecies+2            = temperature
    ! 2*nspecies+3            = electric potential

    integer, intent(in   ) :: phys_bc
    integer, intent(inout) :: bc_code(1:num_scal_bc)

    ! set bc types for everything except electric potential, which is handled in define_bc_tower.f90
    if ((phys_bc == NO_SLIP_WALL) .or. (phys_bc == SLIP_WALL)) then

       bc_code(1:num_scal_bc-2) = FOEXTRAP  ! Pure Neumann for total density, conctractions, mol/mass fractions
       bc_code(  num_scal_bc-1) = EXT_DIR   ! But temperature is still specified at the boundary (via call to multifab_coefbc)

    else if ((phys_bc == NO_SLIP_RESERVOIR) .or. (phys_bc == SLIP_RESERVOIR)) then

       bc_code(1:num_scal_bc-1) = EXT_DIR   ! Pure Dirichlet

    else if (phys_bc == PERIODIC .or. phys_bc == INTERIOR ) then

       ! retain the default value of INTERIOR

    else

       ! return an error
       bc_code = -999

    end if

  end subroutine scalar_bc

  subroutine transport_bc(phys_bc, bc_code)

    ! we use num_tran_bc=1; same bc for all transport coefficients

    integer, intent(in   ) :: phys_bc
    integer, intent(inout) :: bc_code(1:num_tran_bc)

    if ((phys_bc == NO_SLIP_WALL) .or. (phys_bc == SLIP_WALL)) then

       bc_code = FOEXTRAP  ! Pure Neumann

    else if ((phys_bc == NO_SLIP_RESERVOIR) .or. (phys_bc == SLIP_RESERVOIR)) then

       bc_code = EXT_DIR   ! Pure Dirichlet

    else if (phys_bc == PERIODIC .or. phys_bc == INTERIOR ) then

       ! retain the default value of INTERIOR

    else

       ! return an error
       bc_code = -999

    end if

  end subroutine transport_bc

  function inhomogeneous_bc_val_2d(comp,x,y,time_in) result(val)

    integer        , intent(in   ) :: comp
    real(kind=dp_t), intent(in   ) :: x,y
    real(kind=dp_t), intent(in), optional :: time_in
    real(kind=dp_t)                :: val

    real(kind=dp_t) :: time, Lx, delta_x, s_j   !s_j is translated value in transition regions: s_j = x - (L_pos) or 
                                                     ! s_j = x - (L_pos+L_trans+L_zero) depending on if we're in transition region
                                                     ! number 1 or 2. 

    if (present(time_in)) then
       time = time_in
    else
       time = 0.d0
    end if

    if (comp .eq. vel_bc_comp) then

       ! x-vel
       if (y .eq. prob_lo(2)) then
          if (abs(prob_type) .eq. 12) then
             if (time .le. 0.5d0) then
                val = wallspeed_lo(1,2)* &
                     0.5d0*(1.d0 + sin(2.d0*M_PI*x - 0.5d0*M_PI)) &
                     *0.5d0*(1.d0 + sin(2.d0*M_PI*time - 0.5d0*M_PI))
             else
                val = wallspeed_lo(1,2)* &
                     0.5d0*(1.d0 + sin(2.d0*M_PI*x - 0.5d0*M_PI))
             end if
          else
             val = wallspeed_lo(1,2)
          end if
       else if (y .eq. prob_hi(2)) then
          if (abs(prob_type) .eq. 12) then
             if (time .le. 0.5d0) then
                val = wallspeed_hi(1,2)* &
                     0.5d0*(1.d0 + sin(2.d0*M_PI*x - 0.5d0*M_PI)) &
                     *0.5d0*(1.d0 + sin(2.d0*M_PI*time - 0.5d0*M_PI))
             else
                val = wallspeed_hi(1,2)* &
                     0.5d0*(1.d0 + sin(2.d0*M_PI*x - 0.5d0*M_PI))
             end if
          else
             val = wallspeed_hi(1,2)
          end if
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

    else if (comp .eq. pres_bc_comp) then

       ! pressure
       val = 0.d0

    else if (comp .eq. scal_bc_comp) then

       ! full density (only for boussinesq algorithm)
       if (algorithm_type .eq. 6) then

          ! full rho boundary condition
          if (x .eq. prob_lo(1)) then
             val = rhotot_bc(1,1)
          else if (x .eq. prob_hi(1)) then
             val = rhotot_bc(1,2)
          else if (y .eq. prob_lo(2)) then
             val = rhotot_bc(2,1)
          else if (y .eq. prob_hi(2)) then
             val = rhotot_bc(2,2)
          else
             val = 0.d0
          end if

       else
          call bl_error("calling inhomogeneous_bc_val_2d with scal_bc_comp (full rho) - only works with algorithm_type=6")
       end if

    else if (comp .ge. c_bc_comp .and. comp .le. c_bc_comp+nspecies-1) then

       ! c_i boundary condition
       if (x .eq. prob_lo(1)) then
          val = c_bc(1,1,comp-c_bc_comp+1)
       else if (x .eq. prob_hi(1)) then
          val = c_bc(1,2,comp-c_bc_comp+1)
       else if (y .eq. prob_lo(2)) then
          val = c_bc(2,1,comp-c_bc_comp+1)
       else if (y .eq. prob_hi(2)) then
          val = c_bc(2,2,comp-c_bc_comp+1)
       else
          val = 0.d0
       end if

    else if (comp .eq. Epot_bc_comp) then
       ! electric potential 
       if (bc_function_type.eq.0) then     ! Constant potential
          if (x .eq. prob_lo(1)) then
             val = Epot_wall(1,1)
          else if (x .eq. prob_hi(1)) then
             val = Epot_wall(2,1)
          else if (y .eq. prob_lo(2)) then

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! if we're doing ICEO, we want to phi_tot = phi_periodic + phi_ext = 0 on the boundary 
             ! Since phi_ext = -x*E_0 + B, where E_0 is the magnitude of the external field and 
             ! B is some (irrelevant) constant, we can rig phi_tot to vanish on the portion of the 
             ! boundary we want (the metallic strip aka a conductor) by prescribing that
             ! phi_periodic = x*E_0. 
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             if (induced_charge_eo) then 
                if (ac_iceo) then 
                   val = 1.d0*(x-(prob_lo(1)+prob_hi(1))/2.d0)*alternating_current_efield(time)
                else
                   val = 1.d0*(x-(prob_lo(1)+prob_hi(1))/2.d0)*E_ext_value(1)  
                endif 
             else 
                val = Epot_wall(1,2) 
             end if 

             if (zero_charge_on_wall_type .eq. 1) then
                Lx = prob_hi(1)-prob_lo(1)
                if (x .le. Lx/4.d0 .or. x .ge. 3.d0*Lx/4.d0) then
                   val = 0.d0
                end if
             end if

          else if (y .eq. prob_hi(2)) then
             val = Epot_wall(2,2)
          else
             val = 0.d0
          end if
       else if (bc_function_type.eq.1) then                        ! Piecewise cubic potential--note we were lazy and assumed prob_lo(1) = 0.d0
          if (y .eq. prob_lo(2)) then                              ! lower boundary 
             if (x .lt. L_pos) then                                ! we're in first positive region
                val = Epot_wall(1,2)
             else if (x .lt. (L_pos + L_trans)) then               ! we're in first transition region
                delta_x = prob_hi(1)/n_cells(1)
                s_j = x - L_pos
                val = 1.d0/delta_x*(0.5d0*Epot_wall(1,2)/L_trans**3*((s_j + delta_x/2.d0)**4 - (s_j - delta_x/2.d0)**4) - & 
                               Epot_wall(1,2)/L_trans**2*((s_j + delta_x/2.d0)**3 - (s_j - delta_x/2.d0)**3)) + Epot_wall(1,2)
             else if (x .lt. (L_pos + L_trans + L_zero)) then      ! we're in zero region
                val = 0.d0 
             else if (x .lt. (L_pos + 2.d0*L_trans + L_zero)) then !we're in second transition region
                delta_x = prob_hi(1)/n_cells(1)
                s_j = x - (L_pos + L_trans + L_zero)
                val = 1.d0/delta_x*(-0.5d0*Epot_wall(1,2)/L_trans**3*((s_j + delta_x/2.d0)**4 - (s_j - delta_x/2.d0)**4) + & 
                               Epot_wall(1,2)/L_trans**2*((s_j + delta_x/2.d0)**3 - (s_j - delta_x/2.d0)**3))
             else                                                  ! we're in second positive region
                val = Epot_wall(1,2)
             end if 
          else if (y .eq. prob_hi(2)) then                         ! upper boundary
             if (x .lt. L_pos) then                                ! we're in first positive region
                val = Epot_wall(2,2)
             else if (x .lt. (L_pos + L_trans)) then               ! we're in first transition region
                delta_x = prob_hi(1)/n_cells(1)
                s_j = x - L_pos
                val = 1.d0/delta_x*(0.5d0*Epot_wall(2,2)/L_trans**3*((s_j + delta_x/2.d0)**4 - (s_j - delta_x/2.d0)**4) - & 
                               Epot_wall(2,2)/L_trans**2*((s_j + delta_x/2.d0)**3 - (s_j - delta_x/2.d0)**3)) + Epot_wall(2,2)
             else if (x .lt. (L_pos + L_trans + L_zero)) then      ! we're in zero region
                val = 0.d0 
             else if (x .lt. (L_pos + 2.d0*L_trans + L_zero)) then !we're in second transition region
                delta_x = prob_hi(1)/n_cells(1)
                s_j = x - (L_pos + L_trans + L_zero)
                val = 1.d0/delta_x*(-0.5d0*Epot_wall(2,2)/L_trans**3*((s_j + delta_x/2.d0)**4 - (s_j - delta_x/2.d0)**4) + & 
                               Epot_wall(2,2)/L_trans**2*((s_j + delta_x/2.d0)**3 - (s_j - delta_x/2.d0)**3))
             else                                                  ! we're in second positive region
                val = Epot_wall(2,2)
             end if 
          else 
             val = 0.d0
          end if 
       else 
          call bl_error("Only bc_function_type 0 and 1 currently supported.")
       end if
    else

       print*,'comp=',comp
       call bl_error("calling inhomogeneous_bc_val_2d with invalid comp")

    end if


  end function inhomogeneous_bc_val_2d

  ! This function allows us to specify an external electric field that oscillates in time. 
  ! The current field is a mollified square wave made by stitching together hyperbolic tangents. 
  !
  ! Currently the period and frequency of the wave is hard-coded in this function. 
  function alternating_current_efield(time) result(val)

    ! inputs/output
    real(kind=dp_t), intent(in)    :: time
    real(kind=dp_t)                :: val

    ! local variables 
    real(kind=dp_t) :: omega, T, shift ! omega controls how quickly the applied field goes from 
                                       ! positive to negative, and vice-versa. The larger the value,
                                       ! the less quickly it transitions.
                                       ! T is the period of oscillation. 
    integer         :: branch                

    ! initialize T and omega 
    T = 1.0d-4               ! these values were selected for a proof of concept AC-ICEO simulation. 
    omega = 0.05d0*T        


    branch = modulo(int(floor(time/T)), 2) !tells us if we use tanh or -1*tanh
    shift = floor(time/T)*T

    if (branch.eq.0) then 
       val = -1.d0*E_ext_value(1)*tanh(((time-shift) - T/2.d0)/omega)
    else 
       val =       E_ext_value(1)*tanh(((time-shift) - T/2.d0)/omega)
    endif  

  end function alternating_current_efield

  function inhomogeneous_bc_val_3d(comp,x,y,z,time_in) result(val)

    integer        , intent(in   ) :: comp
    real(kind=dp_t), intent(in   ) :: x,y,z
    real(kind=dp_t), intent(in), optional :: time_in
    real(kind=dp_t)                :: val

    real(kind=dp_t) :: time

    if (present(time_in)) then
       time = time_in
    else
       time = 0.d0
    end if

    if (comp .eq. vel_bc_comp) then

       ! x-vel
       if (y .eq. prob_lo(2)) then
          if (abs(prob_type) .eq. 12) then
             if (time .le. 0.5d0) then
                val = wallspeed_lo(1,2)* &
                     0.5d0*(1.d0 + sin(2.d0*M_PI*x - 0.5d0*M_PI)) &
                     *0.5d0*(1.d0 + sin(2.d0*M_PI*z - 0.5d0*M_PI)) &
                     *0.5d0*(1.d0 + sin(2.d0*M_PI*time - 0.5d0*M_PI))
             else
                val = wallspeed_lo(1,2)* &
                     0.5d0*(1.d0 + sin(2.d0*M_PI*x - 0.5d0*M_PI)) &
                     *0.5d0*(1.d0 + sin(2.d0*M_PI*z - 0.5d0*M_PI))
             end if
          else
             val = wallspeed_lo(1,2)
          end if
       else if (y .eq. prob_hi(2)) then
          if (abs(prob_type) .eq. 12) then
             if (time .le. 0.5d0) then
                val = wallspeed_hi(1,2)* &
                     0.5d0*(1.d0 + sin(2.d0*M_PI*x - 0.5d0*M_PI)) &
                     *0.5d0*(1.d0 + sin(2.d0*M_PI*z - 0.5d0*M_PI)) &
                     *0.5d0*(1.d0 + sin(2.d0*M_PI*time - 0.5d0*M_PI))
             else
                val = wallspeed_hi(1,2)* &
                     0.5d0*(1.d0 + sin(2.d0*M_PI*x - 0.5d0*M_PI)) &
                     *0.5d0*(1.d0 + sin(2.d0*M_PI*z - 0.5d0*M_PI))
             end if
          else
             val = wallspeed_hi(1,2)
          end if
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
          if (abs(prob_type) .eq. 12) then
             if (time .le. 0.5d0) then
                val = wallspeed_lo(2,2)* &
                     0.5d0*(1.d0 + sin(2.d0*M_PI*z - 0.5d0*M_PI)) &
                     *0.5d0*(1.d0 + sin(2.d0*M_PI*x - 0.5d0*M_PI)) &
                     *0.5d0*(1.d0 + sin(2.d0*M_PI*time - 0.5d0*M_PI))
             else
                val = wallspeed_lo(2,2)* &
                     0.5d0*(1.d0 + sin(2.d0*M_PI*z - 0.5d0*M_PI)) &
                     *0.5d0*(1.d0 + sin(2.d0*M_PI*x - 0.5d0*M_PI))
             end if
          else
             val = wallspeed_lo(2,2)
          end if
       else if (y .eq. prob_hi(2)) then
          if (abs(prob_type) .eq. 12) then
             if (time .le. 0.5d0) then
                val = wallspeed_hi(2,2)* &
                     0.5d0*(1.d0 + sin(2.d0*M_PI*z - 0.5d0*M_PI)) &
                     *0.5d0*(1.d0 + sin(2.d0*M_PI*x - 0.5d0*M_PI)) &
                     *0.5d0*(1.d0 + sin(2.d0*M_PI*time - 0.5d0*M_PI))
             else
                val = wallspeed_hi(2,2)* &
                     0.5d0*(1.d0 + sin(2.d0*M_PI*z - 0.5d0*M_PI)) &
                     *0.5d0*(1.d0 + sin(2.d0*M_PI*x - 0.5d0*M_PI))
             end if
          else
             val = wallspeed_hi(2,2)
          end if
       else
          val = 0.d0
       end if

    else if (comp .eq. pres_bc_comp) then

       ! pressure
       val = 0.d0

    else if (comp .eq. scal_bc_comp) then

       ! full density (only for boussinesq algorithm)
       if (algorithm_type .eq. 6) then

          ! full rho boundary condition
          if (x .eq. prob_lo(1)) then
             val = rhotot_bc(1,1)
          else if (x .eq. prob_hi(1)) then
             val = rhotot_bc(1,2)
          else if (y .eq. prob_lo(2)) then
             val = rhotot_bc(2,1)
          else if (y .eq. prob_hi(2)) then
             val = rhotot_bc(2,2)
          else if (z .eq. prob_lo(3)) then
             val = rhotot_bc(3,1)
          else if (z .eq. prob_hi(3)) then
             val = rhotot_bc(3,2)
          else
             val = 0.d0
          end if
       else
          call bl_error("calling inhomogeneous_bc_val_3d with scal_bc_comp (full rho) - only works with algorithm_type=6")
       end if

    else if (comp .ge. c_bc_comp .and. comp .le. c_bc_comp+nspecies-1) then

       ! c_i boundary condition
       if (x .eq. prob_lo(1)) then
          val = c_bc(1,1,comp-c_bc_comp+1)
       else if (x .eq. prob_hi(1)) then
          val = c_bc(1,2,comp-c_bc_comp+1)
       else if (y .eq. prob_lo(2)) then
          val = c_bc(2,1,comp-c_bc_comp+1)
       else if (y .eq. prob_hi(2)) then
          val = c_bc(2,2,comp-c_bc_comp+1)
       else if (z .eq. prob_lo(3)) then
          val = c_bc(3,1,comp-c_bc_comp+1)
       else if (z .eq. prob_hi(3)) then
          val = c_bc(3,2,comp-c_bc_comp+1)
       else
          val = 0.d0
       end if

    else if (comp .eq. Epot_bc_comp) then 
       ! electric potential
       if (x .eq. prob_lo(1)) then
          val = Epot_wall(1,1)
       else if (x .eq. prob_hi(1)) then
          val = Epot_wall(2,1)
       else if (y .eq. prob_lo(2)) then
          val = Epot_wall(1,2)
       else if (y .eq. prob_hi(2)) then
          val = Epot_wall(2,2)
       else if (z .eq. prob_lo(3)) then
          val = Epot_wall(1,3)
       else if (z .eq. prob_hi(3)) then
          val = Epot_wall(2,3)
       else
          val = 0.d0
       end if

    else

       print*,'comp=',comp
       call bl_error("calling inhomogeneous_bc_val_3d with invalid comp")

    end if

  end function inhomogeneous_bc_val_3d

end module inhomogeneous_bc_val_module
