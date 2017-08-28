module inhomogeneous_bc_val_module

  use bl_types
  use bc_module
  use bl_error_module
  use bl_constants_module
  use probin_multispecies_module, only: c_bc, rhotot_bc
  use probin_common_module, only: prob_lo, prob_hi, wallspeed_lo, wallspeed_hi, prob_type, &
                                  nspecies, algorithm_type, rho0
  use probin_charged_module, only: Epot_wall, Epot_wall_bc_type

  implicit none

  private

  public :: scalar_bc, transport_bc, inhomogeneous_bc_val_2d, inhomogeneous_bc_val_3d

contains

  subroutine scalar_bc(phys_bc, bc_code)

    ! variable ordering for charged fluid code (num_scal_bc = 2*nspecies+3):
    ! 1                       = total density
    ! 2:nspecies+1            = concentrations
    ! nspecies+2:2*nspecies+1 = molfrac or massfrac (dimensionless fractions
    ! 2*nspecies+2            = temperature
    ! 2*nspecies+3            = electric potential

    integer, intent(in   ) :: phys_bc
    integer, intent(inout) :: bc_code(1:num_scal_bc)

    if ((phys_bc == NO_SLIP_WALL) .or. (phys_bc == SLIP_WALL)) then

       bc_code(1:num_scal_bc-2) = FOEXTRAP  ! Pure Neumann for total density, conctractions, mol/mass fractions
       bc_code(  num_scal_bc-1) = EXT_DIR   ! But temperature is still specified at the boundary (via call to multifab_coefbc)
       
       if (Epot_wall_bc_type .eq. 1) then
          bc_code(  num_scal_bc  ) = EXT_DIR   ! Electric potential fixed on walls
       else if (Epot_wall_bc_type .eq. 2) then
          bc_code(  num_scal_bc  ) = FOEXTRAP  ! Charge density fixed on walls
       end if

    else if ((phys_bc == NO_SLIP_RESERVOIR) .or. (phys_bc == SLIP_RESERVOIR)) then

       bc_code = EXT_DIR   ! Pure Dirichlet

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
          call bl_error("calling inhomogeneous_bc_val_2d with scal_bc_comp (full rho")
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
       if (x .eq. prob_lo(1)) then
          val = Epot_wall(1,1)
       else if (x .eq. prob_hi(1)) then
          val = Epot_wall(2,1)
       else if (y .eq. prob_lo(2)) then
          val = Epot_wall(1,2)
       else if (y .eq. prob_hi(2)) then
          val = Epot_wall(2,2)
       else
          val = 0.d0
       end if

    else

       print*,'comp=',comp
       call bl_error("calling inhomogeneous_bc_val_2d with invalid comp")

    end if


  end function inhomogeneous_bc_val_2d

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
          call bl_error("calling inhomogeneous_bc_val_3d with scal_bc_comp (full rho")
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
