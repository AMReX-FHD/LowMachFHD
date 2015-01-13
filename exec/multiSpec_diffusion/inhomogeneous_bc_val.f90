module inhomogeneous_bc_val_module

  use bl_types
  use bc_module
  use bl_error_module
  use probin_multispecies_module, only: nspecies, rho_bc, rho_part_bc_comp
  use probin_common_module, only: prob_lo, prob_hi

  implicit none

  private

  public :: scalar_bc, transport_bc, inhomogeneous_bc_val_2d, inhomogeneous_bc_val_3d

contains

  subroutine scalar_bc(phys_bc, bc_code)

    integer, intent(in   ) :: phys_bc
    integer, intent(inout) :: bc_code(1:num_scal_bc)

    if ((phys_bc == NO_SLIP_WALL) .or. (phys_bc == SLIP_WALL)) then

       bc_code(1:num_scal_bc-1) = FOEXTRAP  ! Pure Neumann for densities / mass fractions
       bc_code(num_scal_bc) = EXT_DIR ! But temperature is still specified at the boundary (via call to multifab_coefbc)

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

    if (comp .lt. rho_part_bc_comp .or. &
        comp .gt. rho_part_bc_comp+nspecies-1) then

       print*,'comp=',comp
       call bl_error("calling inhomogeneous_bc_val_2d with invalid comp")

    else
       
       ! rho_i boundary condition
       if (x .eq. prob_lo(1)) then
          val = rho_bc(1,1,comp-rho_part_bc_comp+1)
       else if (x .eq. prob_hi(1)) then
          val = rho_bc(1,2,comp-rho_part_bc_comp+1)
       else if (y .eq. prob_lo(2)) then
          val = rho_bc(2,1,comp-rho_part_bc_comp+1)
          !write(*,*) "bottom:", val, comp-rho_part_bc_comp+1
       else if (y .eq. prob_hi(2)) then
          val = rho_bc(2,2,comp-rho_part_bc_comp+1)
          !write(*,*) "top:", val, comp-rho_part_bc_comp+1
       else
          val = 0.d0
       end if

    end if

  end function inhomogeneous_bc_val_2d

  function inhomogeneous_bc_val_3d(comp,x,y,z,time_in) result(val)

    integer        , intent(in   ) :: comp
    real(kind=dp_t), intent(in   ) :: x,y,z
    real(kind=dp_t), intent(in), optional :: time_in
    real(kind=dp_t)                :: val

    if (comp .lt. rho_part_bc_comp .or. &
        comp .gt. rho_part_bc_comp+nspecies-1) then

       print*,'comp=',comp
       call bl_error("calling inhomogeneous_bc_val_3d with invalid comp")

    else
       
       ! rho_i boundary condition
       if (x .eq. prob_lo(1)) then
          val = rho_bc(1,1,comp-rho_part_bc_comp+1)
       else if (x .eq. prob_hi(1)) then
          val = rho_bc(1,2,comp-rho_part_bc_comp+1)
       else if (y .eq. prob_lo(2)) then
          val = rho_bc(2,1,comp-rho_part_bc_comp+1)
       else if (y .eq. prob_hi(2)) then
          val = rho_bc(2,2,comp-rho_part_bc_comp+1)
       else if (z .eq. prob_lo(3)) then
          val = rho_bc(3,1,comp-rho_part_bc_comp+1)
       else if (z .eq. prob_hi(3)) then
          val = rho_bc(3,2,comp-rho_part_bc_comp+1)
       else
          val = 0.d0
       end if

    end if

  end function inhomogeneous_bc_val_3d

end module inhomogeneous_bc_val_module
