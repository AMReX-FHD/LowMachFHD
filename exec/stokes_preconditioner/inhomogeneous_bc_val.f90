module inhomogeneous_bc_val_module

  use bl_types
  use bl_error_module
  use bc_module
  use probin_module       , only: prob_sol, coeff_mag, ABC_coefs
  use probin_common_module, only: prob_lo, prob_hi
  

  implicit none

  private

  public :: scalar_bc, transport_bc, inhomogeneous_bc_val_2d, inhomogeneous_bc_val_3d

contains

  subroutine scalar_bc(phys_bc, bc_code)

    integer, intent(in   ) :: phys_bc
    integer, intent(inout) :: bc_code(1:num_scal_bc)

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
  
  function inhomogeneous_bc_val_2d(comp,x,y) result(val)

    use bl_constants_module, only: M_PI

    integer        , intent(in   ) :: comp
    real(kind=dp_t), intent(in   ) :: x,y
    real(kind=dp_t)                :: val

    ! local
    real(kind=dp_t) :: velx,vely,Lx,Ly,time, visc_coef

    visc_coef=coeff_mag(2)/coeff_mag(1)

    select case (prob_sol)
    case (1,2)

       ! Taylor vortex
       velx = 1.d0
       vely = 1.d0
       Lx = prob_hi(1)-prob_lo(1)
       Ly = prob_hi(2)-prob_lo(2)
       time = 0.d0

       if (comp .eq. vel_bc_comp) then
          val = velx - 2.d0*exp(-8.d0*M_PI**2*visc_coef*time/(Lx*Ly)) &
               *cos(2.d0*M_PI*(x-velx*time)/Lx)*sin(2.d0*M_PI*(y-vely*time)/Ly)
       else if (comp .eq. vel_bc_comp+1) then
          val = vely + 2.d0*exp(-8.d0*M_PI**2*visc_coef*time/(Lx*Ly)) &
               *sin(2.d0*M_PI*(x-velx*time)/Lx)*cos(2.d0*M_PI*(y-vely*time)/Ly)
       else
          val = 0.d0
       end if

    case (3)
       if (comp .eq. vel_bc_comp) then
          val = exp(-8.0d0*M_PI**2*time)*sin(2.d0*M_PI*x)*sin(2.d0*M_PI*y)
       else if (comp .eq. vel_bc_comp+1) then
          val = exp(-8.0d0*M_PI**2*time)*cos(2.d0*M_PI*x)*cos(2.d0*M_PI*y)
       else
          val = 0.d0
       end if

    case (20)

       ! lid driven cavity moving hi-y wall with +x velocity
       if (y .eq. prob_hi(2) .and. comp .eq. vel_bc_comp) then
          val = sin(M_PI*(x-prob_lo(1))/(prob_hi(1)-prob_lo(1)))
       else
          val = 0.d0
       end if

    case (31,32) ! div not free, for testing accuracy

       if (comp .eq. vel_bc_comp) then
          val = sin(2.d0*M_PI*x)*sin(2.d0*M_PI*y)
       elseif (comp .eq. vel_bc_comp+1) then
          val = cos(2.d0*M_PI*x)*sin(2.d0*M_PI*y)
       else 
          call bl_error('inhomogeneous_bc_val_3d comp should not be greater than  3')
       end if

    case default

       val = 0.d0

    end select

  end function inhomogeneous_bc_val_2d

  function inhomogeneous_bc_val_3d(comp,x,y,z) result(val)

    use bl_constants_module, only: M_PI

    integer        , intent(in   ) :: comp
    real(kind=dp_t), intent(in   ) :: x,y,z
    real(kind=dp_t)                :: val

    ! local
    real(kind=dp_t) :: velx,vely,Lx,Ly,time, visc_coef
    real(kind=dp_t) :: ufac,vfac,pfac,hx,hy,hz,freq

    visc_coef=coeff_mag(2)/coeff_mag(1)

    select case (prob_sol)
    case (1)

       ! Taylor vortex, quasi-2D
       velx = 1.d0
       vely = 1.d0
       Lx = prob_hi(1)-prob_lo(1)
       Ly = prob_hi(2)-prob_lo(2)
       time = 0.d0

       if (comp .eq. vel_bc_comp) then
          val = velx - 2.d0*exp(-8.d0*M_PI**2*visc_coef*time/(Lx*Ly)) &
               *cos(2.d0*M_PI*(x-velx*time)/Lx)*sin(2.d0*M_PI*(y-vely*time)/Ly)
       else if (comp .eq. vel_bc_comp+1) then
          val = vely + 2.d0*exp(-8.d0*M_PI**2*visc_coef*time/(Lx*Ly)) &
               *sin(2.d0*M_PI*(x-velx*time)/Lx)*cos(2.d0*M_PI*(y-vely*time)/Ly)
       else
          val = 0.d0
       end if

    case (2) 

       !ABC flow    
       freq  = 2.d0*M_PI

       ufac = dexp(-freq*freq*visc_coef*time)
       pfac = dexp(-2.0d0*freq*freq*visc_coef*time)

       hx = ABC_coefs(1)
       hy = ABC_coefs(2)
       hz = ABC_coefs(3)

       if (comp .eq. vel_bc_comp) then
          val = 1.0d0 + ufac*(hz*cos(freq*(y-time))+hx*sin(freq*(z-time)))
       else if (comp .eq. vel_bc_comp+1) then
          val = 1.0d0 + ufac*(hy*sin(freq*(x-time))+hx*cos(freq*(z-time)))
       else if (comp .eq. vel_bc_comp+2) then
          val = 1.0d0 + ufac*(hy*cos(freq*(x-time))+hz*sin(freq*(y-time)))
       else 
          call bl_error('inhomogeneous_bc_val_3d comp should not be greater than  3')
       end if

    case (20)

       ! lid driven cavity moving hi-z wall with +x and +y velocity
       if (z .eq. prob_hi(3) .and. (comp .eq. vel_bc_comp .or. comp .eq. vel_bc_comp+1) ) then
          val = sin(M_PI*(x-prob_lo(1))/(prob_hi(1)-prob_lo(1))) &
               *sin(M_PI*(y-prob_lo(2))/(prob_hi(2)-prob_lo(2)))
       else
          val = 0.d0
       end if

    case (31,32) ! div not free, for testing accuracy

       if (comp .eq. vel_bc_comp) then
          val = sin(2.d0*M_PI*x)*sin(2.d0*M_PI*y)*sin(2.d0*M_PI*z)
       else if (comp .eq. vel_bc_comp+1) then
          val = cos(2.d0*M_PI*x)*sin(2.d0*M_PI*y)*sin(2.d0*M_PI*z)
       else if (comp .eq. vel_bc_comp+2) then
          val = sin(2.d0*M_PI*x)*sin(2.d0*M_PI*y)*cos(2.d0*M_PI*z)
       else 
          call bl_error('inhomogeneous_bc_val_3d comp should not be greater than  3')
       end if

    case default

       val = 0.d0

    end select

  end function inhomogeneous_bc_val_3d

end module inhomogeneous_bc_val_module
