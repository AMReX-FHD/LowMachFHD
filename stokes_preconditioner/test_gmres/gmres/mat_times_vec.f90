!**********************************************************************
! User defined subroutine for calculating matrix times a vector
!**********************************************************************

module mat_times_vec
  use working_precision
  use prob_size_module

  implicit none 
  
  contains 

  ! A simple 1D case, artificial
  subroutine MatVec_1D(x, y)     ! calculate Ax=y
    real(wp), intent(in), dimension(:), target :: x
    real(wp), intent(inout), dimension(:), target :: y
!    integer, intent(out) :: status
    
    !local
    integer :: i, nunk
    nunk=size(x)

    do i=1, nunk, 1
       if (i.eq.1) then
          y(i)=4.0d0*x(i)-1.0d0*x(i+1)
       elseif(i.eq.nunk) then
          y(i)=4.0d0*x(i)-1.0d0*x(i-1)
       else
          y(i)=4.0d0*x(i)-1.0d0*x(i-1)-1.0d0*x(i+1)
       end if
    end do

  end subroutine MatVec_1D
!---------------------------------------------------------------------------


  ! 2D case. Periodic BC, Poisson problem
  subroutine MatVec(x, y)     ! calculate Ax=y
  real(wp), intent(in), dimension(:), target :: x
  real(wp), intent(inout), dimension(:), target :: y
! integer, intent(out) :: status
    
  !local
  integer :: i,j,nunk
  real(wp) :: dx, dy
 
!  dx=1.0d0
!  dy=1.0d0
  dx=(prob_hi_x-prob_lo_x)/nx
  dy=(prob_hi_y-prob_lo_y)/ny
  nunk=size(x)
  
  !initialized to 0  
  y=0.0d0

  ! second order derivative---------x component 
  do j=1,ny
    do i=1, nx
      if (i.eq.1) then
        y((j-1)*nx+i)=y((j-1)*nx+i)+(2.0d0*x((j-1)*nx+i)-1.0d0*x((j-1)*nx+i+1)-1.0d0*x((j-1)*nx+nx))/(dx**2)
      elseif(i.eq.nx) then
        y((j-1)*nx+i)=y((j-1)*nx+i)+(2.0d0*x((j-1)*nx+i)-1.0d0*x((j-1)*nx+i-1)-1.0d0*x((j-1)*nx+1))/(dx**2)
      else
        y((j-1)*nx+i)=y((j-1)*nx+i)+(2.0d0*x((j-1)*nx+i)-1.0d0*x((j-1)*nx+i-1)-1.0d0*x((j-1)*nx+i+1))/(dx**2)
      end if
    end do
  end do

  ! second order derivative---------y component 
  do j=1,ny
    do i=1, nx
      if (j.eq.1) then
        y((j-1)*nx+i)=y((j-1)*nx+i)+(2.0d0*x((j-1)*nx+i)-1.0d0*x(j*nx+i)-1.0d0*x((ny-1)*nx+i))/(dx**2)
      elseif(j.eq.ny) then
        y((j-1)*nx+i)=y((j-1)*nx+i)+(2.0d0*x((j-1)*nx+i)-1.0d0*x((1-1)*nx+i)-1.0d0*x((j-2)*nx+i))/(dx**2)
      else
        y((j-1)*nx+i)=y((j-1)*nx+i)+(2.0d0*x((j-1)*nx+i)-1.0d0*x(j*nx+i)-1.0d0*x((j-2)*nx+i))/(dx**2)
      end if
    end do
  end do

  ! Alternative way of implementation can be first reshape x according mesh grid, then take finite difference.

  end subroutine MatVec


end module mat_times_vec
