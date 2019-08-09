program main

use gmres_module
use fgmres_module
use mat_times_vec
use pre_inv_vec
use working_precision
use prob_size_module                      ! control the problem size: nx, ny


  ! local variables 
  !integer, parameter :: wp=kind(0.0d0)
  real(wp), dimension(:), allocatable, target :: x, true_sol
  real(wp), dimension(:), allocatable, target :: b, y
  integer :: status
  integer :: i, nunk
  real(wp) :: dx,dy, L, rel_diff  

  nunk=nx*ny          ! from prob_size

  allocate(x(nunk),b(nunk), true_sol(nunk))
  !allocate(y(nunk))
  
  !----Step 1, evaluate right hand side------------------
  ! Analytic solution \phi = sin(2\pi x/L)*cos(2\pi y/L)
  ! right hand side, -Lap \phi = 8\pi^2/L^2 *sin(x/L)*cos(y/L)
  dx=(prob_hi_x-prob_lo_x)/nx
  dy=(prob_hi_y-prob_lo_y)/ny
  L=prob_hi_x-prob_lo_x    ! we also assume prob_hi_y-prob_lo_y=L
  do j=1,ny
    do i=1,nx
      indx=(j-1)*nx+i    
      x_coor=prob_lo_x+(i-0.5d0)*dx
      y_coor=prob_lo_y+(j-0.5d0)*dy
      b(indx)=8.0d0*pi**2/(L**2)*sin(2.0d0*pi*x_coor/L)*cos(2.0d0*pi*y_coor/L)
      true_sol(indx)=sin(2.0d0*pi*x_coor/L)*cos(2.0d0*pi*y_coor/L)
!     x(indx)=5.0d0*true_sol(indx)      ! initial solution, just for simplicity
    end do 
  end do
!  x(2)=5.0d0      ! do some perturbations ,for initial solution
!  x(20)=2.5d0     

   call random_number(x)   ! initial solution 
 !------------------------------------------------------- 
! for testing 
!  allocate(y(nunk))
!  call MatVec(true_sol, y)
!  do i=1, nunk
   ! print*, y(i)-b(i)       ! for periodic BC A*true=b
!    print*, (y(i)-b(i))/sqrt(dot_product(b,b))    ! it is better to measure realitve diff
!  end do
!  deallocate(y)
!  stop

!  call ApplyPrecond (b, x)  ! b=P^{-1}x, here P=A but using cholesky fact to invert 
!  do i=1, nunk
!    print*, x(i)-true_sol(i)   ! check P^{-1} is also good enough
!    print*, (x(i)-true_sol(i))/sqrt(dot_product(true_sol,true_sol))  ! it is better to measure realitve diff
!  end do
!  stop

 !-----Step 2, call GMRES code----------------------------
 ! call FGMRES (x, b, MatVec_1D, ApplyPrecond_1D)   ! need to modify a little bit
 ! call GMRES (x, b, MatVec_1D, ApplyPrecond_1D)     
 !  call GMRES (x, b, MatVec, ApplyPrecond)
  call GMRES (x, b, MatVec, ApplyPrecond_1D)      ! using diagonal preconditioner

 !--------------------------------------------------------

 !-----Step 3, normalize the numerical solution ?, not necessary-----------
 !do i=1, nunk
 !  x(i)=x(i)-sum(x)/nunk
 !end do
 !-----Step 3, compare with analytic solution-------------
 print*, 'the difference between the sol and the true sol is constant, to see this uncomment the following'
 do i=1, nunk
!   print*, x(i)
   print*, x(i)-true_sol(i)
 end do
 !--------------------------------------------------------

  write (*, '(a)') 'gmres called'
  deallocate(x,b)
 ! deallocate(y)

end program
