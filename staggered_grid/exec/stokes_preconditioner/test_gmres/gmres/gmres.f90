!**********************************************************************
! Compute the solution to a nonsymmetric system of linear equations of the
! form Ax=b by the restarted Generalized Minimal Residual method, 
! This code is based on Givens tranform, not Housholder transform, we refer 
! to the Matlab code in 
! http://www.netlib.org/templates/matlab/gmres.m
! and 
! http://people.sc.fsu.edu/~jburkardt/f_src/mgmres/mgmres.html
!**********************************************************************
module gmres_module
  use working_precision
  
   implicit none 
  ! constants
  real(wp), parameter :: tol_rel=1.0e-8        ! tolerance for checking convergence,relative 
  logical, parameter :: verbose = .true.
  integer, parameter :: itr_max=4              ! the max number of outer iterations
  integer, parameter :: max_inner=40           ! the max number of inner iteration, or restart num

contains

  subroutine gmres (x, rhs, MatVec, ApplyPrecond)
    
    interface
      subroutine MatVec(x, y)     ! calculate Ax=y
      use working_precision
       real(wp), intent(in), dimension(:), target :: x
       real(wp), intent(inout), dimension(:), target :: y
      end subroutine MatVec
    end interface

    interface
      subroutine ApplyPrecond (b, x) ! calculate x=P^{-1} b
      use working_precision
        real(wp), intent(in), dimension(:), target :: b
        real(wp), intent(inout), dimension(:), target :: x
      end subroutine ApplyPrecond
    end interface

    ! Arguments:
    optional :: ApplyPrecond
    real(wp), intent(in), dimension(:) :: rhs    ! use multifab
    real(wp), intent(inout), dimension(:) :: x   ! use multifab
  
    ! Local:
    integer :: n, k_copy, i, j, k                ! for looping iteration
    integer :: itr_used, itr                     ! number of iterations used
    real(wp) :: norm_resid, rel_resid, norm_init_resid, mu, av, htmp, error, norm_rhs, normw

    real(wp), allocatable, dimension(:)    :: r,w,v_tmp     ! use multifab
    real(wp), allocatable, dimension(:, :) :: v       ! use multifab
    real(wp), allocatable, dimension(:)    :: g,y     ! for triangular solve
    real(wp), allocatable, dimension(:)    :: c,s     ! for Givens transform
    real(wp), allocatable, dimension(:, :) :: h       ! for saving Hessenberg mat

    n=size(x)
    if (verbose) then
      write ( *, '(a,i4)' ) '  Number of unknowns = ', n
    end if
    allocate(r(n),w(n),v_tmp(n))         ! for storing residual

    itr_used = 0
    norm_rhs=sqrt(dot_product(rhs, rhs))   
    if (norm_rhs .eq. 0.0d0) then 
       print*, 'the RHS is 0, stop the code'
       stop
    end if 
    
    call MatVec(x, r)
    r(1:n)=rhs(1:n)-r(1:n)
    if(present(ApplyPrecond)) then
       call ApplyPrecond(r, w)            ! solve for P^{-1} r
       r(1:n)=w(1:n)
    end if
    norm_resid = sqrt(dot_product(r, r))             ! norm of residual 
    norm_init_resid=norm_resid
  !  rel_resid=norm_resid/norm_rhs
    rel_resid=norm_resid/norm_init_resid
  
    if (rel_resid <= tol_rel) then
       print*, 'convege in 0 iteration, x itself is the solution'
       stop
    end if

    allocate(v(n,max_inner+1))           ! for storing Krylov subspace vec  ! use multifab
    allocate(h(max_inner+1, max_inner))  ! for storing Hessenberg mat
    allocate(c(max_inner))
    allocate(s(max_inner))
    allocate(y(max_inner+1))
    allocate(g(max_inner+1))
    v(1:n,1:max_inner+1) = 0.0d0
    h(1:max_inner+1,1:max_inner) = 0.0d0
    c(1:max_inner)=0.0d0
    s(1:max_inner)=0.0d0
    y(1:max_inner+1)=0.0d0

    do itr = 1, itr_max                    ! outer iteration 
      call MatVec(x, r)                    ! Cal r=Ax
      r(1:n) = rhs(1:n) - r(1:n)           ! r= b-Ax

      if(present(ApplyPrecond)) then
        call ApplyPrecond(r, w)            ! solve for P^{-1} r
        r(1:n)=w(1:n)
      end if

      norm_resid = sqrt(dot_product(r, r))
      if (verbose) then
        write ( *, '(a,i4,a,g14.6)' ) ' Outer ITR = ', itr, ' Preconditioned Residual = ', norm_resid
      end if

      v(1:n,1) = r(1:n)/norm_resid
      g(1) = norm_resid
      g(2:max_inner+1) = 0.0d0

      do k = 1, max_inner
        k_copy = k
        call MatVec(v(1:n,k), v_tmp)     ! v_tmp=A*V(k)
        
        if (present(ApplyPrecond)) then
           call ApplyPrecond(v_tmp, w)     !
        end if

        do j = 1, k
          h(j,k) = dot_product(w, v(1:n,j))
          w(1:n) = w(1:n) - v(1:n,j) * h(j,k)
        end do
        normw = sqrt(dot_product(w, w))
        h(k+1,k) = normw     

        if ( h(k+1,k) /= 0.0d0 ) then
          v(1:n,k+1)=w/h(k+1,k)         ! w/norm(w)
        else
          print*, 'error in orthogonalization'
          stop
        end if 

        call least_squares(k, h, c, s, g)                     ! solve least square problem

 !       rel_resid = abs(g(k+1))/norm_rhs                            ! for checking convergence
         rel_resid = abs(g(k+1))/norm_init_resid 
        itr_used = itr_used + 1

        if (verbose) then
          write ( *, '(a,i4,a,g14.6)' ) ' Inner itr: K = ', k, ' Residual/norm_b = ', rel_resid
        end if
        if (rel_resid .le. tol_rel) then
          exit
        end if
      end do                                              ! end of inner loop

      ! update the solution 
      k = k_copy-1
      call SolveUTriangular(k, h, g, y)                   ! y= H(1:k, 1:k)^{-1} g

      do i = 1, n
        x(i) = x(i) + dot_product(v(i,1:k+1), y(1:k+1))  ! update solution
      end do

      if (rel_resid .le. tol_rel) then
        exit
      end if

    end do                                              ! end of outer loop

    deallocate(h)
    deallocate(v)
    deallocate(g,y) 
    deallocate(c,s)
    deallocate(r,w,v_tmp)
   
    if (verbose ) then
      write ( *, '(a)' ) 'Preconditioned GMRES:'
      write ( *, '(a,i6)' ) '  Iterations = ', itr_used
      write ( *, '(a,g14.6)' ) '  Final relative residual = ', rel_resid
    end if

    return

    contains 
      subroutine least_squares(k, H, c, s, g)
        integer, intent(in) :: k
        real(wp), intent(inout), dimension(:, :) :: H
        real(wp), intent(inout), dimension(:) :: c
        real(wp), intent(inout), dimension(:) :: s
        real(wp), intent(inout), dimension(:) :: g
        ! local variable       
        integer :: i
        real(wp) :: y(1:k+1)
        real(wp) :: temp

        do i = 1,k-1                              ! apply Givens rotation
          temp     =  c(i)*H(i,k) + s(i)*H(i+1,k)
          H(i+1,k) = -s(i)*H(i,k) + c(i)*H(i+1,k)
          H(i,k)   = temp
	end do
        call rotmat(H(k,k), H(k+1,k), c(k), s(k)) ! form i-th rotation matrix

        temp   = c(k)*g(k)                        ! approximate residual norm
        g(k+1) = -s(k)*g(k)
	g(k)   = temp
        H(k,k) = c(k)*H(k,k) + s(k)*H(k+1,k)
        H(k+1,k) = 0.0d0

      end subroutine least_squares


      subroutine rotmat(a, b, c, s )
      ! Compute the Givens rotation matrix parameters for a and b.
      real(wp), intent(in) :: a,b
      real(wp), intent(inout) :: c,s
      !local
      real(wp) :: temp

      if ( b .eq. 0.0d0 ) then
         c = 1.0d0
         s = 0.0d0
      elseif(abs(b) > abs(a)) then
         temp = a/b
         s = 1.0d0/sqrt(1.0d0 + temp**2)
         c = temp * s
      else
         temp = b / a;
         c = 1.0d0/sqrt(1.0d0 + temp**2)
         s = temp * c
      end if
      end subroutine rotmat
      
      subroutine SolveUTriangular(k, h, g, y)
        integer, intent(in) :: k
        real(wp), intent(in), dimension(:, :) :: h
        real(wp), intent(in), dimension(:) :: g
        real(wp), intent(inout), dimension(:)  :: y
        !local 
        integer :: i
        
        y(k+1) = g(k+1)/h(k+1,k+1)
        do i = k, 1, -1
          y(i) = (g(i)-dot_product(h(i,i+1:k+1), y(i+1:k+1)))/h(i,i)
        end do
      end subroutine SolveUTriangular

  end subroutine gmres



end module gmres_module

