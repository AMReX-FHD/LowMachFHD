!**********************************************************************
! Compute the solution to a nonsymmetric system of linear equations of the
! form Ax=b by the restarted Flexible Generalized Minimal Residual method, 
! This code is based on Givens tranform, not Housholder transform, we refer 
! to the paper by VAL´ERIE FRAYSS´E
! Kvasar Technology LLC  LUC GIRAUD  ENSEEIHT-IRIT
! and  SERGE GRATTON  CNES
!**********************************************************************
module fgmres_module
  use working_precision

  implicit none 
  ! constants
!  integer, parameter :: wp=kind(0.0d0)
  real(wp), parameter :: tol_rel=1.0e-6         ! tolerance for checking convergence,relative 
  logical, parameter :: verbose = .true.        ! whether the code print out more information 
  integer, parameter :: itr_max=4               ! the max number of outer iterations
  integer, parameter :: max_inner=40             ! the max number of inner iteration, or restart num

contains

  subroutine FGMRES (x, rhs, MatVec, ApplyPrecond)

    interface
      subroutine MatVec(x, y)        ! calculate Ax=y
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
    real(wp), intent(in), dimension(:) :: rhs             ! The RHS
    real(wp), intent(inout), dimension(:) :: x            ! solution and initial solution  
  
    ! Local:
    integer :: n, k_copy, i, j, k,k_tmp                     ! for looping iteration
    integer :: itr_used, itr                          ! number of iteration used
    real(wp) :: rho, rho_tol, norm_rhs, normw, beta
    real(wp) :: h1, h2, mu
    real(wp) :: g1, g2  

    real(wp), allocatable, dimension(:, :) :: H
    real(wp), allocatable, dimension(:, :) :: Z, V
    real(wp), allocatable, dimension(:) :: v_tmp
    real(wp), allocatable, dimension(:) :: g, y
    real(wp), allocatable, dimension(:) :: r, w
    real(wp), allocatable, dimension(:) :: c, s
    !real(wp), allocatable, dimension(:) :: additive, y_tmp
   
    print*, 'solve 2D Press Poisson equation with pure periodic BC, Cholesky factorization as preconditioner'
    n=size(x)
    if (verbose) then
      write ( *, '(a,i4)' ) '  Number of unknowns = ', n
    end if

    itr_used = 0
    norm_rhs=sqrt(dot_product(rhs, rhs))   
    if (norm_rhs .eq. 0.0) then 
      print*, 'the RHS is 0, stop the code'
      stop
    end if 
    rho_tol = norm_rhs*tol_rel

    allocate(r(n))                           ! for storing residual
    allocate(w(n))                           ! for storing temp vector
    allocate(Z(max_inner, n))                ! for storing preconditioned vectors
    allocate(H(max_inner+1, max_inner))      ! for storing Hessenberg mat

    allocate(c(max_inner+1))                 ! for Givens transform
    allocate(s(max_inner+1))                 ! for Givens transform
    allocate(y(max_inner+1))                 ! for solution of least square
    allocate(g(max_inner+1))
    
   ! allocate(v(n))
    allocate(V(n,max_inner+1))
    allocate(v_tmp(n))
   ! v=r/sqrt(dot_product(r, r))

    do itr = 1, itr_max                         ! outer iteration
      call MatVec(x, r)
      r(1:n)=rhs(1:n)-r(1:n)                    ! Calculate the residual 
      
      rho=sqrt(dot_product(r, r))
      V(1:n,1) = r(1:n)/rho
      g(1) = rho
      g(2:max_inner+1) = 0.0D+00
      ! if we do not need to use restart, we can initialize the following outside the outer iteraton
      Z(1:max_inner, 1:n)=0.0D+00 
      H(1:max_inner+1,1:max_inner) = 0.0D+00
      c(1:max_inner+1)=0.0D+00
      s(1:max_inner+1)=0.0D+00
      y(1:max_inner+1)=0.0D+00

      do k = 1, max_inner
        k_copy = k

        if (present(ApplyPrecond)) then
          call ApplyPrecond(V(1:n, k), v_tmp)     !
        end if
        Z(k, 1:n)=v_tmp

        call MatVec(v_tmp, w)
        do j = 1, k                               ! Calculate H
          H(j,k) = dot_product(w, V(1:n,j))
          w = w - V(1:n,j) * H(j,k)
        end do
        normw = sqrt(dot_product(w, w))
        H(k+1,k) = normw      

        if ( H(k+1,k) /= 0.0D+00 ) then           ! Cal the basis V(k+1)
          V(1:n,k+1)=w/H(k+1,k)                   ! w/norm(w)
        end if 
             
        ! solve the least square problem 
        call least_squares(k, H, c, s, g)

        ! check convergence
        rho = abs(g(k+1))
        itr_used = itr_used + 1

        if (verbose) then
          write ( *, '(a,i4,a,i4)' ) ' Outer_iter = ', itr, '  Inner iter = ', k
         ! write ( *, '(a,i4,a,g14.6)' ) '  K = ', k, '  Residual = ', rho
          write ( *, '(a,i4,a,g14.6)' ) '  K = ', k, '  Relative Residual = ', rho/norm_rhs
        end if
   
        if ( rho <= rho_tol ) then
          exit
        end if
      end do                                              ! end of inner loop

      ! solve the triangular system H(k+1, k+1) y = g, use do loop     
      k = k_copy -1
      call SolveUTriangular(k, H, g, y)

      do i = 1, n                                  ! update solution
        x(i) = x(i) + dot_product(Z(1:k+1, i), y(1:k+1))  
      end do
      
      if (rho<=rho_tol ) then
        exit
      end if

    end do                                         ! end of outer loop

    deallocate(c,s)          ! for Givens transform
    deallocate(Z)          ! for storing preconditioned vectors
    deallocate(H)          ! for storing Hessenberg mat
    deallocate(g,y)
    deallocate(r,w)
    deallocate(V)
    deallocate(v_tmp)

    if (verbose ) then
 !     write ( *, '(a)' ) 'Preconditioned GMRES:'
      write ( *, '(a,i6)' ) '  Iterations = ', itr_used
      write ( *, '(a,g14.6)' ) '  Final relative residual = ', rho/norm_rhs
    end if

  contains
    subroutine least_squares(k, H, c, s, g)
      integer, intent(in) :: k
      real(wp), intent(inout), dimension(:, :) :: H
      real(wp), intent(inout), dimension(:) :: c
      real(wp), intent(inout), dimension(:) :: s
      real(wp), intent(inout), dimension(:) :: g
    
      !local
      integer :: j
      real(wp) :: y(1:k+1)
      real(wp) :: mu, g1, g2 
     ! print*, H(k+1, k)

      if ( 1 < k ) then                         ! apply Givens rotation to H
        y(1:k+1) = H(1:k+1,k)                   ! y is temparily used to do Givens transform 
        do j = 1, k - 1
          call mult_givens(c(j), s(j), j, y )
        end do
        H(1:k+1,k) = y(1:k+1)
      end if

      mu = sqrt(H(k,k)**2 + H(k+1,k)**2 )
      c(k) = H(k,k)/mu
      s(k) = H(k+1,k)/mu
      g1= c(k) * H(k,k)+ s(k) * H(k+1,k)
      g2= -s(k) * H(k,k)+ c(k) * H(k+1,k)
      H(k,k) = g1
      H(k+1, k)=g2
      !H(k+1,k) =0.0D+00
        
      ! apply Givens rotation to g
      call mult_givens (c(k), s(k), k, g )
    end subroutine least_squares

    subroutine SolveUTriangular(k, H, g, y)
      integer, intent(in) :: k
      real(wp), intent(in), dimension(:, :) :: H
      real(wp), intent(in), dimension(:) :: g
      real(wp), intent(inout), dimension(:)  :: y
      !local 
      integer :: i
    
      y(k+1) = g(k+1)/H(k+1,k+1)
      do i = k, 1, -1
        y(i) = (g(i)-dot_product(H(i,i+1:k+1), y(i+1:k+1)))/H(i,i)
      end do
    end subroutine SolveUTriangular

  end subroutine FGMRES

!*****************************************************************************80
  subroutine mult_givens (c, s, k, g)
!! MULT_GIVENS applies a Givens rotation to two successive entries of a vector.
!  Discussion:
!    Input, real (wp) C, S, the cosine and sine of a Givens rotation.
!    Input, integer, K, indicates the location of the first
!    vector entry.
!    Input/output, real(wp) G(1:K+1), the vector to be modified.
!    On output, the Givens rotation has been applied to entries G(K) and G(K+1).

    implicit none

    integer, intent(in) :: k
    real(wp), intent(in) :: c
    real(wp), intent(in) :: s
    real(wp), intent(inout) :: g(1:k+1)
    !local
    real(wp) :: g1
    real(wp) :: g2

    g1 = c * g(k) + s * g(k+1)
    g2 =-s * g(k) + c * g(k+1)
    g(k)   = g1
    g(k+1) = g2
  end subroutine mult_givens
end module fgmres_module

