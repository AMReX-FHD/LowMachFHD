!**********************************************************************
! User defined subroutine for calculating inverse of precondtioner 
! times a vector
!**********************************************************************
module pre_inv_vec
  use working_precision
  use prob_size_module

  implicit none 
  
  contains 

  subroutine ApplyPrecond_1D (b, x) ! calculate x=P^{-1} b   ! diagonal preconditioner
       real(wp), intent(in), dimension(:), target :: b
       real(wp), intent(inout), dimension(:), target :: x

    !local
    integer :: i, nunk
    real(wp) :: dx, dy

    nunk=size(b)

    dx=(prob_hi_x-prob_lo_x)/nx
    dy=(prob_hi_y-prob_lo_y)/ny

    do i=1, nunk
  !     if (i.eq.1) then
  !        x(i)=b(i)*(dx**2)/4.0d0
  !     elseif(i.eq.nunk) then
  !        x(i)=b(i)*(dx**2)/4.0d0
  !     else
          x(i)=b(i)*(dx**2)/4.0d0
  !     end if
    end do

   ! normalize, because we have mean value 0
   ! do i=1,nunk
   !   x(i)=x(i)-sum(x)/nunk
   ! end do 
   ! print*, 'this subroutine is calledddddddddd'

  end subroutine ApplyPrecond_1D


  ! Cholesky factorization for the system
  subroutine ApplyPrecond (b, x) ! calculate x=P^{-1} b
       real(wp), intent(in), dimension(:), target :: b
       real(wp), intent(inout), dimension(:), target :: x

    !local
    integer :: i, j, k, nunk, rc
    real(wp), allocatable, dimension(:, :) :: mat
    real(wp) :: dx, dy
 
    dx=(prob_hi_x-prob_lo_x)/nx
    dy=(prob_hi_y-prob_lo_y)/ny
    nunk=size(b)

    allocate(mat(nunk,nunk))
    mat=0.0d0
    
    ! initialize 
    x=0.0d0
    
    ! evaluate coeff mat, first part of coeff mat
    do j=1,ny
      do i=1,nx
        do k=1,nx
          if (i.eq.k) then                                                     ! diagonal
            mat((j-1)*nx+i,(j-1)*nx+k)=mat((j-1)*nx+i,(j-1)*nx+k)+2.0d0/(dx**2)  
          elseif((i.eq. k-1).or.(i.eq.k+1)) then                               ! non-diagonal
            mat((j-1)*nx+i,(j-1)*nx+k)=mat((j-1)*nx+i,(j-1)*nx+k)-1.0d0/(dx**2)
          elseif(((i.eq.1).and.(k.eq.nx)).or.((i.eq.nx).and.(k.eq.1))) then    ! circulant structure
            mat((j-1)*nx+i,(j-1)*nx+k)=mat((j-1)*nx+i,(j-1)*nx+k)-1.0d0/(dx**2)
          end if
        end do
      end do 
    end do
    
    ! second part of coeff mat
    do j=1,ny
      do k=1,ny
        do i=1,nx
          if (j.eq.k) then                                                       ! diagonal
            mat((j-1)*nx+i,(k-1)*nx+i)=mat((j-1)*nx+i,(k-1)*nx+i)+2.0d0/(dy**2)  
          elseif((j.eq.k-1).or.(j.eq.k+1)) then                                  ! non-diagonal
            mat((j-1)*nx+i,(k-1)*nx+i)=mat((j-1)*nx+i,(k-1)*nx+i)-1.0d0/(dy**2)
          elseif(((j.eq.1).and.(k.eq.nx)) .or. ((j.eq.nx).and.(k.eq.1))) then    ! circulant structure
            mat((j-1)*nx+i,(k-1)*nx+i)=mat((j-1)*nx+i,(k-1)*nx+i)-1.0d0/(dy**2)
          end if
        end do
      end do 
    end do

!    do i=1, nunk   ! check mat is correct
!      do j=1, nunk
!        print*, mat(i,j)
!      end do
!    end do 
!    stop 

    ! rank 1 correction to avoid singularity
    mat(nunk,nunk)=mat(nunk,nunk)+1.0d0/(dy**2)

    ! use chol_dec, 
    call chol_dec(nunk, mat, rc)
    !print*, 'if the chol fact works? rc should be 0'    
    call chol_sol(nunk, mat, b, x, rc)
    print*, 'if Cholesky works properly, rc should be 0'
    print*, rc

    ! normalize, because we have mean value 0
    ! do i=1,nunk
    !   x(i)=x(i)-sum(x)/nunk
    ! end do 
    deallocate(mat)
    return

  contains 

    subroutine chol_dec(n, a, rc)
    !======================================================================
    !*  chol_dec decomposes the symmetric positive definite matrix mat.     *
    !*   Input parameters:                                                *
    !*   ================                                                 *
    !*      n        integer (n > 0)  Dimension of matrix a               *
    !*      a        REAL matrix (n,n)                                    *
    !*               Matrix of left hand coefficients                     *
    !*                                                                    *
    !*   Output parameter:                                                *
    !*   ================                                                 *
    !*      a        REAL matix (n,n)                                     *
    !*               Cholesky decomposition in the lower triangle         *
    !*                                                                    *
    !*   Return value rc:                                                 *
    !*      = 0      all ok                                               *
    !*      = 1      n < 1                                                *
    !*      = 2      Matrix not  positive definite                        *
    !*====================================================================*
    !*   Functions in use:   dsqrt (square root in double precision)      *
    !*====================================================================*
    !*                                                                    *
    !*   Constants in use:   EPSQUAD (small number)                       *
    !======================================================================
    integer, intent(in) :: n
    integer, intent(inout) :: rc
    real(wp), intent(inout) ::  a(0:n-1,0:n-1)
    !local 
    real(wp) :: EPSQUAD, sum_val
    integer :: i, j, k
    

    EPSQUAD = 1.0d-14

    if (n < 1) then
      rc=1                           ! n < 1  error
      return
    end if

    if (a(0,0) < EPSQUAD) then       ! matrix a not positive definite
      rc=2
      return
    end if
    a(0,0) = dsqrt(a(0,0))
    do j = 1, n-1
      a(j,0) = a(j,0) / a(0,0)
    end do

    do i = 1, n-1
      sum_val = a(i,i)
      do j = 0, i-1
       sum_val = sum_val - a(i,j)*a(i,j)
      end do
      if (sum_val < EPSQUAD) then         ! matrix a not positive definite
        rc=2
        return
      end if
      a(i,i) = dsqrt(sum_val)
      do j = i+1, n-1
        sum_val = a(j,i)
        do k = 0, i-1
          sum_val = sum_val - a(i,k) * a(j,k)
        end do
        a(j,i) = sum_val / a(i,i)
      end do
    end do
    rc = 0  ! all OK
    return
  end subroutine chol_dec

  subroutine chol_sol(n, a, b, x, rc)
  !======================================================================
  !*  chol_sol finds the solution x of the linear system B' *  B * x = b  *
  !*  for a lower triangular nonsingular matrix B as supplied in chol_dec.*  
  !*   Input parameters:                                                *
  !*      n        integer  (n > 0)                                     *
  !*               Dimension of matrix a, size of vectors b and x       *
  !*      a        REAL  matrix (n,n)                                   *
  !*               lower triangular matrix as supplied by  chol_dec       *
  !*      b        REAL  vector (n)   Right hand side                   *
  !*   Output parameter:                                                *
  !*      x        REAL  vector (n)   solution vector                   *
  !*   Return value rc:                                                 *
  !*      = 0      all ok                                               *
  !*      = 1      improper lwer triangular matrix or  n < 1            *
  !======================================================================
  integer, intent(in) :: n
  integer, intent(inout) :: rc
  real(wp), intent(in) :: a(0:n-1,0:n-1),b(0:n-1)
  real(wp), intent(inout) :: x(0:n-1)
  real(wp) :: ZERO
  ! local 
  real(wp) :: sum_val
  integer :: i, j, k
  ZERO = 0.d0

  if (n < 1) then   ! n < 1 error
    rc=1
    return
  end if

  if (a(0,0).eq.ZERO) then               ! improper factor matrix
    rc=1
    return
  end if

  x(0) = b(0) / a(0,0)                   ! update right hand side
  do k = 1, n-1
    sum_val = ZERO
    do j = 0, k-1
      sum_val = sum_val + a(k,j) * x(j)
    end do
    if (a(k,k).eq.ZERO) then
      rc=1
      return
    end if
    x(k) = (b(k) - sum_val) / a(k,k)
  end do

  x(n-1) = x(n-1) / a(n-1,n-1)            ! back substitution
  do k=n-2, 0, -1
    sum_val = ZERO
    do j = k+1, n-1
      sum_val = sum_val + a(j,k) * x(j)
    end do
    x(k) = (x(k) - sum_val) / a(k,k)
  end do

  rc = 0  !all ok
  return

  end subroutine chol_sol
!=======================end of cholesky part===========================

  end subroutine ApplyPrecond


end module pre_inv_vec



