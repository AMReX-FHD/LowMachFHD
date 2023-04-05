module matrix_utilities

  use bl_types
  use bl_prof_module

  implicit none

  private
  
  integer, parameter, private :: dp=dp_t, wp=dp_t ! Working precision

  public :: Dbar2chi_iterative, choldc, mat_inv, mat_svd, mat_chol
  
  interface mat_svd
     module procedure DGESVD_F95
  end interface   
  
  interface mat_chol
     module procedure DPOTRF_F95
  end interface

  INTERFACE ! F77 LAPACK routines needed

      SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT,     &
     &                   LDVT, WORK, LWORK, INFO )
         IMPORT
         CHARACTER(LEN=1), INTENT(IN) :: JOBU, JOBVT
         INTEGER, INTENT(IN) :: M, N, LDA, LDU, LDVT, LWORK
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: S(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: U(LDU,*), VT(LDVT,*), WORK(*)
      END SUBROUTINE DGESVD

      SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
         IMPORT
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
      END SUBROUTINE DPOTRF
   
      FUNCTION DLANSY( NORM, UPLO, N, A, LDA, WORK )
         IMPORT
         REAL(WP) :: DLANSY
         CHARACTER(LEN=1), INTENT(IN) :: NORM, UPLO
         INTEGER, INTENT(IN) :: LDA, N
         REAL(WP), INTENT(IN) :: A( LDA, * )
         REAL(WP), INTENT(OUT) :: WORK( * )
      END FUNCTION DLANSY

       SUBROUTINE DPOCON( UPLO, N, A, LDA, ANORM, RCOND, WORK, IWORK,   &
     &                    INFO )
         IMPORT
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: ANORM
         REAL(WP), INTENT(OUT) :: RCOND
         INTEGER, INTENT(OUT) :: IWORK( * )
         REAL(WP), INTENT(IN) :: A( LDA, * )
         REAL(WP), INTENT(OUT) :: WORK( * )
      END SUBROUTINE DPOCON

   END INTERFACE
  
contains

    ! nspecies_local is number of species
    ! num_iterations is the number of terms in the sum to use: 3-5 are reasonable values
    ! D_bar is matrix of Maxwell-Stefan binary diffusion coefficient
    ! chi is the multispecies diffusion matrix
    ! Xk is mole fractions --- MUST NOT BE ZERO
    subroutine Dbar2chi_iterative(nspecies_local,num_iterations,D_bar,Xk,molmass_local,chi)
      integer, intent(in) :: nspecies_local
      integer, intent(in) :: num_iterations
      real(kind=wp), intent(in) :: D_bar(1:nspecies_local,1:nspecies_local)
      real(kind=wp), intent(in) :: Xk(1:nspecies_local), molmass_local(1:nspecies_local)
      real(kind=wp), intent(out) :: chi(1:nspecies_local,1:nspecies_local)
      
      ! Local variables
      real(kind=wp) :: term1, term2, MWmix
      real(kind=wp) :: Di(1:nspecies_local)
      real(kind=wp) :: Deltamat(1:nspecies_local,1:nspecies_local), Zmat(1:nspecies_local,1:nspecies_local)
      real(kind=wp), dimension(1:nspecies_local,1:nspecies_local) :: Pmat, Jmat
      real(kind=wp), dimension(1:nspecies_local) :: Minv, Mmat
      real(kind=wp), dimension(1:nspecies_local,1:nspecies_local) :: PJ, matrix1, matrix2
      real(kind=wp) :: scr
      real(kind=wp) :: Ykp(1:nspecies_local), Xkp(1:nspecies_local)

      integer :: i, j, k, ii, jj

      type(bl_prof_timer), save :: bpt

      call build(bpt,"Dbar2chi_iterative")
      
      ! mole fractions correction
      ! Turned this off since it should be done in the caller
      do ii = 1, nspecies_local
       Xkp(ii) = Xk(ii)
      end do

      ! molecular weight of mixture - EGLIB
      Mwmix = 0.0d0
      do ii = 1, nspecies_local
       MWmix = MWmix + Xkp(ii)*molmass_local(ii)
      end do

      ! mass fractions correction - EGLIB
      do ii = 1, nspecies_local
       Ykp(ii) = molmass_local(ii)/MWmix*Xkp(ii)
      end do

      ! Find Di matrix 
      do i = 1, nspecies_local
       term2 = 0.0d0
       do j = 1, nspecies_local
        if(j.ne.i) then
          term2 = term2 + Xkp(j)/D_bar(i,j)
        end if
       end do   
       Di(i) = (1.d0-Ykp(i))/term2 
      end do   

      ! Compute Mmat and Minv
      do i = 1, nspecies_local
       Mmat(i) = Xkp(i)/Di(i)
       Minv(i) = Di(i)/Xkp(i)
      end do
      

      ! Compute P matrix
      Pmat = 0.0d0
      do i = 1, nspecies_local
       do j = 1, nspecies_local
         Pmat(i,j) = - Ykp(j) 
         if(i.eq.j) then
          Pmat(i,j) =  Pmat(i,j) + 1.0d0  
         end if
       end do
      end do

      ! Compute Deltamat
      Deltamat = 0.0d0 
      do i = 1, nspecies_local
       do j = 1, nspecies_local
         if(i.eq.j) then
          term1 = 0.0d0
          do k = 1, nspecies_local
           if(k.ne.i) then
            term1 = term1 + Xkp(i)*Xkp(k)/D_bar(i,k)
           end if
          end do  
          Deltamat(i,i) = term1
         else
          Deltamat(i,j) = -Xkp(i)*Xkp(j)/D_bar(i,j) 
         end if  
          Zmat(i,j) = -Deltamat(i,j)
       end do
      end do  

      ! Compute Zmat
      do i = 1, nspecies_local
        Zmat(i,i) = Zmat(i,i) + Mmat(i)
      end do  

      ! Compute Jmat
      do i = 1, nspecies_local
       do j = 1, nspecies_local
         Jmat(i,j) = Minv(i)*Zmat(i,j)
        end do
       end do

      ! Compute PJ
      PJ = 0.0d0
      do i = 1, nspecies_local
       do j = 1, nspecies_local
        do k = 1, nspecies_local
         PJ(i,j) = PJ(i,j) + Pmat(i,k)*Jmat(k,j)
        end do
       end do
      end do

      ! Compute P M^-1 Pt; store it in matrix2
      do i = 1, nspecies_local
       do j = 1, nspecies_local
        scr = 0.d0
        do k = 1, nspecies_local
         scr = scr + Pmat(i,k)*Minv(k)*Pmat(j,k) 
            ! notice the change in indices for Pmat to represent Pmat^t
        end do
         matrix2(i,j) = scr
         chi(i,j) = scr
       end do
      end do


      do jj = 1,num_iterations
       do i = 1, nspecies_local
        do j = 1, nspecies_local
         scr = 0.d0
         do k = 1, nspecies_local
            scr = scr + PJ(i,k)*chi(k,j)
         end do
          matrix1(i,j) = scr+matrix2(i,j)
        end do
       end do 
       chi=matrix1
      end do

      call destroy(bpt)

  end subroutine

   ! a is input matrix.  
   ! upon return the lower triangle and diagonal are overwritten by the cholesky factor
   subroutine choldc(a,np)
       integer :: np
       real(kind=wp), intent(inout) :: a(np,np)

       real(kind=wp) :: p(np), dij(np,np)
       real(kind=wp) :: dd(np,np)
       real(kind=wp) :: yy(np), mwmix 

       integer :: i, j, k, ii, jj
       real(kind=wp) :: sum1
       real(kind=wp) :: small_number = 0.0d0 ! Some tolerance

       integer :: idiag,ising

       type(bl_prof_timer), save :: bpt

       call build(bpt,"choldc")

       do i = 1, np

           ising = 0

        do j = i, np

           sum1 = a(i,j)

           do k = i-1, 1, -1

              sum1 = sum1 - a(i,k)*a(j,k)

           end do

           if(i.eq.j) then

             if(sum1.le.small_number) then

             p(i) = 0.d0

             ising = 1

             else

             p(i) = sqrt(sum1)

             end if

           else

             if(ising.eq.0)then

                a(j,i) = sum1/p(i)

             else

                a(j,i) = 0.d0

             end if

           end if

        end do

       end do


       do i = 1, np

          do j = i+1, np

           a(i,j) = 0.0d0 ! Zero upper triangle

          end do
          
          a(i,i) = p(i)

       end do

       call destroy(bpt)

    end subroutine

!-------------------------------------------------
! Code adopted by Aleks Donev from 
! https://fortranwiki.org/fortran/show/Matrix+inversion
!-------------------------------------------------

function mat_inv(A) result(Ainv)
  real(wp), dimension(:,:), intent(in) :: A
  real(wp), dimension(size(A,1),size(A,2)) :: Ainv
  
  integer :: n
  
  if(size(A,1)/=size(A,2)) then
    stop "Matrix must be square"
  end if
  
  n = size(A,1)
  select case(n)
  case(1)
    Ainv(1,1)=1.0_wp/A(1,1)
  case(2)
    Ainv=matinv2(A)
  case(3)
    Ainv=matinv3(A)
  case(4)
    Ainv=matinv4(A)
  case default
    Ainv=matinv2(A)
  end select      

end function
  
! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
function inv(A) result(Ainv_wp)
  real(wp), dimension(:,:), intent(in) :: A
  real(wp), dimension(size(A,1),size(A,2)) :: Ainv_wp
  
  real(dp), dimension(size(A,1),size(A,2)) :: Ainv ! Use LAPACK double precision routines

  real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
  
  Ainv_wp = Ainv ! Convert precisions if needed
  
end function inv

  pure function matinv2(A) result(B)
    !! Performs a direct calculation of the inverse of a 2×2 matrix.
    real(wp), intent(in) :: A(2,2)   !! Matrix
    real(wp)             :: B(2,2)   !! Inverse matrix
    real(wp)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * A(2,2)
    B(2,1) = -detinv * A(2,1)
    B(1,2) = -detinv * A(1,2)
    B(2,2) = +detinv * A(1,1)
  end function

  pure function matinv3(A) result(B)
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
    real(wp), intent(in) :: A(3,3)   !! Matrix
    real(wp)             :: B(3,3)   !! Inverse matrix
    real(wp)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
              - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
              + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
  end function

  pure function matinv4(A) result(B)
    !! Performs a direct calculation of the inverse of a 4×4 matrix.
    real(wp), intent(in) :: A(4,4)   !! Matrix
    real(wp)             :: B(4,4)   !! Inverse matrix
    real(wp)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = &
      1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
       - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
       + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
       - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

    ! Calculate the inverse of the matrix
    B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
    B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
    B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
    B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
    B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
    B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
    B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
    B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
    B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
    B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
    B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
    B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
  end function

!-----------------------------------------------
! A. Donev adopted from LAPACK90 interface
!-----------------------------------------------
SUBROUTINE DGESVD_F95( A, S, U, VT, WW, JOB, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   !USE LA_PRECISION, ONLY: WP => DP
   !USE LA_AUXMOD, ONLY: ERINFO, LSAME
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE      
!  .. SCALAR ARGUMENTS ..
   CHARACTER(LEN=1), OPTIONAL, INTENT(IN) :: JOB
   INTEGER, INTENT(OUT), OPTIONAL :: INFO
!  .. ARRAY ARGUMENTS ..
   REAL(WP), INTENT(INOUT) :: A(:,:)
   REAL(WP), INTENT(OUT) :: S(:)
   REAL(WP), INTENT(OUT), OPTIONAL :: WW(:)
   REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: U(:,:), VT(:,:)
!----------------------------------------------------------------------
! 
!  Purpose
!  =======
!  
!       LA_GESVD and LA_GESDD compute the singular values and, 
! optionally, the left and/or right singular vectors from the singular
! value decomposition  (SVD) of a real or complex m by n matrix A. The
! SVD of A is written
!                    A = U * SIGMA * V^H
! where SIGMA is an  m by n matrix which is zero except for its 
! min(m, n) diagonal elements, U is an m by m orthogonal (unitary) 
! matrix, and V is an n by n orthogonal (unitary) matrix. The diagonal
! elements of SIGMA , i.e., the values 
! 
!      sigma(i)= SIGMA(i,i), i = 1, 2,..., min(m, n)
! are the singular values of A; they are real and non-negative, and are
! returned in descending order. The first min(m, n) columns of U and V 
! are the left and right singular vectors of A, respectively.
! LA_GESDD solves the same problem as LA_GESVD but uses a divide and 
! conquer method if singular vectors are desired. For large matrices it
! is usually much faster than LA_GESVD when singular vectors are 
! desired, but uses more workspace.
! 
! Note: The routine returns V^H , not V .
! 
! ========
! 
!    SUBROUTINE LA_GESVD / LA_GESDD( A, S, U=u, VT=vt, &
!              WW=ww, JOB=job, INFO=info )  
!      <type>(<wp>), INTENT(INOUT) :: A(:,:)
!      REAL(<wp>), INTENT(OUT) :: S(:)
!      <type>(<wp>), INTENT(OUT), OPTIONAL :: U(:,:), VT(:,:)
!      REAL(<wp>), INTENT(OUT), OPTIONAL :: WW(:)
!      CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: JOB
!      INTEGER, INTENT(OUT), OPTIONAL :: INFO
!      where
!      <type> ::= REAL | COMPLEX
!      <wp>   ::= KIND(1.0) | KIND(1.0D0)
!  
! Arguments
! =========
! 
! A      (input/output) REAL or COMPLEX array, shape (:, :) with 
!        size(A, 1) = m and size(A, 2) = n.
!        On entry, the matrix A.
!        On exit, if JOB = 'U' and U is not present, then A is 
!        overwritten with the first min(m, n) columns of U (the left
!        singular vectors, stored columnwise).
!        If JOB = 'V' and VT is not present, then A is overwritten with
!        the first min(m, n) rows of V^H (the right singular vectors, 
!        stored rowwise).
!        In all cases the original contents of A are destroyed.
! S      (output) REAL array, shape (:) with size(S) = min(m, n).
!        The singular values of A, sorted so that S(i) >= S(i+1).
! U      Optional (output) REAL or COMPLEX array, shape (:, :) with 
!        size(U, 1) = m  and size(U, 2) = m or min(m, n).
!        If size(U, 2) = m, U contains the m by m matrix U .
!        If size(U; 2) = min(m, n), U contains the first min(m, n) 
!        columns of U (the left singular vectors, stored columnwise).
! VT     Optional (output) REAL or COMPLEX array, shape (:, :) with 
!        size(VT, 1) = n or min(m, n) and size(VT, 2) = n.
!        If size(VT, 1) = n , VT contains the n by n matrix V^H .
!        If size(VT, 1) = min(m, n), VT contains the first min(m, n)
!        rows of V^H (the right singular vectors, stored rowwise).
! WW     Optional (output) REAL array, shape (:) with size(WW) = 
!        min(m, n) - 1
!        If INFO > 0, WW contains the unconverged superdiagonal elements
!        of an upper bidiagonal matrix B whose diagonal is in SIGMA (not
!        necessarily sorted). B has the same singular values as A.
!        Note: WW is a dummy argument for LA_GESDD.
! JOB    Optional (input) CHARACTER(LEN=1).
!        = 'N': neither columns of U nor rows of V^H are returned in 
!          array A.
!        = 'U': if U is not present, the first min(m, n) columns of U 
!          (the left singular vectors) are returned in array A;
!        = 'V': if VT is not present, the first min(m, n) rows of V^H 
!          (the right singular vectors) are returned in array A;
!        Default value: 'N'.
! INFO   Optional (output) INTEGER.
!        = 0: successful exit.
!        < 0: if INFO = -i, the i-th argument had an illegal value.
!        > 0: The algorithm did not converge.
!        If INFO is not present and an error occurs, then the program is
!        terminated with an error message.
!----------------------------------------------------------------------
!  .. LOCAL PARAMETERS ..
   CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_GESVD'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LJOB
   CHARACTER(LEN=1) :: LJOBU, LJOBVT
   INTEGER, SAVE :: LWORK = 0
   INTEGER :: N, M, LINFO, LD, ISTAT, ISTAT1, S1U, S2U, S1VT, S2VT, &
              NN, MN, SWW
!  .. LOCAL ARRAYS ..
   REAL(WP), TARGET :: LLU(1,1), LLVT(1,1)
   REAL(WP), POINTER :: WORK(:)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC MIN, MAX, PRESENT, SIZE
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; ISTAT = 0; M = SIZE(A,1); N = SIZE(A,2)
   LD = MAX(1,M); MN = MIN(M,N)
   IF( PRESENT(JOB) )THEN; LJOB = JOB; ELSE; LJOB = 'N'; ENDIF
   IF( PRESENT(U) )THEN; S1U = SIZE(U,1); S2U = SIZE(U,2)
   ELSE; S1U = 1; S2U = 1; END IF
   IF( PRESENT(VT) )THEN; S1VT = SIZE(VT,1); S2VT = SIZE(VT,2)
   ELSE; S1VT = 1; S2VT = 1; END IF
   IF( PRESENT(WW) )THEN; SWW = SIZE(WW); ELSE; SWW = MN-1; ENDIF
!  .. TEST THE ARGUMENTS
   IF( M < 0 .OR. N < 0 )THEN; LINFO = -1
   ELSE IF( SIZE( S ) /= MN )THEN; LINFO = -2
   ELSE IF( PRESENT(U) .AND. ( S1U /= M .OR. &
            ( S2U /= M .AND. S2U /= MN ) ) )THEN; LINFO = -3
   ELSE IF( PRESENT(VT) .AND. ( ( S1VT /= N .AND. S1VT /= MN ) &
            .OR. S2VT /= N ) )THEN; LINFO = -4
   ELSE IF( SWW /= MN-1 .AND. MN > 0 ) THEN; LINFO = -5
   ELSE IF( PRESENT(JOB) .AND. ( .NOT. ( LSAME(LJOB,'U') .OR. &
            LSAME(LJOB,'V') .OR. LSAME(LJOB,'N') ) .OR. &
            LSAME(LJOB,'U') .AND. PRESENT(U) .OR. &
            LSAME(LJOB,'V') .AND. PRESENT(VT)) )THEN; LINFO = -6
   ELSE
      IF( PRESENT(U) )THEN
         IF( S2U == M )THEN; LJOBU = 'A'; ELSE; LJOBU = 'S'; ENDIF
      ELSE; IF( LSAME(LJOB,'U') ) THEN; LJOBU = 'O'
         ELSE; LJOBU = 'N'; ENDIF
      ENDIF
      IF( PRESENT(VT) )THEN
         IF( S1VT == N )THEN; LJOBVT = 'A'; ELSE; LJOBVT = 'S'; ENDIF
      ELSE; IF( LSAME(LJOB,'V') )THEN; LJOBVT = 'O'
         ELSE; LJOBVT = 'N'; ENDIF
      ENDIF
      IF( ISTAT == 0 )THEN
         NN = MAX( 5, 3*MN + MAX(M,N), 5*MN )
	 LWORK = MAX( 1, NN, LWORK ); ALLOCATE(WORK(LWORK), STAT=ISTAT)
	 IF( ISTAT /= 0 )THEN; DEALLOCATE(WORK,STAT=ISTAT1)
	  LWORK = MAX( 1, NN); ALLOCATE(WORK(LWORK), STAT=ISTAT)
	  IF( ISTAT == 0) CALL ERINFO( -200, SRNAME, LINFO )
         END IF
       END IF
       IF( ISTAT == 0 ) THEN 
          IF( PRESENT(U) ) THEN
            IF ( PRESENT(VT) )THEN
              CALL DGESVD( LJOBU, LJOBVT, M, N, A, LD, S, U, MAX(1,S1U), &
             VT, MAX(1,S1VT), WORK, LWORK, LINFO )
           ELSE
             CALL DGESVD( LJOBU, LJOBVT, M, N, A, LD, S, U, MAX(1,S1U), &
 &             LLVT, MAX(1,S1VT), WORK, LWORK, LINFO )
           ENDIF
         ELSE
           IF ( PRESENT(VT) )THEN
             CALL DGESVD( LJOBU, LJOBVT, M, N, A, LD, S, LLU, MAX(1,S1U), &
 &             VT, MAX(1,S1VT), WORK, LWORK, LINFO )
           ELSE
             CALL DGESVD( LJOBU, LJOBVT, M, N, A, LD, S, LLU, MAX(1,S1U), &
 &             LLVT, MAX(1,S1VT), WORK, LWORK, LINFO )
           ENDIF
         ENDIF        
         LWORK = INT(WORK(1)+1)
         IF( LINFO > 0 .AND. PRESENT(WW) ) WW(1:MN-1) = WORK(2:MN)
      ELSE; LINFO = -100; ENDIF
      DEALLOCATE(WORK, STAT=ISTAT1)
   ENDIF
   CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE DGESVD_F95

SUBROUTINE DPOTRF_F95( A, UPLO, RCOND, NORM, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   !USE LA_PRECISION, ONLY: WP => DP
   !USE LA_AUXMOD, ONLY: ERINFO, LSAME
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. CHARACTER ARGUMENTS ..
   CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: NORM, UPLO
!  .. SCALAR ARGUMENTS ..
   INTEGER, INTENT(OUT), OPTIONAL  :: INFO
   REAL(WP), INTENT(OUT), OPTIONAL :: RCOND
!  .. ARRAY ARGUMENTS ..
   REAL(WP), INTENT(INOUT) :: A(:,:)
!-----------------------------------------------------------------
!
! Purpose
! =======
!
! LA_POTRF computes the Cholesky factorization of a real symmetric or
! complex Hermitian positive definite matrix A.
!
! The factorization has the form
!    A = U**H * U,  if UPLO = 'U', or
!    A = L * L**H,  if UPLO = 'L',
! where U is an upper triangular matrix and L is lower triangular.
!
! This is the block version of the algorithm, calling Level 3 BLAS.
!
! LA_POTRF optionally estimates the reciprocal of the condition number
! (in the 1-norm) of a real symmetric or complex Hermitian positive 
! definite matrix A.
! An estimate is obtained for norm(inv(A)), and the reciprocal of the
! condition number is computed as RCOND = 1 / (norm(A) * norm(inv(A))).
!
! =======
!
!    SUBROUTINE LA_POTRF( A, UPLO, RCOND, NORM, INFO )
!       <type>(<wp>), INTENT(INOUT) :: A(:,:)
!       CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
!       REAL(<wp>), INTENT(OUT), OPTIONAL :: RCOND
!       CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: NORM
!       INTEGER, INTENT(OUT), OPTIONAL :: INFO
!    where
!       <type> ::= REAL | COMPLEX
!       <wp>   ::= KIND(1.0) | KIND(1.0D0)
!
! Defaults
! ========
!
! 1. If UPLO is not present then UPLO = 'U' is assumed.
!
! Arguments
! =========
!
! A       (input/output) either REAL or COMPLEX square array, 
!         shape (:,:), size(A,1) == size(A,2) >= 0.
!         On entry, the symmetric (Hermitian) matrix A.  
!            If UPLO = 'U', the upper triangular part of A contains
!               the upper triangular part of the matrix A, and the 
!               strictly lower triangular part of A is not referenced.
!            If UPLO = 'L', the lower triangular part of A contains
!               the lower triangular part of the matrix A, and the
!               strictly upper triangular part of A is not referenced.
!         On exit, if INFO = 0, the factor U or L from the Cholesky
!            factorization A = U**H*U or A = L*L**H.
!
! UPLO    Optional, (input) CHARACTER*1
!         If UPLO is present then:
!            = 'U':  Upper triangle of A is stored;
!            = 'L':  Lower triangle of A is stored.
!         otherwise UPLO = 'U' is assumed.
!
! RCOND   Optional (output) REAL
!         The reciprocal of the condition number of the matrix A 
!         computed as RCOND = 1/(norm(A) * norm(inv(A))).
! NORM    Optional (input) CHARACTER*1
!         Specifies whether the 1-norm condition number or the
!         infinity-norm condition number is required:
!           If NORM is present then:
!              = '1', 'O' or 'o': 1-norm;
!              = 'I' or 'i': infinity-norm.
!           otherwise NORM = '1' is used.
!
! INFO    Optional, (output) INTEGER
!         If INFO is present:
!            = 0: successful exit
!            < 0: if INFO = -i, the i-th argument had an illegal value
!            > 0: if INFO = i, the leading minor of order i is not
!               positive definite, and the factorization could not be
!               completed.
!         If INFO is not present and an error occurs, then the program
!            is terminated with an error message.
!
! --------------------------------------
!  .. LOCAL PARAMETERS ..
   CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_POTRF'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LNORM, LUPLO
   INTEGER :: LINFO, N, ISTAT, ISTAT1, LD
   REAL(WP) :: ANORM
!  .. LOCAL POINTERS ..
   INTEGER, POINTER :: IWORK(:)
   REAL(WP), POINTER :: WORK(:)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC PRESENT, MAX
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; N = SIZE(A,1); LD = MAX(1,N); ISTAT = 0
   IF( PRESENT(UPLO) ) THEN; LUPLO = UPLO; ELSE; LUPLO = 'U'; END IF
   IF( PRESENT(NORM) ) THEN; LNORM = NORM; ELSE; LNORM = '1'; END IF
!  .. TEST THE ARGUMENTS
   IF( SIZE( A, 2 ) /= N .AND. N < 0 )THEN; LINFO = -1
   ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN; LINFO = -2
   ELSE IF( ( .NOT.PRESENT(RCOND) .AND. PRESENT(NORM) ) .OR. &
            ( .NOT.LSAME(LNORM,'I') .AND. .NOT.LSAME(LNORM,'O') &
              .AND. LNORM /= '1' ) ) THEN; LINFO = -4
   ELSE IF(  N > 0 )THEN
      IF( PRESENT(RCOND) ) THEN
!     .. COMPUTE THE NORM OF THE MATRIX A
         ALLOCATE(WORK(N), STAT=ISTAT)
         IF( ISTAT == 0 )THEN; ANORM = DLANSY( LNORM, LUPLO, LD, A, N, WORK )
         ELSE; LINFO = -100; END IF
         DEALLOCATE(WORK, STAT=ISTAT1)
      END IF
!
      IF( LINFO == 0 ) THEN
!     .. COMPUTE THE CHOLESKY FACTORS OF THE MATRIX A
         CALL DPOTRF( LUPLO, N, A, LD, LINFO )
!
         IF( PRESENT(RCOND) .AND. LINFO == 0 ) THEN
!        .. COMPUTE THE RECIPROCAL OF THE CONDITION NUMBER OF A
            IF( ANORM == 0.0_WP )THEN; RCOND = 0.0_WP
            ELSE; ALLOCATE(WORK(3*N), IWORK(N), STAT=ISTAT)
               IF( ISTAT == 0 )THEN
                  CALL DPOCON( LUPLO, N, A, LD, ANORM, RCOND, &
                                  WORK, IWORK, LINFO )
               ELSE; LINFO = -100; END IF
               DEALLOCATE(WORK, IWORK, STAT=ISTAT1)
            END IF
         END IF
      END IF
   ELSE IF( PRESENT(RCOND) ) THEN; RCOND = 1.0_WP; ENDIF
   CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE DPOTRF_F95


!------------------------------------------
      LOGICAL FUNCTION LSAME( CA, CB )
!
!  PURPOSE
!  =======
!
!  LSAME  TESTS IF CA IS THE SAME LETTER AS CB REGARDLESS OF CASE.
!
!  PARAMETERS
!  ==========
!
!  CA      (INPUT) CHARACTER*1
!  CB      (INPUT) CHARACTER*1
!          CHARACTERS TO BE COMPARED.
!
!  .. SCALAR ARGUMENTS ..
      CHARACTER*1, INTENT(IN) :: CA, CB
!  .. PARAMETERS ..
      INTEGER, PARAMETER      :: IOFF=32
!  .. LOCAL SCALARS ..
      INTEGER                 :: INTA, INTB, ZCODE
!  .. INTRINSIC FUNCTIONS ..
      INTRINSIC                  ICHAR
!
!  .. EXECUTABLE STATEMENTS ..
!
!  TEST IF THE CHARACTERS ARE EQUAL
!
      LSAME = CA == CB
!
!  NOW TEST FOR EQUIVALENCE
!
      IF( .NOT.LSAME )THEN
!
!     USE 'Z' RATHER THAN 'A' SO THAT ASCII CAN BE DETECTED ON PRIME
!     MACHINES, ON WHICH ICHAR RETURNS A VALUE WITH BIT 8 SET.
!     ICHAR('A') ON PRIME MACHINES RETURNS 193 WHICH IS THE SAME AS
!     ICHAR('A') ON AN EBCDIC MACHINE.
!
         ZCODE = ICHAR( 'Z' )
!
         INTA = ICHAR( CA )
         INTB = ICHAR( CB )
!
         IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 )THEN
!
!        ASCII IS ASSUMED - ZCODE IS THE ASCII CODE OF EITHER LOWER OR
!        UPPER CASE 'Z'.
!
            IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
            IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
!
         ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 )THEN
!
!        EBCDIC IS ASSUMED - ZCODE IS THE EBCDIC CODE OF EITHER LOWER OR
!        UPPER CASE 'Z'.
!
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.                         &
!    &       INTA.GE.145 .AND. INTA.LE.153 .OR.                         &
     &       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.                         &
     &       INTB.GE.145 .AND. INTB.LE.153 .OR.                         &
     &       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
!
         ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 )THEN
!
!        ASCII IS ASSUMED, ON PRIME MACHINES - ZCODE IS THE ASCII CODE
!        PLUS 128 OF EITHER LOWER OR UPPER CASE 'Z'.
!
            IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
         ENDIF
         LSAME = INTA == INTB
      ENDIF
      END FUNCTION LSAME

!------------------------------------------------
      SUBROUTINE ERINFO(LINFO, SRNAME, INFO, ISTAT)
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. IMPLICIT STATEMENT ..
         IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
         CHARACTER( LEN = * ), INTENT(IN)              :: SRNAME
         INTEGER             , INTENT(IN)              :: LINFO
         INTEGER             , INTENT(OUT), OPTIONAL   :: INFO
         INTEGER             , INTENT(IN), OPTIONAL    :: ISTAT
!  .. EXECUTABLE STATEMENTS ..
!         IF( ( LINFO < 0 .AND. LINFO > -200 ) .OR.                     &
!    &       ( LINFO > 0 .AND. .NOT.PRESENT(INFO) ) )THEN
      IF( ( ( LINFO < 0 .AND. LINFO > -200 ) .OR. LINFO > 0 )           &
     &           .AND. .NOT.PRESENT(INFO) )THEN
        WRITE (*,*) 'Program terminated in LAPACK95 subroutine ',SRNAME
        WRITE (*,*) 'Error indicator, INFO = ',LINFO
        IF( PRESENT(ISTAT) )THEN
          IF( ISTAT /= 0 ) THEN
            IF( LINFO == -100 )THEN
              WRITE (*,*) 'The statement ALLOCATE causes STATUS = ',    &
     &                    ISTAT
            ELSE
              WRITE (*,*) 'LINFO = ', LINFO, ' not expected'
            END IF
          END IF   
        END IF
        STOP
         ELSE IF( LINFO <= -200 ) THEN
           WRITE(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++'
           WRITE(*,*) '*** WARNING, INFO = ', LINFO, ' WARNING ***'
           IF( LINFO == -200 )THEN
             WRITE(*,*)                                                 &
     &        'Could not allocate sufficient workspace for the optimum'
             WRITE(*,*)                                                 &
     &        'blocksize, hence the routine may not have performed as'
             WRITE(*,*) 'efficiently as possible'
         ELSE
           WRITE(*,*) 'Unexpected warning'
         END IF
           WRITE(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++'
        END IF
        IF( PRESENT(INFO) ) THEN
          INFO = LINFO
        END IF
      END SUBROUTINE ERINFO

end module
