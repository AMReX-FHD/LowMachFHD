module matvec_mul_module

  use multifab_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use ml_layout_module

  implicit none

  private

  public :: matvec_mul

  contains

  subroutine matvec_mul(mla, x, A) ! Computes x=A*x, where A is an nxn matrix and x is an nx1 vector  
   
    type(ml_layout), intent(in  )  :: mla
    type(multifab),  intent(inout) :: x
    type(multifab),  intent(in)    :: A
 
    ! local variables
    real(kind=dp_t), pointer       :: xp(:,:,:,:) 
    real(kind=dp_t), pointer       :: ap(:,:,:,:)
    integer                        :: lo(mla%dim), hi(mla%dim), i, dm, nc 
   
    dm = mla%dim 
    nc = x%nc    
  
    do i=1,nfabs(x)
       xp  => dataptr(x, i)   
       ap  => dataptr(A, i)  
       lo(1) = lbound(ap,1)  ! this adjusts lo & hi ghost cells on nodes & faces,
       lo(2) = lbound(ap,2)  ! so in the subroutine loop from lo to hi without ng.
       hi(1) = ubound(ap,1)
       hi(2) = ubound(ap,2)
       !print*, 'lo(1), hi(1)=',lo(1),hi(1), 'lo(2), hi(2)=',lo(2),hi(2)

       select case (dm)
         case (2)
             call matvec_mul_2d(xp(:,:,1,:), ap(:,:,1,:), lo, hi, nc)
         case (3)
             stop "3d matvec_mul not yet implemented"
             lo(3) = lbound(ap,3) 
             hi(3) = ubound(ap,3)
             call matvec_mul_3d(xp(:,:,:,:), ap(:,:,:,:), lo, hi, nc)
       end select
    end do

  end subroutine matvec_mul

  subroutine matvec_mul_2d(xp, ap, lo, hi, nc)

    integer                        :: lo(:), hi(:) 
    real(kind=dp_t), intent(inout) :: xp(lo(1):,lo(2):,:) ! last dimension for nc
    real(kind=dp_t), intent(in)    :: ap(lo(1):,lo(2):,:)
    integer                        :: i,j,nc

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          call matvec_mul_comp(xp(i,j,:), ap(i,j,:))
       end do
    end do

    contains 
    
    ! Use contained subroutine to do rank conversion and mat-vec mult
    subroutine matvec_mul_comp(xp_ij, ap_ij)

        real(kind=dp_t), dimension(nc),    intent(inout) :: xp_ij
        real(kind=dp_t), dimension(nc,nc), intent(in)    :: ap_ij  
        
        ! local variables
        ! use dummy matrix cp_ij to store the matrix-vector multiplication
        real(kind=dp_t), dimension(nc)  :: cp_ij  
        cp_ij=0.d0
        
        cp_ij = matmul(ap_ij, xp_ij)
        ! Donev: In Fortran you can write directly xp_ij = matmul(ap_ij, xp_ij)
        ! This is better as it may be optimized by the compiler better
        xp_ij = cp_ij
 
    end subroutine matvec_mul_comp 

  end subroutine matvec_mul_2d
  
end module matvec_mul_module
