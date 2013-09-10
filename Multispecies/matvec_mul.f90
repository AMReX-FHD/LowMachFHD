module matvec_mul_module

  use multifab_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use ml_layout_module
  use probin_multispecies_module 

  implicit none

  private

  public :: matvec_mul

  contains

subroutine matvec_mul(x, A) ! Computes x=A*x, where A is an nxn matrix and x is an nx1 vector    
    type(multifab), intent(inout) :: x
    type(multifab), intent(in)    :: A

    integer                       :: lo(2), hi(2), dm, ng, nc 
    
    ! local
    real(kind=dp_t), pointer      :: ap(:,:,:,:)
    real(kind=dp_t), pointer      :: xp(:,:,:,:) 
    
    do i=1,nfabs(x(n))

       ! Donev: Unfinished
       ap  => dataptr(a(n), i)
       xp  => dataptr(x(n), i)
       lo(1) = lbound(pp,1)
       lo(2) = lbound(pp,2)
       hi(1) = ubound(pp,1)
       hi(2) = ubound(pp,2)

       ! Do loops should go from lo:hi 


          select case (dm)
          case (2)
             call matvec_mul_2d(ap, xp, nc, lo, hi)
          case (3)
             stop "3d matvec_mul not implemented"
          end select
    end do

end subroutine matvec_mul

  end subroutine matvec_mul

  subroutine matvec_mul_2d(ap, bp, lo, hi, ng_f, dim_in)

    real(kind=dp_t), pointer       :: ap(:,:,:,:)
    real(kind=dp_t), pointer       :: bp(:,:,:,:)
    !real(kind=dp_t), intent(inout) :: ap(lo(1)-ng_f:,lo(2)-ng_f:,:)
    !real(kind=dp_t), intent(in)    :: bp(lo(1)-ng_f:,lo(2)-ng_f:,:)
    integer :: lo(2), hi(2), ng_f, i, j, k, m, n, dim_in

    do k = lbound(ap,dim=3), ubound(ap,dim=3)       ! 1:Lz
       !do j = lbound(ap,dim=2), ubound(ap,dim=2)    ! 1:Ly
       !   do i = lbound(ap,dim=1), ubound(ap,dim=1) ! 1:Lx

    select case(dim_in)
    case(1) 
    do j=lo(2)-ng_f,hi(2)+ng_f
       do i=lo(1)-ng_f,hi(1)+ng_f+1
          !call matvec_mul_comp(ap(i,j,:), bp(i,j,:), i, j)
          call matvec_mul_comp(ap(i,j,k,:), bp(i,j,k,:), i, j)
       end do
    end do
    case(2)
     do j=lo(2)-ng_f,hi(2)+ng_f+1
       do i=lo(1)-ng_f,hi(1)+ng_f
          !call matvec_mul_comp(ap(i,j,:), bp(i,j,:), i, j)
          call matvec_mul_comp(ap(i,j,k,:), bp(i,j,k,:), i, j)
       end do
    end do
    end select 

    end do

    contains 
    
    ! Use contained (internal) subroutine to do rank conversion and
    ! matrix-vector multiplication 
    subroutine matvec_mul_comp(ap_ij, bp_ij, i, j)

        real(kind=dp_t), dimension(nspecies),       intent(inout) :: ap_ij
        real(kind=dp_t), dimension(nspecies,nspecies), intent(in) :: bp_ij  
        integer,                                       intent(in) :: i,j 
        
        ! local variables
        ! use dummy matrix cp_ij to store the matrix-vector multiplication
        real(kind=dp_t), dimension(nspecies)  :: cp_ij  
        real(kind=dp_t)                       :: mvprod

        do n=1, nspecies 
           mvprod=0.d0
           do m=1, nspecies
              mvprod = mvprod + bp_ij(n,m)*ap_ij(m)
              if(i.eq.9 .and. j.eq.11) then
                 !print*, i, j, bp_ij(n,m)
              endif
           enddo
           cp_ij(n) = mvprod
        enddo      
        ! populate ap_ij with cp_ij 
        ap_ij = cp_ij

     end subroutine matvec_mul_comp 

  end subroutine matvec_mul_2d
  
end module matvec_mul_module
