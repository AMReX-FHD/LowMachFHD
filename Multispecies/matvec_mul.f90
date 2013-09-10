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

subroutine matvec_mul(a, targ, b, src, nc, ng, lo, hi, dm, ng_f)
    
    integer, intent(in)           :: targ, src
    integer, intent(in)           :: nc
    integer, intent(in), optional :: ng
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)    :: b
    integer                       :: lo(2), hi(2), dm, ng_f 
    
    ! local
    real(kind=dp_t), pointer      :: ap(:,:,:,:)
    real(kind=dp_t), pointer      :: bp(:,:,:,:)  ! the last entry is nspecies^2.
    integer                       :: i,lng
    
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= ng,"not enough ghost cells in matvec_mul")
    do i = 1, nlocal(a%la)
       if ( lng > 0 ) then
          ap => dataptr(a, i, grow(get_ibox(a, i),lng), targ, nc)
          bp => dataptr(b, i, grow(get_ibox(b, i),lng), src, nc**2)
       else
          ap => dataptr(a, i, get_ibox(a, i), targ, nc)
          bp => dataptr(b, i, get_ibox(b, i), src, nc**2)
!          ap => dataptr(a(1),1)
!          bp => dataptr(b(1),1)
       end if

       select case (dm)     
       case (2)
          !call matvec_mul_2d(ap(:,:,1,:), bp(:,:,1,:), lo, hi, ng_f)
          call matvec_mul_2d(ap, bp, lo, hi, ng_f)
       end select
    end do
end subroutine matvec_mul

 subroutine matvec_mul_2d(ap, bp, lo, hi, ng_f)

    real(kind=dp_t), pointer       :: ap(:,:,:,:)
    real(kind=dp_t), pointer       :: bp(:,:,:,:)
    !real(kind=dp_t), intent(inout) :: ap(lo(1)-ng_f:,lo(2)-ng_f:,:)
    !real(kind=dp_t), intent(in)    :: bp(lo(1)-ng_f:,lo(2)-ng_f:,:)
    integer :: lo(2), hi(2), ng_f, i, j, k, m, n

    do k = lbound(ap,dim=3), ubound(ap,dim=3)       ! 1:Lz
       do j = lbound(ap,dim=2), ubound(ap,dim=2)    ! 1:Ly
          do i = lbound(ap,dim=1), ubound(ap,dim=1) ! 1:Lx
     
    !do j=lo(2)-ng_f,hi(2)+ng_f+1
    !   do i=lo(1)-ng_f,hi(1)+ng_f+1
          !call matvec_mul_comp(ap(i,j,:), bp(i,j,:), i, j)
          call matvec_mul_comp(ap(i,j,k,:), bp(i,j,k,:), i, j)
       end do
    end do
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
