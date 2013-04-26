module exact_solutions_module

  use bl_types
       
  implicit none

  private

  public :: exact_2d, exact_3d

contains

  subroutine exact_2d(mx,my,s,lo,hi,ng_m,ng_s,dx,time)

    integer        , intent(in   ) :: lo(:), hi(:), ng_m, ng_s
    real(kind=dp_t), intent(  out) :: mx(lo(1)-ng_m:,lo(2)-ng_m:)
    real(kind=dp_t), intent(  out) :: my(lo(1)-ng_m:,lo(2)-ng_m:)
    real(kind=dp_t), intent(  out) ::  s(lo(1)-ng_s:,lo(2)-ng_s:,:)  
    real(kind=dp_t), intent(in   ) :: dx(:),time

    ! local
    integer :: i,j

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          s(i,j,1) = 1.d0
          s(i,j,2) = 1.d0
       enddo
    enddo

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1
          mx(i,j) = 0.d0
       enddo
    enddo

    do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          my(i,j) = 0.d0
       enddo
    enddo

  end subroutine exact_2d

  subroutine exact_3d(mx,my,mz,s,lo,hi,ng_m,ng_s,dx,time)

    integer        , intent(in   ) :: lo(:), hi(:), ng_m, ng_s
    real(kind=dp_t), intent(  out) :: mx(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real(kind=dp_t), intent(  out) :: my(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real(kind=dp_t), intent(  out) :: mz(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real(kind=dp_t), intent(  out) ::  s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)  
    real(kind=dp_t), intent(in   ) :: dx(:),time

    ! local
    integer :: i,j,k

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             s(i,j,k,1) = 1.d0
             s(i,j,k,2) = 1.d0
          enddo
       enddo
    end do

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
             mx(i,j,k) = 0.d0
          enddo
       enddo
    enddo

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
             my(i,j,k) = 0.d0
          enddo
       enddo
    enddo

    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             mz(i,j,k) = 0.d0
          enddo
       enddo
    enddo

  end subroutine exact_3d

end module exact_solutions_module
