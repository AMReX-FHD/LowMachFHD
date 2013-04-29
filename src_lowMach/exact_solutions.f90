module exact_solutions_module

  use bl_types
  use bl_error_module
  use probin_common_module , only: prob_lo, prob_hi
  use probin_lowmach_module, only: prob_type, rhobar
       
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
    real(kind=dp_t) :: x,y,r
    real(kind=dp_t) :: c_init(2),smoothing_width

    select case (prob_type)
    case (1)

       ! lo density spherical bubble

       mx = 0.d0
       my = 0.d0

       c_init(1) = 0.1d0
       c_init(2) = 0.9d0

       smoothing_width = 5.d0

       do j=lo(2),hi(2)
          y = prob_lo(2) + dx(2) * (dble(j)+0.5d0) - 0.5d0*(prob_lo(2)+prob_hi(2))
          do i=lo(1),hi(1)
             x = prob_lo(1) + dx(1) * (dble(i)+0.5d0) - 0.5d0*(prob_lo(1)+prob_hi(1))

             r = sqrt (x**2 + y**2)

             ! tanh smoothing
             s(i,j,2) = c_init(1) + 0.5d0*(c_init(2)-c_init(1))* &
                  (1.d0 + tanh((r-2.5d0*smoothing_width*dx(1))/(smoothing_width*dx(1))))

             s(i,j,1) = 1.0d0/(s(i,j,2)/rhobar(1)+(1.0d0-s(i,j,2))/rhobar(2))
             s(i,j,2) = s(i,j,1)*s(i,j,2)
          enddo
       enddo

    case default

       call bl_error("exact_2d: invalid prob_type")

    end select

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

    select case (prob_type)
    case (1)
       
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                s(i,j,k,1) = 1.d0
                s(i,j,k,2) = 1.d0
             enddo
          enddo
       end do
       
       mx = 0.d0
       my = 0.d0
       mz = 0.d0

    case default

       call bl_error("exact_3d: invalid prob_type")
          
    end select

  end subroutine exact_3d

end module exact_solutions_module
