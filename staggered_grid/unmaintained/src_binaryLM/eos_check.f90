module eos_check_module

  use multifab_module
  use ml_layout_module
  use probin_common_module, only: rhobar

  implicit none

  private

  public :: eos_check

contains

  subroutine eos_check(mla,s)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: s(:)

    ! local
    integer i,n,dm,nlevs,ng_s
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: sp(:,:,:,:)

    real(kind=dp_t) :: eos_error, eos_error_grid, eos_error_proc

    nlevs = mla%nlevel
    dm = mla%dim

    ng_s = s(1)%ng
    eos_error_proc = -1.d20
    do n=1,nlevs
       do i=1,nfabs(s(n))
         sp => dataptr(s(n), i)
         lo =  lwb(get_box(s(n), i))
         hi =  upb(get_box(s(n), i))
         eos_error_grid = -1.d20
         select case (dm)
         case (2)
            call eos_check_2d(sp(:,:,1,:),ng_s,eos_error_grid,lo,hi)
         case (3)
            call eos_check_3d(sp(:,:,:,:),ng_s,eos_error_grid,lo,hi)
         end select
         eos_error_proc = max(eos_error_grid, eos_error_proc)
      end do
   end do

   ! This sets eos_error to be the max of eos_error_proc over all processors.
   call parallel_reduce(eos_error, eos_error_proc, MPI_MAX)

   if (parallel_IOProcessor()) then
      print*,"EOS ERROR in L1 norm: ",eos_error
      print*,""
   end if

  end subroutine eos_check

  subroutine eos_check_2d(s,ng_s,eos_error,lo,hi)

    integer        , intent(in   ) :: lo(:), hi(:), ng_s
    real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(inout) :: eos_error

    ! local
    integer :: i,j

    real(kind=dp_t) :: error

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          error = abs( (s(i,j,2)/rhobar(1) + (s(i,j,1)-s(i,j,2))/rhobar(2) ) - 1.d0)
          eos_error = max(eos_error,error)

       end do
    end do

  end subroutine eos_check_2d

  subroutine eos_check_3d(s,ng_s,eos_error,lo,hi)

    integer        , intent(in   ) :: lo(:), hi(:), ng_s
    real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(inout) :: eos_error

    ! local
    integer :: i,j,k

    real(kind=dp_t) :: error

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             error = abs( (s(i,j,k,2)/rhobar(1) + (s(i,j,k,1)-s(i,j,k,2))/rhobar(2) ) - 1.d0)
             eos_error = max(eos_error,error)
             
          end do
       end do
    end do

  end subroutine eos_check_3d

end module eos_check_module
