module debug_module

  use bl_types
  use multifab_module

  implicit none

  private

  public :: print_edge

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_edge(mf,face,ival,jval,comp)

    type(multifab) , intent(in   ) :: mf(:,:)
    integer        , intent(in   ) :: face,ival,jval,comp

    real(kind=dp_t), pointer:: dp(:,:,:,:)

    integer :: lo(mf(1,face)%dim),hi(mf(1,face)%dim)
    integer :: n,i,dm,nlevs,ng

    dm = mf(1,face)%dim
    nlevs = 1
    ng = mf(1,face)%ng
    do n = 1, nlevs
       do i = 1, nfabs(mf(n,face))
          dp => dataptr(mf(n,face), i)
          lo = lwb(get_box(mf(n,face), i))
          hi = upb(get_box(mf(n,face), i))
          call print_edge_2d(dp(:,:,1,:), ng, face, ival, jval, comp, lo, hi)
       end do
    enddo

  end subroutine print_edge
  
  subroutine print_edge_2d(mf, ng, face, ival, jval, comp, lo, hi)
    
    integer          :: ng, face, ival, jval, comp, lo(2), hi(2)
    real (kind=dp_t) :: mf(lo(1)-ng:,lo(2)-ng:,:)

    ! local
    integer :: i,j

    integer :: ioff, joff

    ioff = 0
    joff = 0
    
    if (face .eq. 1) then
       ioff = 1
    end if
    
    if (face .eq. 2) then
       joff = 1
    end if
    
    do j=lo(2),hi(2)+joff
    do i=lo(1),hi(1)+ioff

       if (i .eq. ival .and. j .eq. jval) then
          print*,'VALUE',i,j,comp,mf(i,j,comp)
       end if
             
    end do
    end do
    
  end subroutine print_edge_2d

end module debug_module
