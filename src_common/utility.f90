module utility_module

  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: invert_multifab

contains

  ! currently invert a cell-centered multifab
  ! only works for cell-centered at the moment
  subroutine invert_multifab(mla,phi,comp,ncomp,nghost)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: phi(:)
    integer        , intent(in   ) :: comp, ncomp, nghost

    ! local
    type(bl_prof_timer), save :: bpt

    type(multifab) :: temp(mla%nlevel)

    integer :: n, nlevs
    integer :: i

    call build(bpt,"invert_multifab")
    
    nlevs = mla%nlevel

    do n=1,nlevs
       call multifab_build(temp(n),mla%la(n),1,nghost)
    end do

    do n=1,nlevs
       do i=comp,comp+ncomp-1
          call setval(temp(n),1.d0,all=.true.)
          call multifab_div_div_c(temp(n),1,phi(n),i,1,nghost)
          call multifab_copy_c(phi(n),i,temp(n),1,1,nghost)
       end do
    end do

    do n=1,nlevs
       call multifab_destroy(temp(n))
    end do
    
    call destroy(bpt)

  end subroutine invert_multifab

end module utility_module
