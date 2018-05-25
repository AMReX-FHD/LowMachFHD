module utility_module

  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: invert_multifab

contains

  ! currently invert a cell-centered multifab
  ! only works for cell-centered at the moment
  subroutine invert_multifab(la,phi,comp,ncomp,nghost,numerator)

    type(layout)   , intent(in   ) :: la
    type(multifab) , intent(inout) :: phi
    integer        , intent(in   ) :: comp, ncomp, nghost
    real(kind=dp_t), intent(in), optional :: numerator

    ! local
    type(bl_prof_timer), save :: bpt

    type(multifab) :: temp

    integer :: i
    real(kind=dp_t) :: num

    call build(bpt,"invert_multifab")
    
    if(present(numerator)) then
      num=numerator
    else
      num=1.0d0
    end if    
    
    call multifab_build(temp,la,1,nghost) ! Temporary multifab
    do i=comp,comp+ncomp-1
       call setval(temp,numerator,all=.true.)
       call multifab_div_div_c(temp,1,phi,i,1,nghost)
       call multifab_copy_c(phi,i,temp,1,1,nghost)
    end do

    call multifab_destroy(temp)
    
    call destroy(bpt)

  end subroutine invert_multifab

end module utility_module
