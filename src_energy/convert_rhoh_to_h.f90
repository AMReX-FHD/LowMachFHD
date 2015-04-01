module convert_rhoh_to_h_module

  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: convert_rhoh_to_h
  
contains

  subroutine convert_rhoh_to_h(mla,rhoh,rhotot,h,rhoh_to_h)
    
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rhoh(:)
    type(multifab) , intent(in   ) :: rhotot(:)
    type(multifab) , intent(inout) :: h(:)
    logical        , intent(in   ) :: rhoh_to_h

    ! local
    integer :: n,nlevs,i

    nlevs = mla%nlevel

    if (rhoh_to_h) then

       ! rhoh to h - NO GHOST CELLS
       do n=1,nlevs
          call multifab_copy_c(h(n),1,rhoh(n),1,1,0)
          call multifab_div_div_c(h(n),1,rhotot(n),1,1,0)
       end do

    else

       ! h to rhoh- VALID + GHOST (CAN CHANGE TO DO ONLY GHOST TO SAVE COMPUTATION)
       do n=1,nlevs
          call multifab_copy_c(rhoh(n),1,h(n),1,1,rhoh(n)%ng)
          call multifab_mult_mult_c(rhoh(n),1,rhotot(n),1,1,rhoh(n)%ng)
       end do

    end if

  end subroutine convert_rhoh_to_h

end module convert_rhoh_to_h_module
