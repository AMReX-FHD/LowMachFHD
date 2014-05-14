module convert_variables_module

  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: convert_cons_to_prim
  
contains

  subroutine convert_cons_to_prim(mla,s,prim,cons_to_prim)
    
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) ::    s(:)
    type(multifab) , intent(inout) :: prim(:)
    logical        , intent(in   ) :: cons_to_prim

    ! local
    integer :: nlevs,n,i

    nlevs = mla%nlevel

    if (cons_to_prim) then

       ! cons to prim - NO GHOST CELLS
       do n=1,nlevs
          call multifab_copy_c(prim(n),1,s(n),1,2,0)
          do i=2,2
             call multifab_div_div_c(prim(n),i,s(n),1,1,0)
          end do
       end do

    else

       ! prim to cons - INCLUDING GHOST CELLS
       do n=1,nlevs
          call multifab_copy_c(s(n),1,prim(n),1,2,s(n)%ng)
          do i=2,2
             call multifab_mult_mult_c(s(n),i,s(n),1,1,s(n)%ng)
          end do
       end do

    end if
    
  end subroutine convert_cons_to_prim

end module convert_variables_module
