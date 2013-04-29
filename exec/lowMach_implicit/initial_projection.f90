module initial_projection_module

  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: initial_projection

contains

  subroutine initial_projection(mla,mold,umac,sold)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: mold(:,:)
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(in   ) :: sold(:)

    ! local
    type(multifab) :: rho_edge(mla%nlevel,mla%dim)

    ! create rho on faces

    ! convert m to u

    ! build rhs = S - div(u)

    ! project to solve for phi - use the 'full' solver

    ! subtract off pressure gradient

    ! convert u to m


  end subroutine initial_projection

end module initial_projection_module
