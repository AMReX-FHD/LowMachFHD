module advance_timestep_module

  use ml_layout_module
  use multifab_module

  implicit none

  private

  public :: advance_timestep

contains

  subroutine advance_timestep(mla,mold,mnew,umac,sold,snew,eta,chi)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: mold(:,:)
    type(multifab) , intent(inout) :: mnew(:,:)
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(in   ) :: eta(:)
    type(multifab) , intent(in   ) :: chi(:)

    ! local

  end subroutine advance_timestep

end module advance_timestep_module
