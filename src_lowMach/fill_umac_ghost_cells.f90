module fill_umac_ghost_cells_module

  use ml_layout_module
  use multifab_module
  use define_bc_module
  use bc_module
  use reservoir_bc_fill_module
  use multifab_physbc_stag_module

  implicit none

  private

  public :: fill_umac_ghost_cells

contains

  subroutine fill_umac_ghost_cells(mla,umac,eta_ed,dx,time,the_bc_tower)

    ! fill the ghost cells for the mac velocity (used for restarts to restore ghost values)
    ! does not modify the domain boundary values, i.e., it does not call multifab_physbc_domainvel
    ! it assumes that normal velocity on physical domain boundaries are properly set (e.g., read from restart)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(inout) :: eta_ed(:,:) ! nodal (2d); edge-centered (3d)
    real(kind=dp_t), intent(in   ) :: dx(:,:),time
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: n,nlevs,i,dm

    type(bl_prof_timer),save :: bpt

    call build(bpt,"fill_umac_ghost_cells")

    nlevs = mla%nlevel
    dm = mla%dim
    ! vel_bc_n here is never used: the normal velocities are assumed to be set already in umac
    call set_inhomogeneous_vel_bcs(mla,vel_bc_n,vel_bc_t,eta_ed,dx,time, &
                                   the_bc_tower%bc_tower_array)

    do n=1,nlevs
       do i=1,dm
          ! set transverse velocity behind physical boundaries
          call multifab_physbc_macvel(umac(n,i),vel_bc_comp+i-1, &
                                      the_bc_tower%bc_tower_array(n), &
                                      dx(n,:),vel_bc_t(n,:))
          ! fill periodic and interior ghost cells
          call multifab_fill_boundary(umac(n,i))
       end do
    end do

    call destroy(bpt)

  end subroutine fill_umac_ghost_cells

end module fill_umac_ghost_cells_module
