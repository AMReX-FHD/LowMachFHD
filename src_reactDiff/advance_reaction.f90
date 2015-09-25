module advance_reaction_module

  use ml_layout_module
  use bc_module
  use define_bc_module
  use multifab_physbc_module
  use probin_common_module, only: seed
  use probin_reactdiff_module, only: nspecies, reaction_type

  implicit none

  private

  public :: advance_reaction

contains

  subroutine advance_reaction(mla,n_old,n_new,dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: n_old(:)
    type(multifab) , intent(inout) :: n_new(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: n,nlevs

    real(kind=dp_t) :: tau
    real(kind=dp_t) :: n_react

    type(gsl_rng) :: test

    print*,gsl_ran_poisson(test,1.d0)

    nlevs = mla%nlevel

    if (reaction_type .eq. 0 .or. reaction_type .eq. 1) then
       ! first-order tau-leaping

       print*,n_react

       if (reaction_type .eq. 1) then
          ! second-order tau-leaping corrector

       end if

    else if (reaction_type .eq. 2 .or. reaction_type .eq. 3) then
       ! first-order CLE

       if (reaction_type .eq. 3) then
          ! second-order CLE corrector

       end if

    else if (reaction_type .eq. 4) then
       ! SSA

    else
       call bl_error("advance_reaction: invalid reaction_type")
    end if

    do n=1,nlevs
       call multifab_fill_boundary(n_new(n))
       call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                            the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
    end do

  end subroutine advance_reaction

end module advance_reaction_module
