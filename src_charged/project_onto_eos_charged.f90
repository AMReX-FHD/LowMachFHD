module project_onto_eos_charged_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use probin_common_module, only: rhobar
  use probin_multispecies_module, only: nspecies
  use probin_charged_module, only: charge_per_mass

  implicit none

  private

  public :: project_onto_eos_charged

contains

  subroutine project_onto_eos_charged(mla,rho)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rho(:)

    ! local
    integer i,n,nlevs,dm,ng_r
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: sp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"project_onto_eos_charged")

    if (charge_per_mass(nspecies) .ne. 0.d0) then
       call bl_error("project_onto_eos_charged assumes last species is neutral")
    end if

    nlevs = mla%nlevel
    dm    = mla%dim

    ng_r = rho(1)%ng

    do n=1,nlevs
       do i=1,nfabs(rho(n))
          sp  => dataptr(rho(n), i)
          lo = lwb(get_box(rho(n), i))
          hi = upb(get_box(rho(n), i))
          select case (dm)
          case (2)
             call project_onto_eos_charged_2d(sp(:,:,1,:), ng_r, lo, hi)
          case (3)
             call project_onto_eos_charged_3d(sp(:,:,:,:), ng_r, lo, hi)
          end select
       end do
    end do

    call destroy(bpt)

  end subroutine project_onto_eos_charged

  subroutine project_onto_eos_charged_2d(rho,ng_r,lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:),ng_r
    real(kind=dp_t), intent(inout) :: rho(lo(1)-ng_r:,lo(2)-ng_r:,:)

    ! integer
    integer :: i,j,comp

    real(kind=dp_t) :: rhosum

    do j=lo(2),hi(2)
    do i=lo(1),hi(1)

       rhosum = 1.d0
       do comp=1,nspecies-1
          rhosum = rhosum - rho(i,j,comp)/rhobar(comp)
       end do

       rho(i,j,nspecies) = rhosum*rhobar(nspecies)

    end do
    end do

  end subroutine project_onto_eos_charged_2d

  subroutine project_onto_eos_charged_3d(rho,ng_r,lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:),ng_r
    real(kind=dp_t), intent(inout) :: rho(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:,:)

    ! integer
    integer :: i,j,k,comp

    real(kind=dp_t) :: rhosum

    do k=lo(3),hi(3)
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)

       rhosum = 1.d0
       do comp=1,nspecies-1
          rhosum = rhosum - rho(i,j,k,comp)/rhobar(comp)
       end do

       rho(i,j,k,nspecies) = rhosum*rhobar(nspecies)

    end do
    end do
    end do

  end subroutine project_onto_eos_charged_3d

end module project_onto_eos_charged_module
