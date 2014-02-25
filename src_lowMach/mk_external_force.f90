module mk_external_force_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use ml_restriction_module
  use multifab_fill_ghost_module
  use bl_constants_module
  use probin_lowmach_module, only: diff_coef
  use probin_common_module, only: prob_type

  implicit none

  private

  public :: mk_external_s_force

contains

  subroutine mk_external_s_force(mla,force,dx,time,comp)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: force(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),time
    integer,         intent(in   ) :: comp

    ! local
    integer :: i,n,ng_u,dm,nlevs
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: fp(:,:,:,:)

    nlevs = mla%nlevel
    dm = mla%dim

    ng_u  = force(1)%ng

    do n=1,nlevs
       do i = 1, nfabs(force(n))
          fp  => dataptr(force(n),i)
          lo = lwb(get_box(force(n), i))
          hi = upb(get_box(force(n), i))
          select case (dm)
          case (2)
             call mk_external_s_force_2d(fp(:,:,1,comp), &
                                         ng_u, lo, hi, dx(n,:), time)
          case (3)
          end select
       end do
    enddo

  end subroutine mk_external_s_force

  subroutine mk_external_s_force_2d(force,ng_u,lo,hi,dx,time)

    integer        , intent(in   ) :: lo(:),hi(:),ng_u
    real(kind=dp_t), intent(inout) :: force(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t), intent(in   ) :: dx(:),time
    

    ! local
    integer :: i,j
    real(kind=dp_t) :: x, pfac, ucst, vcst, freq, pfreq

    select case (prob_type)
    case (6) ! Taylor (traveling wave) vortices

       ucst = 0.75d0
       vcst = 0.75d0
       freq  = 2.d0*M_PI
       pfreq  = 4.d0*M_PI
       pfac = dexp(-4.0d0*freq*freq*time*diff_coef)/64.0d0

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             x = dx(1) * (dble(i) + half) - ucst*time
             force(i,j) = force(i,j) - &
                       pfac*pfreq*sin(pfreq*x)
          enddo
       enddo

    end select

  end subroutine mk_external_s_force_2d

end module mk_external_force_module
