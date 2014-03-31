module mk_external_force_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use bl_constants_module
  use probin_binarylm_module, only: diff_coef, visc_coef
  use probin_common_module, only: prob_type

  implicit none

  private

  public :: mk_external_m_force, mk_external_s_force

contains



  ! Important note: For periodic boundaries, the mk_external_m_force routine should fill out
  ! *both* sides of the domain with values, even though this is duplicate information
  ! We ensure the two sides are bitwise identical, but low and/or high sides may win
  subroutine mk_external_m_force(mla,m_force,dx,time)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: m_force(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),time

    ! local
    integer :: i,n,ng_s,ng_u,dm,nlevs
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: fxp(:,:,:,:)
    real(kind=dp_t), pointer :: fyp(:,:,:,:)
    real(kind=dp_t), pointer :: fzp(:,:,:,:)

    nlevs = mla%nlevel
    dm = mla%dim

    ng_u  = m_force(1,1)%ng

    do n=1,nlevs
       do i=1, nfabs(m_force(n,1))
          fxp => dataptr(m_force(n,1), i)
          fyp => dataptr(m_force(n,2), i)
          lo = lwb(get_box(m_force(n,1), i))
          hi = upb(get_box(m_force(n,1), i))
          select case (dm)
          case (2)
             call mk_external_m_force_2d(fxp(:,:,1,1), fyp(:,:,1,1), ng_u, &
                                         lo, hi, dx(n,:), time)
          case (3)
             fzp => dataptr(m_force(n,3), i)
             call mk_external_m_force_3d(fxp(:,:,:,1), fyp(:,:,:,1), fzp(:,:,:,1), ng_u, &
                                         lo, hi, dx(n,:), time)
          end select
       end do

       ! For periodic boundaries, ensure the low and high side are consistent:
       ! Note: multifab_internal_sync compares the box number of the two boxes 
       ! with overlapping values and the data on the box with lower number wins. 
       do i=1,dm
          call multifab_internal_sync(m_force(n,i))
       end do

    enddo

  end subroutine mk_external_m_force

  subroutine mk_external_m_force_2d(m_forcex,m_forcey,ng_u,lo,hi,dx,time) 

    integer        , intent(in   ) :: lo(:),hi(:),ng_u
    real(kind=dp_t), intent(inout) :: m_forcex(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t), intent(inout) :: m_forcey(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t), intent(in   ) :: dx(:),time

    ! local
    integer i,j
    real(kind=dp_t) :: x,y,ucst,vcst,freq,ufac

    select case (prob_type)
    case (6) ! Taylor (traveling wave) vortices

       ucst = 0.75d0
       vcst = 0.75d0

       freq = 2.d0*M_PI

       ufac = dexp(-2.d0*freq*freq*time*visc_coef)

       do j=lo(2),hi(2)
          y = dx(2)*(dble(j)+half)
          do i=lo(1),hi(1)+1
             x = dx(1)*dble(i)
             m_forcex(i,j) = m_forcex(i,j) - (M_PI/8.d0)*cos(freq*(time*ucst-x))*ufac &
                  *(16.d0*M_PI*visc_coef*sin(freq*(time*vcst-y)) + ufac*sin(freq*(time*ucst-x)))
          end do
       end do

       do j=lo(2),hi(2)+1
          y = dx(2)*dble(j)
          do i=lo(1),hi(1)
             x = dx(1)*(dble(i)+half)
             m_forcey(i,j) = m_forcey(i,j) - (M_PI/8.d0)*cos(freq*(time*vcst-y))*ufac &
                  *(-16.d0*M_PI*visc_coef*sin(freq*(time*ucst-x)) + ufac*sin(freq*(time*vcst-y)))
          end do
       end do

    end select

  end subroutine mk_external_m_force_2d

  subroutine mk_external_m_force_3d(m_forcex,m_forcey,m_forcez,ng_u,lo,hi,dx,time)

    integer        , intent(in   ) :: lo(:),hi(:),ng_u
    real(kind=dp_t), intent(inout) :: m_forcex(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(inout) :: m_forcey(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(inout) :: m_forcez(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(in   ) :: dx(:),time

    ! local
    integer i,j,k

  end subroutine mk_external_m_force_3d

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
             force(i,j) = force(i,j) - pfac*pfreq*sin(pfreq*x)
          enddo
       enddo

    end select

  end subroutine mk_external_s_force_2d

end module mk_external_force_module
