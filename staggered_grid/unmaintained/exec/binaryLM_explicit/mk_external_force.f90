module mk_external_force_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use multifab_fill_ghost_module

  implicit none

  private

  public :: mk_external_m_force, mk_external_s_force

contains

  ! Important note: For periodic boundaries, the mk_external_m_force routine should fill out
  ! *both* sides of the domain with values, even though this is duplicate information
  ! We ensure the two sides are bitwise identical, but low and/or high sides may win
  subroutine mk_external_m_force(mla,m_update,s,dx,time)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: m_update(:,:)
    type(multifab) , intent(in   ) :: s(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),time

    ! local
    integer :: i,n,ng_s,ng_u,dm,nlevs
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: fxp(:,:,:,:)
    real(kind=dp_t), pointer :: fyp(:,:,:,:)
    real(kind=dp_t), pointer :: fzp(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)

    nlevs = mla%nlevel
    dm = mla%dim

    ng_s  = s(1)%ng
    ng_u  = m_update(1,1)%ng

    do n=1,nlevs
       do i=1, nfabs(s(n))
          fxp => dataptr(m_update(n,1), i)
          fyp => dataptr(m_update(n,2), i)
          sp  => dataptr(s(n),i)
          lo = lwb(get_box(s(n), i))
          hi = upb(get_box(s(n), i))
          select case (dm)
          case (2)
             call mk_external_m_force_2d(fxp(:,:,1,1), fyp(:,:,1,1), ng_u, &
                                sp(:,:,1,:), ng_s, lo, hi, dx(n,:),time)
          case (3)
             fzp => dataptr(m_update(n,3), i)
             call mk_external_m_force_3d(fxp(:,:,:,1), fyp(:,:,:,1), fzp(:,:,:,1), ng_u, &
                                sp(:,:,:,:), ng_s, lo, hi, dx(n,:),time)
          end select
       end do

       ! For periodic boundaries, ensure the low and high side are consistent:
       ! Note: multifab_internal_sync compares the box number of the two boxes 
       ! with overlapping values and the data on the box with lower number wins. 
       do i=1,dm
          call multifab_internal_sync(m_update(n,i))
       end do

    enddo

  end subroutine mk_external_m_force

  subroutine mk_external_m_force_2d(m_updatex,m_updatey,ng_u,s,ng_s,lo,hi,dx,time)
    use exact_solutions_module, only: exact_m_force_2d
 
    use probin_module, only: grav

    integer        , intent(in   ) :: lo(:),hi(:),ng_u,ng_s
    real(kind=dp_t), intent(inout) :: m_updatex(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t), intent(inout) :: m_updatey(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t), intent(in   ) ::         s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),time

    ! local
    integer i,j

    do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          m_updatey(i,j) = m_updatey(i,j) - grav*(s(i,j,1)+s(i,j-1,1))/2.d0
       end do
    end do

    call exact_m_force_2d(m_updatex,m_updatey,ng_u,lo,hi,dx,time)

  end subroutine mk_external_m_force_2d

  subroutine mk_external_m_force_3d(m_updatex,m_updatey,m_updatez,ng_u,s, &
                       ng_s,lo,hi,dx,time)

    use exact_solutions_module, only: exact_m_force_3d
 
    use probin_module, only: grav

    integer        , intent(in   ) :: lo(:),hi(:),ng_u,ng_s
    real(kind=dp_t), intent(inout) :: m_updatex(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(inout) :: m_updatey(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(inout) :: m_updatez(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(in   ) ::         s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),time

    ! local
    integer i,j,k

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
             m_updatey(i,j,k) = m_updatey(i,j,k) - grav*(s(i,j,k,1)+s(i,j-1,k,1))/2.d0
          end do
       end do
    end do

    call exact_m_force_3d(m_updatex,m_updatey,m_updatez,ng_u,lo,hi,dx,time)

  end subroutine mk_external_m_force_3d

  subroutine mk_external_s_force(mla,s_update,s,dx,time)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: s_update(:)
    type(multifab) , intent(in   ) :: s(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),time

    ! local
    integer :: i,n,ng_s,ng_u,dm,nlevs
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)

    nlevs = mla%nlevel
    dm = mla%dim

    ng_s  = s(1)%ng
    ng_u  = s_update(1)%ng

    do n=1,nlevs
       do i = 1, nfabs(s(n))
          fp  => dataptr(s_update(n),i)
          sp  => dataptr(s(n),i)
          lo = lwb(get_box(s(n), i))
          hi = upb(get_box(s(n), i))
          select case (dm)
          case (2)
             call mk_external_s_force_2d(fp(:,:,1,:), sp(:,:,1,:), &
                                         ng_u, ng_s, lo, hi, dx(n,:), time)
          case (3)
             call mk_external_s_force_3d(fp(:,:,:,:), sp(:,:,:,:), &
                                         ng_u, ng_s, lo, hi, dx(n,:), time)
          end select
       end do
    enddo

  end subroutine mk_external_s_force

  subroutine mk_external_s_force_2d(s_update,s,ng_u,ng_s,lo,hi,dx,time)
    use exact_solutions_module, only: exact_s_force_2d
 
    integer        , intent(in   ) :: lo(:),hi(:),ng_u,ng_s
    real(kind=dp_t), intent(inout) :: s_update(lo(1)-ng_u:,lo(2)-ng_u:,:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),time
    
    call exact_s_force_2d(s_update,s,ng_u,ng_s,lo,hi,dx,time)

  end subroutine mk_external_s_force_2d

  subroutine mk_external_s_force_3d(s_update,s,ng_u,ng_s,lo,hi,dx,time)

    use exact_solutions_module, only: exact_s_force_3d

    integer        , intent(in   ) :: lo(:),hi(:),ng_u,ng_s
    real(kind=dp_t), intent(inout) :: s_update(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),time
    
    call exact_s_force_3d(s_update,s,ng_u,ng_s,lo,hi,dx,time)

  end subroutine mk_external_s_force_3d

end module mk_external_force_module
