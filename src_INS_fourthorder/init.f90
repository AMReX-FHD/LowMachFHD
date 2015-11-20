module init_module

  use bl_types
  use bl_constants_module
  use ml_layout_module
  use probin_common_module, only: prob_lo, prob_hi
  
  implicit none

  private

  public :: init

contains

  subroutine init(mla,umac,dx,time)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:), time
 
    ! local variables
    integer                        :: lo(mla%dim), hi(mla%dim)
    integer                        :: i, dm, n, nlevs, ng_u
    real(kind=dp_t), pointer       :: up(:,:,:,:)
    real(kind=dp_t), pointer       :: vp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "init")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_u = umac(1,1)%ng

    do n=1,nlevs
       ! looping over boxes 
       do i=1,nfabs(umac(n,1))
          up => dataptr(umac(n,1),i)
          vp => dataptr(umac(n,2),i)
          lo = lwb(get_box(umac(n,1),i))
          hi = upb(get_box(umac(n,1),i))
          select case (dm)
          case (2)
             call init_2d(up(:,:,1,1),vp(:,:,1,1),ng_u,lo,hi,dx(n,:),time)
          case (3)
          end select
       end do
    end do

    do n=1,nlevs
       do i=1,dm
          call multifab_fill_boundary(umac(n,i))
       end do
    end do

    call destroy(bpt)

  end subroutine init

  subroutine init_2d(umac,vmac,ng_u,lo,hi,dx,time)

    integer          :: lo(:), hi(:), ng_u
    real(kind=dp_t)  :: umac(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t)  :: vmac(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t)  :: dx(:), time

    ! local varables
    integer         :: i,j
    real(kind=dp_t) :: x,y,xlo,xhi,ylo,yhi,L(2)

    L(1:2) = prob_hi(1:2) - prob_lo(1:2)

    ! x-velocity
    do j=lo(2),hi(2)
       y   = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
       ylo = prob_lo(2) + (dble(j)      )*dx(2)
       yhi = prob_lo(2) + (dble(j)+1.d0 )*dx(2)
       do i=lo(1),hi(1)+1
          x = prob_lo(1) + dble(i)*dx(1)

          umac(i,j) = -2.d0*cos(2.d0*M_PI*x/L(1))*sin(2.d0*M_PI*y/L(2))

       end do
    end do

    ! y-velocity
    do j=lo(2),hi(2)+1
       y = prob_lo(2) + dble(j)*dx(2)
       do i=lo(1),hi(1)
          x   = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
          xlo = prob_lo(1) + (dble(i)      )*dx(1)
          xhi = prob_lo(1) + (dble(i)+1.d0 )*dx(1)

          vmac(i,j) = 2.d0*sin(2.d0*M_PI*x/L(1))*cos(2.d0*M_PI*y/L(2))

       end do
    end do

  end subroutine init_2d

end module init_module
