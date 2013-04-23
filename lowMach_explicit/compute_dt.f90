module compute_dt_module 

  use bl_types
  use multifab_module
  use probin_module, only : verbose

  implicit none

  private

  public :: compute_dt

contains

  subroutine compute_dt(lev, umac, dx, dtold, dt)

    use probin_module, only: max_dt_growth, cflfac

    type(multifab) , intent( in) :: umac(:)
    real(kind=dp_t), intent( in) :: dx(:)
    real(kind=dp_t), intent( in) :: dtold
    real(kind=dp_t), intent(out) :: dt
    integer        , intent( in) :: lev

    real(kind=dp_t), pointer :: ump(:,:,:,:), vmp(:,:,:,:), wmp(:,:,:,:)

    integer :: lo(get_dim(umac(1))),hi(get_dim(umac(1)))
    integer :: i,dm,ng_u

    real(kind=dp_t) :: dt_proc, dt_grid, dt_start

    ng_u = nghost(umac(1))
    dm = get_dim(umac(1))

    dt_proc  = 1.d20
    dt_start = 1.d20

    do i = 1, nfabs(umac(1))
       ump => dataptr(umac(1), i)
       vmp => dataptr(umac(2), i)
       lo =  lwb(get_box(umac(1), i))
       hi =  upb(get_box(umac(1), i))
       dt_grid = 1.d20
       select case (dm)
       case (2)
          call compute_dt_2d(ump(:,:,1,1), vmp(:,:,1,1), ng_u, lo, hi, dx, dt_grid)
       case (3)
          wmp => dataptr(umac(3), i)
          call compute_dt_3d(ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                             ng_u, lo, hi, dx, dt_grid)
       end select

       dt_proc = min(dt_grid, dt_proc)
    end do

    ! This sets dt to be the min of dt_proc over all processors.
    call parallel_reduce(dt ,dt_proc ,MPI_MIN)

    if (dt .eq. dt_start) then
       dt = min(dx(1),dx(2))
       if (dm .eq. 3) dt = min(dt,dx(3))
       if (parallel_IOProcessor() .and. verbose .ge. 1) then
          print*,'setting dt = cflfac*dx'
       end if
    end if

    dt = dt*cflfac

    if (dtold .gt. 0.0d0 .and. max_dt_growth*dtold .lt. dt) then
       dt = max_dt_growth*dtold
       if (parallel_IOProcessor() .and. verbose .ge. 1) then
          print*,'max_dt_growth limiting new dt to',dt
       end if
    end if

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,1000) lev,dt
    end if
1000 format("Computing dt at level ",i2," to be ... ",e15.8)

  end subroutine compute_dt

  subroutine compute_dt_2d(umac,vmac,ng_u,lo,hi,dx,dt)

    integer           , intent(in   ) :: lo(:), hi(:), ng_u
    real (kind = dp_t), intent(in   ) :: umac(lo(1)-ng_u:,lo(2)-ng_u:)  
    real (kind = dp_t), intent(in   ) :: vmac(lo(1)-ng_u:,lo(2)-ng_u:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real (kind = dp_t), intent(inout) :: dt

    !     Local variables
    real (kind = dp_t)  u,v
    real (kind = dp_t)  eps
    integer :: i, j

    eps = 1.0e-8

    u  = 0.0D0 
    v  = 0.0D0 

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          u  = max(u ,abs(umac(i,j)))
          v  = max(v ,abs(vmac(i,j)))
       enddo
    enddo

    if (u .gt. eps) dt = min(dt,dx(1)/u)
    if (v .gt. eps) dt = min(dt,dx(2)/v)

  end subroutine compute_dt_2d

  subroutine compute_dt_3d(umac,vmac,wmac,ng_u,lo,hi,dx,dt)

    integer           , intent(in   ) :: lo(:), hi(:), ng_u
    real (kind = dp_t), intent(in   ) :: umac(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real (kind = dp_t), intent(in   ) :: vmac(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real (kind = dp_t), intent(in   ) :: wmac(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real (kind = dp_t), intent(in   ) :: dx(:)
    real (kind = dp_t), intent(inout) :: dt

    !     Local variables
    real (kind = dp_t)  u,v,w
    real (kind = dp_t)  eps
    integer :: i, j, k

    eps = 1.0e-8

    u  = 0.0D0 
    v  = 0.0D0 
    w  = 0.0D0

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             u  = max(u ,abs(umac(i,j,k)))
             v  = max(v ,abs(vmac(i,j,k)))
             w  = max(w ,abs(wmac(i,j,k)))
          enddo
       enddo
    enddo

    if (u .gt. eps) dt = min(dt,dx(1)/u)
    if (v .gt. eps) dt = min(dt,dx(2)/v)
    if (w .gt. eps) dt = min(dt,dx(3)/w)

  end subroutine compute_dt_3d

end module compute_dt_module
