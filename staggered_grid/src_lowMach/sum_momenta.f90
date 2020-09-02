module sum_momenta_module

  use multifab_module
  use ml_layout_module
  use probin_common_module, only: total_volume

  implicit none

  private

  public :: sum_momenta, sum_kinetic

contains

  subroutine sum_momenta(mla,m,av_m)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: m(:,:)
    real(kind=dp_t), intent(inout), optional :: av_m(:)

    ! local
    integer :: i,dm,nlevs,ng_m
    integer :: lo(mla%dim), hi(mla%dim)

    real(kind=dp_t) :: mom_tot(mla%dim), mom_lev(mla%dim), mom_proc(mla%dim), mom_grid(mla%dim)

    real(kind=dp_t), pointer :: mxp(:,:,:,:)
    real(kind=dp_t), pointer :: myp(:,:,:,:)
    real(kind=dp_t), pointer :: mzp(:,:,:,:)

    type(bl_prof_timer),save :: bpt

    call build(bpt,"sum_momenta")

    mom_tot = 0.d0
    mom_lev = 0.d0
    mom_proc = 0.d0

    ng_m = m(1,1)%ng
    
    nlevs = mla%nlevel
    dm = mla%dim
    
    if (nlevs .gt. 1) then
       call bl_error('sum_momenta not written for multilevel yet')
    end if

    do i=1,nfabs(m(1,1))
       mxp => dataptr(m(1,1), i)
       myp => dataptr(m(1,2), i)
       lo = lwb(get_box(m(1,1), i))
       hi = upb(get_box(m(1,1), i))
       mom_grid  = 0.d0
       select case (dm)
       case (2)
          call sum_momenta_2d(mxp(:,:,1,1), myp(:,:,1,1), ng_m, &
                              lo, hi, mom_grid)
       case (3)
          mzp => dataptr(m(1,3), i)
          call sum_momenta_3d(mxp(:,:,:,1), myp(:,:,:,1), mzp(:,:,:,1), ng_m, &
                              lo, hi, mom_grid)

       end select

       mom_proc(1:dm) = mom_proc(1:dm) + mom_grid(1:dm)

    end do

    call parallel_reduce(mom_lev(1:dm)    , mom_proc(1:dm)    , MPI_SUM)

    if (parallel_IOProcessor()) then
       write(*,"(A,100G17.9)") "CONSERVE: <mom_k>=", mom_lev(1:dm)/total_volume
    end if
        
    if(present(av_m)) av_m=mom_lev(1:dm)/total_volume
    
    call destroy(bpt)

  end subroutine sum_momenta

  subroutine sum_momenta_2d(mx,my,ng_m,lo,hi,mom)

    integer        , intent(in   ) :: lo(:), hi(:), ng_m
    real(kind=dp_t), intent(in   ) :: mx(lo(1)-ng_m:,lo(2)-ng_m:)
    real(kind=dp_t), intent(in   ) :: my(lo(1)-ng_m:,lo(2)-ng_m:)
    real(kind=dp_t), intent(inout) :: mom(:)

    ! local
    integer :: i,j

    ! mx, interior cells
    do j=lo(2),hi(2)
       do i=lo(1)+1,hi(1)
          mom(1) = mom(1) + mx(i,j)
       end do
    end do

    ! mx, boundary cells
    do j=lo(2),hi(2)
       mom(1) = mom(1) + 0.5d0*mx(lo(1)  ,j)
       mom(1) = mom(1) + 0.5d0*mx(hi(1)+1,j)
    end do

    ! my, interior cells
    do j=lo(2)+1,hi(2)
       do i=lo(1),hi(1)
          mom(2) = mom(2) + my(i,j)
       end do
    end do

    ! my, boundary cells
    do i=lo(1),hi(1)
       mom(2) = mom(2) + 0.5d0*my(i,lo(2)  )
       mom(2) = mom(2) + 0.5d0*my(i,hi(2)+1)
    end do

  end subroutine sum_momenta_2d

  subroutine sum_momenta_3d(mx,my,mz,ng_m,lo,hi,mom)

    integer        , intent(in   ) :: lo(:), hi(:), ng_m
    real(kind=dp_t), intent(in   ) :: mx(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real(kind=dp_t), intent(in   ) :: my(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real(kind=dp_t), intent(in   ) :: mz(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real(kind=dp_t), intent(inout) :: mom(:)

    ! local
    integer :: i,j,k

    ! mx, interior cells
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1)+1,hi(1)
             mom(1) = mom(1) + mx(i,j,k)
          end do
       end do
    end do

    ! mx, boundary cells
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          mom(1) = mom(1) + 0.5d0*mx(lo(1)  ,j,k)
          mom(1) = mom(1) + 0.5d0*mx(hi(1)+1,j,k)
       end do
    end do

    ! my, interior cells
    do k=lo(3),hi(3)
       do j=lo(2)+1,hi(2)
          do i=lo(1),hi(1)
             mom(2) = mom(2) + my(i,j,k)
          end do
       end do
    end do

    ! my, boundary cells
    do k=lo(3),hi(3)
       do i=lo(1),hi(1)
          mom(2) = mom(2) + 0.5d0*my(i,lo(2)  ,k)
          mom(2) = mom(2) + 0.5d0*my(i,hi(2)+1,k)
       end do
    end do

    ! mz, interior cells
    do k=lo(3)+1,hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             mom(3) = mom(3) + mz(i,j,k)
          end do
       end do
    end do

    ! mz, boundary cells
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          mom(3) = mom(3) + 0.5d0*mz(i,j,lo(3)  )
          mom(3) = mom(3) + 0.5d0*mz(i,j,hi(3)+1)
       end do
    end do

  end subroutine sum_momenta_3d

  subroutine sum_kinetic(mla,umac,rho_fc)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: umac(:,:)
    type(multifab) , intent(in   ) :: rho_fc(:,:)

    ! local
    integer :: i,dm,nlevs,ng_u,ng_r
    integer :: lo(mla%dim), hi(mla%dim)

    real(kind=dp_t) :: ke_tot(mla%dim), ke_lev(mla%dim), ke_proc(mla%dim)

    real(kind=dp_t), pointer :: mxp(:,:,:,:)
    real(kind=dp_t), pointer :: myp(:,:,:,:)
    real(kind=dp_t), pointer :: mzp(:,:,:,:)

    real(kind=dp_t), pointer :: rxp(:,:,:,:)
    real(kind=dp_t), pointer :: ryp(:,:,:,:)
    real(kind=dp_t), pointer :: rzp(:,:,:,:)
    
    type(bl_prof_timer),save :: bpt

    call build(bpt,"sum_kinetic")

    ke_tot = 0.d0
    ke_lev = 0.d0
    ke_proc = 0.d0

    ng_u = umac(1,1)%ng
    ng_r = rho_fc(1,1)%ng
    
    nlevs = mla%nlevel
    dm = mla%dim
    
    if (nlevs .gt. 1) then
       call bl_error('sum_kinetic not written for multilevel yet')
    end if

    do i=1,nfabs(umac(1,1))
       mxp => dataptr(umac(1,1), i)
       myp => dataptr(umac(1,2), i)
       rxp => dataptr(rho_fc(1,1), i)
       ryp => dataptr(rho_fc(1,2), i)
       lo = lwb(get_box(umac(1,1), i))
       hi = upb(get_box(umac(1,1), i))
       select case (dm)
       case (2)
          call sum_kinetic_2d(mxp(:,:,1,1), myp(:,:,1,1), ng_u, &
                              rxp(:,:,1,1), ryp(:,:,1,1), ng_r, &
                              lo, hi, ke_proc)
       case (3)
          mzp => dataptr(umac(1,3), i)
          rzp => dataptr(rho_fc(1,3), i)
          call sum_kinetic_3d(mxp(:,:,:,1), myp(:,:,:,1), mzp(:,:,:,1), ng_u, &
                              rxp(:,:,:,1), ryp(:,:,:,1), rzp(:,:,:,1), ng_r, &
                              lo, hi, ke_proc)

       end select

    end do

    call parallel_reduce(ke_lev(1:dm), ke_proc(1:dm), MPI_SUM)

    if (parallel_IOProcessor()) then
       print*,'Total Kinetic Energy',sum(ke_lev(1:dm))
    end if
            
    call destroy(bpt)

  end subroutine sum_kinetic

  subroutine sum_kinetic_2d(ux,uy,ng_u,rx,ry,ng_r,lo,hi,ke)

    integer        , intent(in   ) :: lo(:), hi(:), ng_u, ng_r
    real(kind=dp_t), intent(in   ) :: ux(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t), intent(in   ) :: uy(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t), intent(in   ) :: rx(lo(1)-ng_r:,lo(2)-ng_r:)
    real(kind=dp_t), intent(in   ) :: ry(lo(1)-ng_r:,lo(2)-ng_r:)
    real(kind=dp_t), intent(inout) :: ke(:)

    ! local
    integer :: i,j

    ! x, interior cells
    do j=lo(2),hi(2)
       do i=lo(1)+1,hi(1)
          ke(1) = ke(1) + 0.5d0*rx(i,j)*ux(i,j)**2
       end do
    end do

    ! x, boundary cells
    do j=lo(2),hi(2)
       ke(1) = ke(1) + 0.25d0*rx(lo(1)  ,j)*ux(lo(1)  ,j)**2
       ke(1) = ke(1) + 0.25d0*rx(hi(1)+1,j)*ux(hi(1)+1,j)**2
    end do

    ! y, interior cells
    do j=lo(2)+1,hi(2)
       do i=lo(1),hi(1)
          ke(2) = ke(2) + 0.5d0*ry(i,j)*uy(i,j)**2
       end do
    end do

    ! y, boundary cells
    do i=lo(1),hi(1)
       ke(2) = ke(2) + 0.25d0*ry(i,lo(2)  )*uy(i,lo(2)  )**2
       ke(2) = ke(2) + 0.25d0*ry(i,hi(2)+1)*uy(i,hi(2)+1)**2
    end do

  end subroutine sum_kinetic_2d

  subroutine sum_kinetic_3d(ux,uy,uz,ng_u,rx,ry,rz,ng_r,lo,hi,ke)

    integer        , intent(in   ) :: lo(:), hi(:), ng_u, ng_r
    real(kind=dp_t), intent(in   ) :: ux(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(in   ) :: uy(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(in   ) :: uz(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(in   ) :: rx(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:)
    real(kind=dp_t), intent(in   ) :: ry(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:)
    real(kind=dp_t), intent(in   ) :: rz(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:)
    real(kind=dp_t), intent(inout) :: ke(:)

    ! local
    integer :: i,j,k

    ! x, interior cells
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1)+1,hi(1)
             ke(1) = ke(1) + 0.5d0*rx(i,j,k)*ux(i,j,k)**2
          end do
       end do
    end do

    ! x, boundary cells
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          ke(1) = ke(1) + 0.25d0*rx(lo(1)  ,j,k)*ux(lo(1)  ,j,k)**2
          ke(1) = ke(1) + 0.25d0*rx(hi(1)+1,j,k)*ux(hi(1)+1,j,k)**2
       end do
    end do

    ! y, interior cells
    do k=lo(3),hi(3)
       do j=lo(2)+1,hi(2)
          do i=lo(1),hi(1)
             ke(2) = ke(2) + 0.5d0*ry(i,j,k)*uy(i,j,k)**2
          end do
       end do
    end do

    ! y, boundary cells
    do k=lo(3),hi(3)
       do i=lo(1),hi(1)
          ke(2) = ke(2) + 0.25d0*ry(i,lo(2)  ,k)*uy(i,lo(2)  ,k)**2
          ke(2) = ke(2) + 0.25d0*ry(i,hi(2)+1,k)*uy(i,hi(2)+1,k)**2
       end do
    end do

    ! z, interior cells
    do k=lo(3)+1,hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ke(3) = ke(3) + 0.5d0*rz(i,j,k)*uz(i,j,k)**2
          end do
       end do
    end do

    ! z, boundary cells
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          ke(3) = ke(3) + 0.25d0*rz(i,j,lo(3)  )*uz(i,j,lo(3)  )**2
          ke(3) = ke(3) + 0.25d0*rz(i,j,hi(3)+1)*uz(i,j,hi(3)+1)**2
       end do
    end do

  end subroutine sum_kinetic_3d
  
end module sum_momenta_module
