module analysis_module

  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: print_errors, eos_check, sum_mass_momentum, staggered_norm

  contains

  subroutine print_errors(mla,sold,snew,mold,mnew,dx,time)

    use init_module
    use probin_module, only: nscal

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(in   ) :: mold(:,:)
    type(multifab) , intent(inout) :: mnew(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:), time

    integer :: n,nlevs,i,dm
    integer(kind=ll_t) :: n_cell
    real(kind=dp_t) :: norm,norm_stag(mla%dim,1:2)

    nlevs = mla%nlevel
    dm = mla%dim

    call initdata(mnew,snew,dx,mla,time)

    do n=1,nlevs
       do i=1,dm
          call multifab_sub_sub_c(mnew(n,i),1,mold(n,i),1,1,0)
       end do
       call multifab_sub_sub_c(snew(n),1,sold(n),1,nscal,0)
    end do

    ! compute L1 and L2 norms for staggered data
    call staggered_norm(mla,mnew,norm_stag)

    n_cell = multifab_volume(snew(1)) / snew(1)%nc

    norm = multifab_norm_inf_c(snew(1),1,1,all=.false.)
    if (parallel_IOProcessor()) print*,"L0 RHO  =",norm
    norm = multifab_norm_inf_c(snew(1),2,1,all=.false.)
    if (parallel_IOProcessor()) print*,"L0 RHOC =",norm
    norm = multifab_norm_inf_c(mnew(1,1),1,1,all=.false.)
    if (parallel_IOProcessor()) print*,"L0 MX   =",norm
    norm = multifab_norm_inf_c(mnew(1,2),1,1,all=.false.)
    if (parallel_IOProcessor()) print*,"L0 MY   =",norm
    if (dm .eq. 3) then
       norm = multifab_norm_inf_c(mnew(1,3),1,1,all=.false.)
       if (parallel_IOProcessor()) print*,"L0 MZ   =",norm
    end if

    if (parallel_IOProcessor()) print*,""

    norm = multifab_norm_l1_c(snew(1),1,1,all=.false.) / n_cell
    if (parallel_IOProcessor()) print*,"L1 RHO  =",norm
    norm = multifab_norm_l1_c(snew(1),2,1,all=.false.)
    if (parallel_IOProcessor()) print*,"L1 RHOC =",norm / n_cell
    if (parallel_IOProcessor()) print*,"L1 MX   =",norm_stag(1,1) / n_cell
    if (parallel_IOProcessor()) print*,"L1 MY   =",norm_stag(2,1) / n_cell
    if (dm .eq. 3) then
       if (parallel_IOProcessor()) print*,"L1 MZ   =",norm_stag(3,1) / n_cell
    end if

    if (parallel_IOProcessor()) print*,""

    norm = multifab_norm_l2_c(snew(1),1,1,all=.false.)
    if (parallel_IOProcessor()) print*,"L2 RHO  =",norm / sqrt(dble(n_cell))
    norm = multifab_norm_l2_c(snew(1),2,1,all=.false.)
    if (parallel_IOProcessor()) print*,"L2 RHOC =",norm / sqrt(dble(n_cell))
    if (parallel_IOProcessor()) print*,"L2 MX   =",norm_stag(1,2) / sqrt(dble(n_cell))
    if (parallel_IOProcessor()) print*,"L2 MY   =",norm_stag(2,2) / sqrt(dble(n_cell))
    if (dm .eq. 3) then
       if (parallel_IOProcessor()) print*,"L2 MZ   =",norm_stag(3,2) / sqrt(dble(n_cell))
    end if

    if (parallel_IOProcessor()) print*,""

  end subroutine print_errors

  subroutine eos_check(mla,s)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: s(:)

    ! local
    integer i,n,dm,nlevs,ng_s
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: sp(:,:,:,:)

    real(kind=dp_t) :: eos_error, eos_error_grid, eos_error_proc

    nlevs = mla%nlevel
    dm = mla%dim

    ng_s = s(1)%ng
    eos_error_proc = -1.d20
    do n=1,nlevs
       do i=1,nfabs(s(n))
         sp => dataptr(s(n), i)
         lo =  lwb(get_box(s(n), i))
         hi =  upb(get_box(s(n), i))
         eos_error_grid = -1.d20
         select case (dm)
         case (2)
            call eos_check_2d(sp(:,:,1,:),ng_s,eos_error_grid,lo,hi)
         case (3)
            call eos_check_3d(sp(:,:,:,:),ng_s,eos_error_grid,lo,hi)
         end select
         eos_error_proc = max(eos_error_grid, eos_error_proc)
      end do
   end do


   ! This sets eos_error to be the max of eos_error_proc over all processors.
   call parallel_reduce(eos_error, eos_error_proc, MPI_MAX)

   if (parallel_IOProcessor()) then
      print*,"EOS ERROR in L1 norm: ",eos_error
      print*,""
   end if

  end subroutine eos_check

  subroutine eos_check_2d(s,ng_s,eos_error,lo,hi)

    use probin_module, only: rhobar

    integer        , intent(in   ) :: lo(:), hi(:), ng_s
    real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(inout) :: eos_error

    ! local
    integer :: i,j

    real(kind=dp_t) :: error

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          error = abs( (s(i,j,2)/rhobar(1) + (s(i,j,1)-s(i,j,2))/rhobar(2) ) - 1.d0)
          eos_error = max(eos_error,error)

       end do
    end do

  end subroutine eos_check_2d

  subroutine eos_check_3d(s,ng_s,eos_error,lo,hi)

    use probin_module, only: rhobar

    integer        , intent(in   ) :: lo(:), hi(:), ng_s
    real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(inout) :: eos_error

    ! local
    integer :: i,j,k

    real(kind=dp_t) :: error

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             error = abs( (s(i,j,k,2)/rhobar(1) + (s(i,j,k,1)-s(i,j,k,2))/rhobar(2) ) - 1.d0)
             eos_error = max(eos_error,error)
             
          end do
       end do
    end do

  end subroutine eos_check_3d

  subroutine sum_mass_momentum(mla,cons,m,av_mass,av_momentum)

    use probin_module, only: nscal

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: cons(:)
    type(multifab) , intent(in   ) ::    m(:,:)
    real(kind=dp_t), intent(out), optional :: av_mass(nscal), av_momentum(mla%dim) ! Level 1 only

    ! local
    integer :: i,comp,n,dm,nlevs,ng_c,ng_m
    integer(kind=ll_t) :: n_cell
    integer :: lo(mla%dim), hi(mla%dim)

    real(kind=dp_t) :: mass_tot(nscal) , mass_lev(nscal),  mass_proc(nscal),  mass_grid(nscal)
    real(kind=dp_t) :: mom_tot(mla%dim), mom_lev(mla%dim), mom_proc(mla%dim), mom_grid(mla%dim)

    real(kind=dp_t), pointer ::  cp(:,:,:,:)
    real(kind=dp_t), pointer :: mxp(:,:,:,:)
    real(kind=dp_t), pointer :: myp(:,:,:,:)
    real(kind=dp_t), pointer :: mzp(:,:,:,:)

    mass_tot = 0.d0
    mass_lev = 0.d0
    mass_proc = 0.d0

    mom_tot = 0.d0
    mom_lev = 0.d0
    mom_proc = 0.d0

    ng_c = cons(1)%ng
    ng_m = m(1,1)%ng
    
    nlevs = mla%nlevel
    dm = mla%dim
    
    n_cell = multifab_volume(cons(1)) / cons(1)%ng
    
    if (nlevs .gt. 1) then
       call bl_error('sum_mass_momentum not written for multilevel yet')
    end if

    do i=1,nfabs(cons(1))
       cp  => dataptr(cons(1), i)
       mxp => dataptr(m(1,1), i)
       myp => dataptr(m(1,2), i)
       lo = lwb(get_box(cons(1), i))
       hi = upb(get_box(cons(1), i))
       mass_grid = 0.d0
       mom_grid  = 0.d0
       select case (dm)
       case (2)
          call sum_mass_momentum_2d(cp(:,:,1,:), ng_c, mxp(:,:,1,1), myp(:,:,1,1), ng_m, &
                                    lo, hi, mass_grid, mom_grid)
       case (3)
          mzp => dataptr(m(1,3), i)
          call sum_mass_momentum_3d(cp(:,:,:,:), ng_c, &
                                    mxp(:,:,:,1), myp(:,:,:,1), mzp(:,:,:,1), ng_m, &
                                    lo, hi, mass_grid, mom_grid)

       end select

       mass_proc(1:nscal) = mass_proc(1:nscal) + mass_grid(1:nscal)
       mom_proc(1:dm)     = mom_proc(1:dm)     + mom_grid(1:dm)

    end do

    call parallel_reduce(mass_lev(1:nscal), mass_proc(1:nscal), MPI_SUM)
    call parallel_reduce(mom_lev(1:dm)    , mom_proc(1:dm)    , MPI_SUM)

    if (parallel_IOProcessor()) then
       write(*,"(A,100G17.9)") "CONSERVE: <rho_i>=", mass_lev(1:nscal)/n_cell
       write(*,"(A,100G17.9)") "CONSERVE: <mom_k>=", mom_lev(1:dm)/n_cell
       write(*,*)
    end if
        
    if(present(av_mass)) av_mass=mass_lev(1:nscal)/n_cell
    if(present(av_momentum)) av_momentum=mom_lev(1:dm)/n_cell
    
  end subroutine sum_mass_momentum

  subroutine sum_mass_momentum_2d(cons,ng_c,mx,my,ng_m,lo,hi,mass,mom)

    use probin_module, only: nscal

    integer        , intent(in   ) :: lo(:), hi(:), ng_c, ng_m
    real(kind=dp_t), intent(in   ) :: cons(lo(1)-ng_c:,lo(2)-ng_c:,:)
    real(kind=dp_t), intent(in   ) ::   mx(lo(1)-ng_m:,lo(2)-ng_m:)
    real(kind=dp_t), intent(in   ) ::   my(lo(1)-ng_m:,lo(2)-ng_m:)
    real(kind=dp_t), intent(inout) :: mass(:),mom(:)

    ! local
    integer :: i,j

    ! mass
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          mass(1:nscal) = mass(1:nscal) + cons(i,j,1:nscal)
       end do
    end do

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

  end subroutine sum_mass_momentum_2d

  subroutine sum_mass_momentum_3d(cons,ng_c,mx,my,mz,ng_m,lo,hi,mass,mom)

    use probin_module, only: nscal

    integer        , intent(in   ) :: lo(:), hi(:), ng_c, ng_m
    real(kind=dp_t), intent(in   ) :: cons(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:,:)
    real(kind=dp_t), intent(in   ) ::   mx(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real(kind=dp_t), intent(in   ) ::   my(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real(kind=dp_t), intent(in   ) ::   mz(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real(kind=dp_t), intent(inout) :: mass(:),mom(:)

    ! local
    integer :: i,j,k

    ! mass
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             mass(1:nscal) = mass(1:nscal) + cons(i,j,k,1:nscal)
          end do
       end do
    end do

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

  end subroutine sum_mass_momentum_3d

  subroutine staggered_norm(mla,m,norm)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: m(:,:)
    real(kind=dp_t), intent(inout) :: norm(:,:)

    ! local
    integer :: i,comp,n,dm,nlevs,ng_m
    integer :: lo(mla%dim), hi(mla%dim)

    real(kind=dp_t) :: norm_proc(mla%dim,2), norm_grid(mla%dim,2)

    real(kind=dp_t), pointer :: mxp(:,:,:,:)
    real(kind=dp_t), pointer :: myp(:,:,:,:)
    real(kind=dp_t), pointer :: mzp(:,:,:,:)

    norm_proc = 0.d0

    ng_m = m(1,1)%ng
    
    nlevs = mla%nlevel
    dm = mla%dim

    if (nlevs .gt. 1) then
       call bl_error('staggered_norm not written for multilevel yet')
    end if

    do i=1,nfabs(m(1,1))
       mxp => dataptr(m(1,1), i)
       myp => dataptr(m(1,2), i)
       lo = lwb(get_box(m(1,1), i))
       hi = upb(get_box(m(1,1), i))
       norm_grid = 0.d0
       select case (dm)
       case (2)
          call staggered_norm_2d(mxp(:,:,1,1), myp(:,:,1,1), ng_m, &
                                 lo, hi, norm_grid)
       case (3)
          mzp => dataptr(m(1,3), i)
          call staggered_norm_3d(mxp(:,:,:,1), myp(:,:,:,1), mzp(:,:,:,1), ng_m, &
                                 lo, hi, norm_grid)
       end select

       norm_proc(1:dm,:) = norm_proc(1:dm,:) + norm_grid(1:dm,:)

    end do

    do comp=1,dm
       call parallel_reduce(norm(comp,:), norm_proc(comp,:), MPI_SUM)
       norm(comp,2) = sqrt(norm(comp,2))
    end do
    
  end subroutine staggered_norm

  subroutine staggered_norm_2d(mx,my,ng_m,lo,hi,norm)

    use probin_module, only: nscal

    integer        , intent(in   ) :: lo(:), hi(:), ng_m
    real(kind=dp_t), intent(in   ) ::   mx(lo(1)-ng_m:,lo(2)-ng_m:)
    real(kind=dp_t), intent(in   ) ::   my(lo(1)-ng_m:,lo(2)-ng_m:)
    real(kind=dp_t), intent(inout) :: norm(:,:)

    ! local
    integer :: i,j

    ! mx, interior cells
       do j=lo(2),hi(2)
          do i=lo(1)+1,hi(1)
             norm(1,1) = norm(1,1) + abs(mx(i,j))
             norm(1,2) = norm(1,2) + mx(i,j)**2
          end do
       end do

    ! mx, boundary cells
       do j=lo(2),hi(2)
          norm(1,1) = norm(1,1) + 0.5d0*abs(mx(lo(1),j))
          norm(1,2) = norm(1,2) + 0.5d0*mx(lo(1),j)**2
          norm(1,1) = norm(1,1) + 0.5d0*abs(mx(hi(1)+1,j))
          norm(1,2) = norm(1,2) + 0.5d0*mx(hi(1)+1,j)**2
       end do

    ! my, interior cells
       do j=lo(2)+1,hi(2)
          do i=lo(1),hi(1)
             norm(2,1) = norm(2,1) + abs(my(i,j))
             norm(2,2) = norm(2,2) + my(i,j)**2
          end do
       end do

    ! my, boundary cells
       do i=lo(1),hi(1)
          norm(2,1) = norm(2,1) + 0.5d0*abs(my(i,lo(2)))
          norm(2,2) = norm(2,2) + 0.5d0*my(i,lo(2))**2
          norm(2,1) = norm(2,1) + 0.5d0*abs(my(i,hi(2)+1))
          norm(2,2) = norm(2,2) + 0.5d0*my(i,hi(2)+1)**2
       end do

  end subroutine staggered_norm_2d

  subroutine staggered_norm_3d(mx,my,mz,ng_m,lo,hi,norm)

    use probin_module, only: nscal

    integer        , intent(in   ) :: lo(:), hi(:), ng_m
    real(kind=dp_t), intent(in   ) ::   mx(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real(kind=dp_t), intent(in   ) ::   my(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real(kind=dp_t), intent(in   ) ::   mz(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real(kind=dp_t), intent(inout) :: norm(:,:)

    ! local
    integer :: i,j,k

    ! mx, interior cells
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1)+1,hi(1)
             norm(1,1) = norm(1,1) + abs(mx(i,j,k))
             norm(1,2) = norm(1,2) + mx(i,j,k)**2
          end do
       end do
    end do

    ! mx, boundary cells
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          norm(1,1) = norm(1,1) + 0.5d0*abs(mx(lo(1),j,k))
          norm(1,2) = norm(1,2) + 0.5d0*mx(lo(1),j,k)**2
          norm(1,1) = norm(1,1) + 0.5d0*abs(mx(hi(1)+1,j,k))
          norm(1,2) = norm(1,2) + 0.5d0*mx(hi(1)+1,j,k)**2
       end do
    end do

    ! my, interior cells
    do k=lo(3),hi(3)
       do j=lo(2)+1,hi(2)
          do i=lo(1),hi(1)
             norm(2,1) = norm(2,1) + abs(my(i,j,k))
             norm(2,2) = norm(2,2) + my(i,j,k)**2
          end do
       end do
    end do

    ! my, boundary cells
    do k=lo(3),hi(3)
       do i=lo(1),hi(1)
          norm(2,1) = norm(2,1) + 0.5d0*abs(my(i,lo(2),k))
          norm(2,2) = norm(2,2) + 0.5d0*my(i,lo(2),k)**2
          norm(2,1) = norm(2,1) + 0.5d0*abs(my(i,hi(2)+1,k))
          norm(2,2) = norm(2,2) + 0.5d0*my(i,hi(2)+1,k)**2
       end do
    end do

    ! mz, interior cells
    do k=lo(3)+1,hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             norm(3,1) = norm(3,1) + abs(mz(i,j,k))
             norm(3,2) = norm(3,2) + mz(i,j,k)**2
          end do
       end do
    end do

    ! mz, boundary cells
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          norm(3,1) = norm(3,1) + 0.5d0*abs(mz(i,j,lo(3)))
          norm(3,2) = norm(3,2) + 0.5d0*mz(i,j,lo(3))**2
          norm(3,1) = norm(3,1) + 0.5d0*abs(mz(i,j,hi(3)+1))
          norm(3,2) = norm(3,2) + 0.5d0*mz(i,j,hi(3)+1)**2
       end do
    end do

  end subroutine staggered_norm_3d

end module analysis_module
