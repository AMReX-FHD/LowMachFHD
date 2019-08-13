module utility_module

  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: invert_multifab, planar_cut

contains

  ! currently invert a cell-centered multifab
  ! only works for cell-centered at the moment
  subroutine invert_multifab(la,phi,comp,ncomp,nghost,numerator)

    type(layout)   , intent(in   ) :: la
    type(multifab) , intent(inout) :: phi
    integer        , intent(in   ) :: comp, ncomp, nghost
    real(kind=dp_t), intent(in), optional :: numerator

    ! local
    type(bl_prof_timer), save :: bpt

    type(multifab) :: temp

    integer :: i
    real(kind=dp_t) :: num

    call build(bpt,"invert_multifab")
    
    if(present(numerator)) then
      num=numerator
    else
      num=1.0d0
    end if    
    
    call multifab_build(temp,la,1,nghost) ! Temporary multifab
    do i=comp,comp+ncomp-1
       call setval(temp,numerator,all=.true.)
       call multifab_div_div_c(temp,1,phi,i,1,nghost)
       call multifab_copy_c(phi,i,temp,1,1,nghost)
    end do

    call multifab_destroy(temp)
    
    call destroy(bpt)

  end subroutine invert_multifab

  ! takes a planar cut of a staggered multifab in the dir direction
  ! at the lo (index 0) side of the domain
  subroutine planar_cut(mla,mf,comp,dir)

    use probin_common_module, only: n_cells

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: mf(:)
    integer        , intent(in   ) :: comp, dir

    integer :: dm

    dm = mla%dim

    if (dm .eq. 2) then
       call planar_cut_2d(mla,mf,comp,dir)
    else if (dm .eq. 3) then
       call bl_error("planar_cut_3d not written yet")
    end if

    contains

      subroutine planar_cut_2d(mla,mf,comp,dir)

        type(ml_layout), intent(in   ) :: mla
        type(multifab) , intent(in   ) :: mf(:)
        integer        , intent(in   ) :: comp
        integer        , intent(in   ) :: dir

        ! local
        integer :: lo(mla%dim),hi(mla%dim)

        real(kind=dp_t), allocatable :: cut(:)
        real(kind=dp_t), allocatable :: cut_proc(:)
        real(kind=dp_t), pointer :: mp(:,:,:,:)

        integer :: i,ng

        if (dir .eq. 1) then
           allocate(cut     (0:n_cells(2)-1))
           allocate(cut_proc(0:n_cells(2)-1))
        else if (dir .eq. 2) then
           allocate(cut     (0:n_cells(1)-1))
           allocate(cut_proc(0:n_cells(1)-1))
        end if

        cut(:) = 0.d0
        cut_proc(:) = 0.d0

        ng = mf(1)%ng

        do i=1,nfabs(mf(1))
           mp => dataptr(mf(1),i)
           lo = lwb(get_box(mf(1),i))
           hi = upb(get_box(mf(1),i))
           call planar_cut_fill_2d(lo,hi,mp(:,:,1,comp),ng,dir,cut_proc)
        end do

        call parallel_reduce(cut,cut_proc,MPI_SUM)

        do i=0,size(cut)-1
           if (parallel_IOProcessor()) then
              print*,i,cut(i)
           end if
        end do

      end subroutine planar_cut_2d

      subroutine planar_cut_fill_2d(lo,hi,mf,ng,dir,cut)

        integer        , intent(in   ) :: lo(:),hi(:),ng,dir
        real(kind=dp_t), intent(in   ) :: mf(lo(1)-ng:,lo(2)-ng:)
        real(kind=dp_t), intent(inout) :: cut(0:)

        integer :: i,j

        if (dir .eq. 1) then
           if (lo(1) .eq. 0) then
              do j=lo(2),hi(2)
                 cut(j) = mf(0,j)
              end do
           end if
        else if (dir .eq. 2) then
           if (lo(2) .eq. 0) then
              do i=lo(1),hi(1)
                 cut(i) = mf(i,0)
              end do
           end if
        end if

      end subroutine planar_cut_fill_2d

    end subroutine planar_cut

end module utility_module
