module utility_module

  use multifab_module
  use ml_layout_module
  use user_analysis
  use parallel

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
  subroutine planar_cut(mla,mf,dir,id)

    use probin_common_module, only: n_cells

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: mf
    integer        , intent(in   ) :: dir
    integer, intent(in), optional  :: id

    integer :: dm

    dm = mla%dim

    if (dm .eq. 2) then
       call planar_cut_2d(mla,mf,dir,id)
    else if (dm .eq. 3) then
       call planar_cut_3d(mla,mf,dir,id)
    end if

    contains

      subroutine planar_cut_2d(mla,mf,dir,id)

        type(ml_layout), intent(in   ) :: mla
        type(multifab) , intent(in   ) :: mf
        integer        , intent(in   ) :: dir
        integer, intent(in), optional  :: id

        ! local
        integer :: lo(mla%dim),hi(mla%dim)

        ! Use the same rank for 2D and 3D to simplify (dir1,dir2,comp)
        real(kind=dp_t), allocatable :: cut     (:,:,:)
        real(kind=dp_t), allocatable :: cut_proc(:,:,:)

        
        real(kind=dp_t), pointer :: mp(:,:,:,:)

        integer :: i,j,ng
        integer :: ncell1,ncell2,ncomp


        if (dir .eq. 1) then
           ncell1 = n_cells(2)
        else if (dir .eq. 2) then
           ncell1 = n_cells(1)
        end if
        ncell2 = 1
        ncomp = mf%nc

        if (dir .eq. 1) then
           allocate(cut     (0:ncell1-1,0:ncell2-1,1:ncomp))
           allocate(cut_proc(0:ncell1-1,0:ncell2-1,1:ncomp))
        else if (dir .eq. 2) then
           allocate(cut     (0:ncell1-1,0:ncell2-1,1:ncomp))
           allocate(cut_proc(0:ncell1-1,0:ncell2-1,1:ncomp))
        end if

        cut = 0.d0
        cut_proc = 0.d0

        ng = mf%ng

        do i=1,nfabs(mf)
           mp => dataptr(mf,i)
           lo = lwb(get_box(mf,i))
           hi = upb(get_box(mf,i))
           call planar_cut_fill_2d(lo,hi,mp(:,:,1,:),ng,dir,ncomp,cut_proc)
        end do

        do j=1,ncomp
        do i=0,ncell2-1
           call parallel_reduce(cut(:,i,j),cut_proc(:,i,j),MPI_SUM)
        end do
        end do
        
        call analyze_planar_cut(cut,id)

      end subroutine planar_cut_2d

      subroutine planar_cut_fill_2d(lo,hi,mf,ng,dir,ncomp,cut)

        integer        , intent(in   ) :: lo(:),hi(:),ng,dir,ncomp
        real(kind=dp_t), intent(in   ) :: mf(lo(1)-ng:,lo(2)-ng:,1:)
        real(kind=dp_t), intent(inout) :: cut(0:,0:,1:)

        integer :: i,j,comp

        do comp=1,ncomp
        
           if (dir .eq. 1) then
              if (lo(1) .eq. 0) then
                 do j=lo(2),hi(2)
                    cut(j,0,comp) = mf(0,j,comp)
                 end do
              end if
           else if (dir .eq. 2) then
              if (lo(2) .eq. 0) then
                 do i=lo(1),hi(1)
                    cut(i,0,comp) = mf(i,0,comp)
                 end do
              end if
           end if

        end do
           
      end subroutine planar_cut_fill_2d

      subroutine planar_cut_3d(mla,mf,dir,id)

        type(ml_layout), intent(in   ) :: mla
        type(multifab) , intent(in   ) :: mf
        integer        , intent(in   ) :: dir
        integer, intent(in), optional  :: id

        ! local
        integer :: lo(mla%dim),hi(mla%dim)

        ! Use the same rank for 2D and 3D to simplify (dir1,dir2,comp)
        real(kind=dp_t), allocatable :: cut     (:,:,:)
        real(kind=dp_t), allocatable :: cut_proc(:,:,:)
        
        real(kind=dp_t), pointer :: mp(:,:,:,:)

        integer :: i,j,ng
        integer :: ncell1,ncell2,ncomp

        if (dir .eq. 1) then
           ncell1 = n_cells(2)
           ncell2 = n_cells(3)
        else if (dir .eq. 2) then
           ncell1 = n_cells(1)
           ncell2 = n_cells(3)
        else if (dir .eq. 3) then
           ncell1 = n_cells(1)
           ncell2 = n_cells(2)
        end if
        ncomp = mf%nc

        allocate(cut     (0:ncell1-1,0:ncell2-1,1:ncomp))
        allocate(cut_proc(0:ncell1-1,0:ncell2-1,1:ncomp))

        cut = 0.d0
        cut_proc = 0.d0

        ng = mf%ng

        do i=1,nfabs(mf)
           mp => dataptr(mf,i)
           lo = lwb(get_box(mf,i))
           hi = upb(get_box(mf,i))
           call planar_cut_fill_3d(lo,hi,mp(:,:,:,:),ng,dir,ncomp,cut_proc)
        end do

        do j=1,ncomp
        do i=0,ncell2-1
           call parallel_reduce(cut(:,i,j),cut_proc(:,i,j),MPI_SUM)
        end do
        end do

        call analyze_planar_cut(cut,id)
        
      end subroutine planar_cut_3d

      subroutine planar_cut_fill_3d(lo,hi,mf,ng,dir,ncomp,cut)

        integer        , intent(in   ) :: lo(:),hi(:),ng,dir,ncomp
        real(kind=dp_t), intent(in   ) :: mf(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,1:)
        real(kind=dp_t), intent(inout) :: cut(0:,0:,1:)

        integer :: i,j,k,comp

        do comp=1,ncomp

           if (dir .eq. 1) then
              if (lo(1) .eq. 0) then
                 do k=lo(3),hi(3)
                    do j=lo(2),hi(2)
                       cut(j,k,comp) = mf(0,j,k,comp)
                    end do
                 end do
              end if
           else if (dir .eq. 2) then
              if (lo(2) .eq. 0) then
                 do k=lo(3),hi(3)
                    do i=lo(1),hi(1)
                       cut(i,k,comp) = mf(i,0,k,comp)
                    end do
                 end do
              end if
           else if (dir .eq. 3) then
              if (lo(3) .eq. 0) then
                 do j=lo(2),hi(2)
                    do i=lo(1),hi(1)
                       cut(i,j,comp) = mf(i,j,0,comp)
                    end do
                 end do
              end if
           end if

        end do

      end subroutine planar_cut_fill_3d

    end subroutine planar_cut

end module utility_module
