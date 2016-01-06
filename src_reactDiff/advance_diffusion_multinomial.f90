module advance_diffusion_multinomial_module

  use ml_layout_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use probin_reactdiff_module, only: nspecies, D_Fick

  implicit none

  private

  public :: advance_diffusion_multinomial

contains

  ! advances n_old to n_new using multinomial diffusion
  subroutine advance_diffusion_multinomial(mla,n_old,n_new,dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: n_old(:)
    type(multifab) , intent(inout) :: n_new(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    type(multifab) :: diff_coef_face(mla%nlevel,mla%dim)

    integer :: n,nlevs,i,dm,spec

    type(bl_prof_timer),save :: bpt

    nlevs = mla%nlevel
    dm = mla%dim

    ! do not do diffusion if only one cell (well-mixed system)
    ! there is no restriction on the number of cells
    ! but we can shortcut the single cell case anyway for simplicity
    if((multifab_volume(n_old(1))/nspecies)<=1) then
       do n=1,nlevs
          ! make sure n_new contains the new state
          call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,n_new(n)%ng)
       end do
    end if
    
    call build(bpt,"advance_diffusion")
    
    do n=1,nlevs
       do i=1,dm
          call multifab_build_edge(diff_coef_face(n,i),mla%la(n),nspecies,0,i)
       end do
    end do

    ! compute the diffusion coefficients (for now just setting each to a different constant)
    ! If one wants a space-dependent D or state-dependent D see multispecies code as example
    ! We have a routine average_cc_to_face there that is meant to compute face-averaged values
    do n=1,nlevs
       do i=1,dm
          do spec=1,nspecies
             call multifab_setval_c(diff_coef_face(n,i), D_Fick(spec),spec,1,all=.true.)
          end do
       end do
    end do


    ! set new state to zero everywhere, including ghost cells
    do n=1,nlevs
       call multifab_setval(n_new(n),0.d0,all=.true.)
    end do    

    ! copy old state into new in valid region only
    do n=1,nlevs
       call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
    end do    

    ! update with multinomial diffusion, each grid in isolation



    ! call sum_boundary to deal with grid boundaries
    do n=1,nlevs
       call multifab_sum_boundary(n_new(n),1)
    end do

    ! properly fill n_new ghost cells
    do n=1,nlevs
       call multifab_fill_boundary(n_new(n))
       call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                            the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
    end do

    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(diff_coef_face(n,i))
       end do
    end do

    call destroy(bpt)

  end subroutine advance_diffusion_multinomial

  subroutine multinomial_diffusion_update(mla,n_new,dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: n_new(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: n,nlevs,i,dm,ng_n

    integer      ::  lo(mla%dim),  hi(mla%dim)
    
    real(kind=dp_t), allocatable :: private_np(:,:,:,:)
    real(kind=dp_t), pointer :: np(:,:,:,:)  ! "n" for n_new

    nlevs = mla%nlevel
    dm = mla%dim

    ng_n = n_new(1)%ng

    ! cannot use OpenMP with tiling since each cell is responsible for updating
    ! cells possibly outside of its file.  OpenMP could be added at the k loop level
    ! with reduction tricks
    do n=1,nlevs
       do i=1,nfabs(n_new(n))
        np => dataptr(n_new(n),i)
        lo = lwb(get_box(n_new(n),i))
        hi = upb(get_box(n_new(n),i))
        select case (dm)
        case (2)
          call multinomial_diffusion_update_2d(np(:,:,1,:),ng_n,lo,hi,dx(n,:),dt)
        case (3)
           call bl_error("multinomial_diffusion_update_3d not written yet")
        end select
      end do
    end do

  end subroutine multinomial_diffusion_update

  subroutine multinomial_diffusion_update_2d(n_new,ng_n,lo,hi,dx,dt)

    integer        , intent(in   ) :: lo(:),hi(:),ng_n
    real(kind=dp_t), intent(inout) :: n_new(lo(1)-ng_n:,lo(2)-ng_n:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt

    ! local
    integer :: i,j

    real(kind=dp_t) :: cell_update(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,nspecies)

    cell_update = 0.d0

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

       end do
    end do

    ! increment n_new
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1

       end do
    end do

  end subroutine multinomial_diffusion_update_2d

end module advance_diffusion_multinomial_module
