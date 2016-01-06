module multinomial_diffusion_module

  use ml_layout_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use probin_reactdiff_module, only: nspecies, D_Fick

  implicit none

  private

  public :: multinomial_diffusion

contains

  ! advances n_old to n_new using multinomial diffusion
  subroutine multinomial_diffusion(mla,n_old,n_new,diff_coef_face, &
                                           dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: n_old(:)
    type(multifab) , intent(inout) :: n_new(:)
    type(multifab) , intent(in   ) :: diff_coef_face(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: n,nlevs,i,dm,spec

    type(bl_prof_timer),save :: bpt

    nlevs = mla%nlevel
    dm = mla%dim
    
    call build(bpt,"multinomial_diffusion")

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

    call destroy(bpt)

  end subroutine multinomial_diffusion

  subroutine multinomial_diffusion_update(mla,n_new,diff_coef_face,dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: n_new(:)
    type(multifab) , intent(in   ) :: diff_coef_face(:,:)
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

end module multinomial_diffusion_module
