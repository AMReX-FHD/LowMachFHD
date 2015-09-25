module init_n_module

  use bl_types
  use ml_layout_module
  use multifab_physbc_module
  use define_bc_module
  use bc_module
  use probin_common_module, only: prob_lo, prob_hi, prob_type
  use probin_reactdiff_module, only: nspecies, n_init_in
  
  implicit none

  private

  public :: init_n

  ! prob_type codes for LowMach:
  ! 0=thermodynamic equilibrium, n=n_init_in(1,1:nspecies)
  ! 1=gaussian spreading (order of accuracy testing)
  ! 2=gradient along y, n=n_init_in(1,1:nspecies) on bottom (y=0) and n_init(2,1:nspecies) on top (y=Ly)

contains

  subroutine init_n(mla,n_init,dx,the_bc_tower)

    ! initialize rho_i and umac in the valid region
    ! we first initialize c_i in the valid region
    ! then enforce that sum(c_i)=1 by overwriting the final concentration,
    ! and then use the EOS to compute rho_i

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: n_init(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower ), intent(in   ) :: the_bc_tower
 
    ! local variables
    integer                        :: lo(mla%dim), hi(mla%dim)
    integer                        :: i, dm, n, nlevs, ng_n
    real(kind=dp_t), pointer       :: np(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "init_n")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_n = n_init(1)%ng

    ! looping over boxes 
    do n=1,nlevs
       do i=1,nfabs(n_init(n))
          np => dataptr(n_init(n),i)
          lo = lwb(get_box(n_init(n),i))
          hi = upb(get_box(n_init(n),i))
          select case (dm)
          case (2)
             call init_n_2d(np(:,:,1,:),ng_n,lo,hi,dx(n,:))
          case (3)
             call init_n_3d(np(:,:,:,:),ng_n,lo,hi,dx(n,:))
          end select
       end do
    end do

    call destroy(bpt)

   do n=1,nlevs
      call multifab_fill_boundary(n_init(n))
      call multifab_physbc(n_init(n),1,scal_bc_comp,nspecies, &
                           the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
   end do


  end subroutine init_n

  subroutine init_n_2d(n_init,ng_n,lo,hi,dx)

    integer          :: lo(:), hi(:), ng_n
    real(kind=dp_t)  :: n_init(lo(1)-ng_n:,lo(2)-ng_n:,:)
    real(kind=dp_t)  :: dx(:)
 
    ! local varables
    integer         :: i,j,comp
    real(kind=dp_t) :: x,y,r,cen(2),sum,L(2)

    L(1:2) = prob_hi(1:2)-prob_lo(1:2) ! Domain length
    
    select case (prob_type)

    case(0) 
       !============================================================
       ! Thermodynamic equilibrium
       !============================================================

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             n_init(i,j,1:nspecies) = n_init_in(1,1:nspecies)

          end do
       end do

    case(1) 
       !=============================================================
       ! Initializing from a Gaussian
       !=============================================================

       cen(1:2) = 0.6d0*prob_lo(1:2) + 0.4d0*prob_hi(1:2)

       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)

             r = sqrt((x-cen(1))**2 + (y-cen(2))**2)

             n_init(i,j,1:nspecies) = n_init_in(1,1:nspecies)*exp(-100.d0*r**2)

          end do
       end do

    case(2) 
       !=========================================================
       ! Initializing with constant gradient 
       !=========================================================

       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2) 
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1) 

             ! linear gradient in rho
             n_init(i,j,1:nspecies) = n_init_in(1,1:nspecies) + & 
                  (n_init_in(2,1:nspecies) - n_init_in(1,1:nspecies))*(y-prob_lo(2))/L(2)

          end do
       end do

    case default

       call bl_error("init_n_2d: prob_type not supported")

    end select

  end subroutine init_n_2d

  subroutine init_n_3d(n_init,ng_n,lo,hi,dx)

    integer          :: lo(:), hi(:), ng_n
    real(kind=dp_t)  :: n_init(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:,:)
    real(kind=dp_t)  :: dx(:)
 
    ! local varables
    integer         :: i,j,k,comp
    real(kind=dp_t) :: x,y,z,r,cen(3),sum,L(3)

    L(1:3) = prob_hi(1:3)-prob_lo(1:3) ! Domain length
    
    select case (prob_type)

    case(0) 
       !================================================================================
       ! Thermodynamic equilibrium
       !================================================================================

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)

                n_init(i,j,k,1:nspecies) = n_init_in(1,1:nspecies)

             end do
          end do
       end do

    case(1) 
       !================================================================================
       ! Initializing from a Gaussian
       !================================================================================
       cen(1:3) = 0.6d0*prob_lo(1:3) + 0.4d0*prob_hi(1:3)

       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             do i=lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)

                r = sqrt((x-cen(1))**2 + (y-cen(2))**2 + (z-cen(3))**2)

                n_init(i,j,k,1:nspecies) = n_init_in(1,1:nspecies)*exp(-100.d0*r**2)

             end do
          end do
       end do

    case(2) 
       !========================================================
       ! Initializing with constant gradient
       !========================================================

       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3) 
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             do i=lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)

                n_init(i,j,k,1:nspecies) = n_init_in(1,1:nspecies) + &
                     (n_init_in(2,1:nspecies) - n_init_in(1,1:nspecies))*(y-prob_lo(2))/L(2)

             end do
          end do
       end do

    case default

       call bl_error("init_n_3d: prob_type not supported")

    end select

  end subroutine init_n_3d

end module init_n_module
