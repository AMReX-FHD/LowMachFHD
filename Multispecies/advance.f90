module advance_module

  use multifab_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use div_and_grad_module
  use diffusive_flux_module

  implicit none

  private

  public :: advance

contains
  
  subroutine advance(rho,dx,dt,the_bc_tower)

    type(multifab) , intent(inout) :: rho(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(in   ) :: dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer i, dm, nsp, n, nlevs

    ! an array of multifabs; one for each direction
    ! Amit : Fluxes don't need to remember levels
    type(multifab) :: flux(rho(1)%dim) !, nlevs) 

    nsp=4            ! number of species
    dm = rho(1)%dim  ! dimensionality

    nlevs = size(rho,1)    

    ! build the flux(:) multifabs
    do n=1,nlevs
       do i=1,dm
          ! flux(i) has one component, zero ghost cells, and is nodal in direction i
          call multifab_build_edge(flux(i),rho(1)%la,nsp,0,i)
       end do

       ! compute the face-centered gradients in each direction
       call diffusive_flux(rho(n),flux,dx(n,1),the_bc_tower)
    
       ! update rho using forward Euler discretization
       call update_rho(rho(n),flux,dx(n,1),dt,the_bc_tower)

       ! destroy the multifab to prevent leakage in memory
       do i=1,dm
          call multifab_destroy(flux(i))
       end do
    end do

  end subroutine advance

  subroutine update_rho(rho,flux,dx,dt,the_bc_tower)

    type(multifab) , intent(inout) :: rho
    type(multifab) , intent(in   ) :: flux(:)
    real(kind=dp_t), intent(in   ) :: dx
    real(kind=dp_t), intent(in   ) :: dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: lo(rho%dim), hi(rho%dim)
    integer :: dm, ng_p, ng_f, i, nsp

    real(kind=dp_t), pointer ::  pp(:,:,:,:)
    real(kind=dp_t), pointer :: fxp(:,:,:,:)
    real(kind=dp_t), pointer :: fyp(:,:,:,:)
    real(kind=dp_t), pointer :: fzp(:,:,:,:)

    nsp  = 4
    dm   = rho%dim
    ng_p = rho%ng
    ng_f = flux(1)%ng

    ! Donev: Use compute_divu to compute div(flux) and store it
    ! And then, do update:
    ! rho = rho + dt*div(flux)
    ! Using BoxLib routines
    ! call multifab_mult_mult_s_c(fluxdiv,1,dt,nspecies,0) ! Multiply by a scalar
    ! call multifab_plus_plus_c(rho,1,fluxdiv,1,nspecies,0)
    ! Note: Start with component 1 in rho and fluxdiv, copy nspecies components
    ! or, I think this also works
    ! call multifab_plus_plus(rho,fluxdiv)

    do i=1,nfabs(rho)
       pp  => dataptr(rho,i)
       fxp => dataptr(flux(1),i)
       fyp => dataptr(flux(2),i)
       lo = lwb(get_box(rho,i))
       hi = upb(get_box(rho,i))
       select case(dm)
       case (2)
          call update_rho_2d(pp(:,:,1,nsp), ng_p, &
                             fxp(:,:,1,nsp),  fyp(:,:,1,nsp), ng_f, &
                             lo, hi, dx, dt)
       case (3)
          fzp => dataptr(flux(3),i)
          call update_rho_3d(pp(:,:,:,nsp), ng_p, &
                             fxp(:,:,:,nsp),  fyp(:,:,:,nsp), fzp(:,:,:,nsp), ng_f, &
                             lo, hi, dx, dt)
       end select
    end do
    
    ! fill ghost cells for two adjacent grids at the same level
    ! this includes periodic domain boundary ghost cells
    call multifab_fill_boundary(rho)

    ! fill non-periodic domain boundary ghost cells
    call multifab_physbc(rho,1,1,nsp,the_bc_tower%bc_tower_array(1))  
    ! Amit: not sure about the first 1,1 entry. 

  end subroutine update_rho

  subroutine update_rho_2d(rho, ng_p, fluxx, fluxy, ng_f, lo, hi, dx, dt)

    integer          :: lo(2), hi(2), ng_p, ng_f
    double precision ::   rho(lo(1)-ng_p:,lo(2)-ng_p:)
    double precision :: fluxx(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision :: fluxy(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision :: dx, dt

    ! local variables
    integer i,j

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          rho(i,j) = rho(i,j) + dt * &
               ( fluxx(i+1,j)-fluxx(i,j) + fluxy(i,j+1)-fluxy(i,j) ) / dx
       end do
    end do

  end subroutine update_rho_2d

  subroutine update_rho_3d(rho, ng_p, fluxx, fluxy, fluxz, ng_f, lo, hi, dx, dt)

    integer          :: lo(3), hi(3), ng_p, ng_f
    double precision ::   rho(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    double precision :: fluxx(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: fluxy(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: fluxz(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: dx, dt

    ! local variables
    integer i,j,k

    !$omp parallel do private(i,j,k)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rho(i,j,k) = rho(i,j,k) + dt * &
                  ( fluxx(i+1,j,k)-fluxx(i,j,k) &
                   +fluxy(i,j+1,k)-fluxy(i,j,k) &
                   +fluxz(i,j,k+1)-fluxz(i,j,k) ) / dx
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine update_rho_3d

end module advance_module
