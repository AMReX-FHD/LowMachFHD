module advance_module

  use multifab_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use div_and_grad_module
  use diffusive_flux_module
  use probin_multispecies_module
  use ml_layout_module

  implicit none

  private

  public :: advance

contains
  
  subroutine advance(rho,dx,dt,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rho(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(in   ) :: dt
    type(bc_level) , intent(in   ) :: the_bc_level

    ! local variables
    integer i, dm, n, nlevs

    ! an array of multifabs; one for each direction
    type(multifab) :: flux(mla%nlevel,mla%dim)
    type(multifab) :: fluxdiv(mla%nlevel)
 
    dm = mla%dim  ! dimensionality
    nlevs = mla%nlevel  ! number of levels 

    ! build the flux(:) multifabs
    do n=1,nlevs
       call multifab_build(fluxdiv(n),mla%la(n),nspecies,0)
       do i=1,dm
          ! flux(i) has one component, zero ghost cells, and is nodal in direction i
          call multifab_build_edge(flux(n,i),mla%la(n),nspecies,0,i)
       end do
    end do   
    
    ! Donev: Fix bc_tower vs bc_level
       
    ! compute the face-centered flux in each direction
    call diffusive_flux(rho,flux,dx,the_bc_tower)

    ! compute div of the flux 
    call compute_div(rho(n),flux,div,dx(n,1:dm),1,1,nspecies)

    ! update rho using forward Euler discretization
    !call update_rho(rho(n),flux,dx(n,1),dt,the_bc_tower)
    call update_rho(rho(n),div,dt,the_bc_tower)

    ! destroy the multifab to prevent leakage in memory
    do i=1,dm
       call multifab_destroy(flux(i))
    end do

  end subroutine advance

  subroutine update_rho(rho,div,dt,the_bc_tower)

    type(multifab) , intent(inout) :: rho
    type(multifab) , intent(in   ) :: div
    real(kind=dp_t), intent(in   ) :: dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: lo(rho%dim), hi(rho%dim)
    integer :: dm, ng_p, ng_f, i, nspecies

    real(kind=dp_t), pointer ::   pp(:,:,:,:)
    real(kind=dp_t), pointer :: divp(:,:,:,:)

    nspecies  = 4 ! passing nspecies to this subroutine? 
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
    ! Amit: This we have to do within div_and_grad routine?

    do i=1,nfabs(rho)
       pp  => dataptr(rho,i)
       divp => dataptr(div,i)
       lo = lwb(get_box(rho,i))
       hi = upb(get_box(rho,i))
       select case(dm)
       case (2)
          call update_rho_2d(pp(:,:,1,nspecies), ng_p, &
                             divp(:,:,1,nspecies), lo, hi, dt)
       case (3)
          fzp => dataptr(flux(3),i)
          call update_rho_3d(pp(:,:,:,nspecies), ng_p, &
                             divp(:,:,:,nspecies), lo, hi, dt)
       end select
    end do
    
    ! fill ghost cells for two adjacent grids at the same level
    ! this includes periodic domain boundary ghost cells
    call multifab_fill_boundary(rho)

    ! fill non-periodic domain boundary ghost cells
    call multifab_physbc(rho,1,1,nspecies,the_bc_tower%bc_tower_array(1))  

  end subroutine update_rho

  subroutine update_rho_2d(rho, ng_p, divp, lo, hi, dt)

    integer          :: lo(2), hi(2), ng_p
    real(kind=dp_t) ::  rho(lo(1)-ng_p:,lo(2)-ng_p:)
    real(kind=dp_t) :: divp(lo(1)-ng_p:,lo(2)-ng_p:)
    real(kind=dp_t) :: dt

    rho(:,:) = rho(:,:) + dt * divp(:,:) 

  end subroutine update_rho_2d

  subroutine update_rho_3d(rho, ng_p, divp, lo, hi, dt)

    integer          :: lo(3), hi(3), ng_p
    real(kind=dp_t) ::  rho(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real(kind=dp_t) :: divp(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real(kind=dp_t) :: dt

    rho(:,:,:) = rho(:,:,:) + dt * divp(:,:,:)

  end subroutine update_rho_3d

end module advance_module
