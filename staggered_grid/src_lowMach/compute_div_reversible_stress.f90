module reversible_stress_module

 
  use ml_layout_module
  use convert_stag_module
  use div_and_grad_module
  use define_bc_module
  use bc_module
  use mass_flux_utilities_module
  use matvec_mul_module
  use compute_mixture_properties_module
  use multifab_physbc_module
  use zero_edgeval_module
  use fill_rho_ghost_cells_module
  use convert_rhoc_to_c_module
  use probin_common_module, only: molmass, k_B, rhobar, nspecies, prob_lo, prob_hi
  use probin_charged_module, only: charge_per_mass, dpdt_factor, &
                                   dielectric_const, dielectric_type, &
                                   E_ext_type, E_ext_value, zero_eps_on_wall_type, &
                                   zero_eps_on_wall_left_end, zero_eps_on_wall_right_start

  implicit none

  private

  public :: compute_div_reversible_stress
  
contains

  subroutine compute_div_reversible_stress(mla,div_reversible_stress,rhotot,rho,dx,the_bc_tower)


!!!need to add concentration
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: div_reversible_stress(:,:)
    type(multifab) , intent(in   ) :: rhotot(:)
    type(multifab) , intent(in   ) :: rho   (:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: i,n,dm,nlevs
    integer :: ng_1,ng_2
    integer :: lo(mla%dim),hi(mla%dim)

    type(multifab) :: conc(mla%nlevel)

    ! pointers into multifabs
    real(kind=dp_t), pointer :: dp1x(:,:,:,:)
    real(kind=dp_t), pointer :: dp1y(:,:,:,:)
    real(kind=dp_t), pointer :: dp1z(:,:,:,:)

    real(kind=dp_t), pointer :: cp(:,:,:,:)

    logical :: use_qE_Lorentz
    
    dm = mla%dim
    nlevs = mla%nlevel

    do n=1,nlevs
       call multifab_build(conc(n),mla%la(n),nspecies,2)
    end do

    ! rho to conc - VALID REGION ONLY
    call convert_rhoc_to_c(mla,rho,rhotot,conc,.true.)

    ! fill conc ghost cells
    call fill_c_ghost_cells(mla,conc,dx,the_bc_tower)
    
    ! spatially-varying permittivity
    ! add -(1/2) E^2 grad(eps)
    
    ng_1 = div_reversible_stress(1,1)%ng
    ng_2 = conc(1)%ng

    do n=1,nlevs
       do i=1,nfabs(div_reversible_stress(n,1))
          dp1x => dataptr(div_reversible_stress(n,1),i)
          dp1y => dataptr(div_reversible_stress(n,2),i)
          cp => dataptr(conc(n),i)
          lo = lwb(get_box(div_reversible_stress(n,1),i))
          hi = upb(get_box(div_reversible_stress(n,1),i))
          select case (dm)
          case (2)
             call compute_div_reversible_stress_2d(dp1x(:,:,1,1),dp1y(:,:,1,1),ng_1, &
                                                   cp(:,:,1,1),ng_2,lo,hi,dx(n,:))
          case (3)
             dp1z => dataptr(div_reversible_stress(n,3),i)
             call compute_div_reversible_stress_3d(dp1x(:,:,:,1),dp1y(:,:,:,1),dp1z(:,:,:,1),ng_1, &
                                                   cp(:,:,:,1),ng_2,lo,hi,dx(n,:))
          end select
       end do
    end do


    ! set force on walls to be zero since normal velocity is zero
    do n=1,nlevs
       call zero_edgeval_walls(div_reversible_stress(n,:),1,1,the_bc_tower%bc_tower_array(n))
    end do

    do n=1,nlevs
       call multifab_destroy(conc(n))
       
    end do

  contains

    subroutine compute_div_reversible_stress_2d(forcex,forcey,ng_1,c,ng_2,lo,hi,dx)

      integer         :: lo(:),hi(:),ng_1,ng_2
      real(kind=dp_t) :: forcex(lo(1)-ng_1:,lo(2)-ng_1:)
      real(kind=dp_t) :: forcey(lo(1)-ng_1:,lo(2)-ng_1:)
      real(kind=dp_t) ::      c(lo(1)-ng_2:,lo(2)-ng_2:)
      real(kind=dp_t) :: dx(:)

      ! local variables
      integer :: i,j 
      real(kind=dp_t) :: kc
      real(kind=dp_t) :: cx_local_plus, cy_local_plus, cx_local_minus, cy_local_minus

      real(kind=dp_t) , allocatable :: node_grad_cx(:,:)
      real(kind=dp_t) , allocatable :: node_grad_cy(:,:)

      allocate(node_grad_cx(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2))
      allocate(node_grad_cy(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2))

      kc = 1.d-4 !make parameter

      !grad x & y
      do j=lo(2)-1,hi(2)+2
       do i=lo(1)-1,hi(1)+2 !!does the plus one matter here?

              node_grad_cx(i,j) = (c(i,j)-c(i-1,j)+c(i,j-1)-c(i-1,j-1))/(2*dx(1))
              node_grad_cy(i,j) = (c(i,j)-c(i,j-1)+c(i-1,j)-c(i-1,j-1))/(2*dx(2))

       end do
      end do 

      do j=lo(2),hi(2)
      do i=lo(1),hi(1)+1
        cx_local_plus = 0.25*(node_grad_cx(i,j)+node_grad_cx(i,j+1)+node_grad_cx(i+1,j+1)+node_grad_cx(i+1,j))
        cy_local_plus = 0.25*(node_grad_cy(i,j)+node_grad_cy(i,j+1)+node_grad_cy(i+1,j+1)+node_grad_cy(i+1,j))
        cx_local_minus = 0.25*(node_grad_cx(i,j)+node_grad_cx(i,j+1)+node_grad_cx(i-1,j+1)+node_grad_cx(i-1,j))
        cy_local_minus = 0.25*(node_grad_cy(i,j)+node_grad_cy(i,j+1)+node_grad_cy(i-1,j+1)+node_grad_cy(i-1,j))

!also have not added k_c anywhere--need to fix
        forcex(i,j) = -(node_grad_cx(i,j+1)*node_grad_cy(i,j+1) - node_grad_cx(i,j)*node_grad_cy(i,j))/dx(2) &
           +(0.5*(cy_local_plus*cy_local_plus-cx_local_plus*cx_local_plus) &
           -0.5*(cy_local_minus*cy_local_minus-cx_local_minus*cx_local_minus))/(dx(1)) 
        forcex(i,j) = kc*forcex(i,j)
      end do
      end do


      do j=lo(2),hi(2)+1
      do i=lo(1),hi(1)

        cx_local_plus = 0.25*(node_grad_cx(i,j)+node_grad_cx(i,j+1)+node_grad_cx(i+1,j+1)+node_grad_cx(i+1,j))
        cy_local_plus = 0.25*(node_grad_cy(i,j)+node_grad_cy(i,j+1)+node_grad_cy(i+1,j+1)+node_grad_cy(i+1,j))
        cx_local_minus = 0.25*(node_grad_cx(i,j)+node_grad_cx(i,j-1)+node_grad_cx(i+1,j-1)+node_grad_cx(i+1,j))
        cy_local_minus = 0.25*(node_grad_cy(i,j)+node_grad_cy(i,j-1)+node_grad_cy(i+1,j-1)+node_grad_cy(i+1,j))

        forcey(i,j) = -(node_grad_cx(i+1,j)*node_grad_cy(i+1,j) - node_grad_cx(i,j)*node_grad_cy(i,j))/dx(1) &
           +(0.5*(cy_local_plus*cy_local_plus-cx_local_plus*cx_local_plus) &
           -0.5*(cy_local_minus*cy_local_minus-cx_local_minus*cx_local_minus))/(dx(2)) 
        forcey(i,j) = kc*forcey(i,j)
      end do
      end do

    end subroutine compute_div_reversible_stress_2d

    subroutine compute_div_reversible_stress_3d(forcex,forcey,forcez,ng_1,c,ng_2,lo,hi,dx)

      integer         :: lo(:),hi(:),ng_1,ng_2
      real(kind=dp_t) :: forcex(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
      real(kind=dp_t) :: forcey(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
      real(kind=dp_t) :: forcez(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
      real(kind=dp_t) ::      c(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_2:)
      real(kind=dp_t) :: dx(:)

      ! local variables
      integer :: i,j,k

      real(kind=dp_t) :: kc
      real(kind=dp_t) :: cx_local_plus, cy_local_plus, cz_local_plus, cx_local_minus, cy_local_minus, cz_local_minus

      real(kind=dp_t),  allocatable :: node_grad_cx(:,:,:)
      real(kind=dp_t),  allocatable :: node_grad_cy(:,:,:)
      real(kind=dp_t),  allocatable :: node_grad_cz(:,:,:)

      allocate(node_grad_cx(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2,lo(3)-1:hi(3)+2))
      allocate(node_grad_cy(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2,lo(3)-1:hi(3)+2))
      allocate(node_grad_cz(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2,lo(3)-1:hi(3)+2))
      
      kc = 1.d-4 !make parameter

      !grad x & y & z
      do k=lo(3)-1,hi(3)+2
       do j=lo(2)-1,hi(2)+2
        do i=lo(1)-1,hi(1)+2 !! again, is this correct placement for +1?

              node_grad_cx(i,j,k) = (c(i,j,k)-c(i-1,j,k)+c(i,j-1,k)-c(i-1,j-1,k) &
		+c(i,j,k-1)-c(i-1,j,k-1)+c(i,j-1,k-1)-c(i-1,j-1,k-1))/(4*dx(1))
              node_grad_cy(i,j,k) = (c(i,j,k)-c(i,j-1,k)+c(i-1,j,k)-c(i-1,j-1,k) &
		+c(i,j,k-1)-c(i,j-1,k-1)+c(i-1,j,k-1)-c(i-1,j-1,k-1))/(4*dx(2))
              node_grad_cz(i,j,k) = (c(i,j,k)-c(i,j,k-1)+c(i-1,j,k)-c(i-1,j,k-1) &
		+c(i,j-1,k)-c(i,j-1,k-1)+c(i-1,j-1,k)-c(i-1,j-1,k-1))/(4*dx(3))

        end do
       end do 
      end do

      do k=lo(3),hi(3)
       do j=lo(2),hi(2)
        do i=lo(1),hi(1)+1
         cx_local_plus = (.125)*(node_grad_cx(i,j,k)+node_grad_cx(i,j+1,k)+node_grad_cx(i+1,j+1,k) &
		+node_grad_cx(i+1,j,k)+node_grad_cx(i+1,j+1,k+1)+node_grad_cx(i,j,k+1) &
              +node_grad_cx(i,j+1,k+1)+node_grad_cx(i+1,j,k+1))
         cy_local_plus = (.125)*(node_grad_cy(i,j,k)+node_grad_cy(i,j+1,k)+node_grad_cy(i+1,j+1,k) &
              +node_grad_cy(i+1,j,k)+node_grad_cy(i+1,j+1,k+1)+node_grad_cy(i,j,k+1) &
              +node_grad_cy(i,j+1,k+1)+node_grad_cy(i+1,j,k+1))
         cz_local_plus = (.125)*(node_grad_cz(i,j,k)+node_grad_cz(i,j+1,k)+node_grad_cz(i+1,j+1,k) &
              +node_grad_cz(i+1,j,k)+node_grad_cz(i+1,j+1,k+1)+node_grad_cz(i,j,k+1)+ &
              node_grad_cz(i,j+1,k+1)+node_grad_cz(i+1,j,k+1))
         cx_local_minus = (.125)*(node_grad_cx(i,j,k)+node_grad_cx(i,j+1,k)+node_grad_cx(i-1,j+1,k)+&
                 node_grad_cx(i-1,j,k)+node_grad_cx(i,j,k+1)+node_grad_cx(i,j+1,k+1)+&
                 node_grad_cx(i-1,j+1,k+1)+node_grad_cx(i-1,j,k+1))
         cy_local_minus = (.125)*(node_grad_cy(i,j,k)+node_grad_cy(i,j+1,k)+node_grad_cy(i-1,j+1,k)+&
                 node_grad_cy(i-1,j,k)+node_grad_cy(i,j,k+1)+node_grad_cy(i,j+1,k+1)+&
                 node_grad_cy(i-1,j+1,k+1)+node_grad_cy(i-1,j,k+1))
         cz_local_minus = (.125)*(node_grad_cz(i,j,k)+node_grad_cz(i,j+1,k)+node_grad_cz(i-1,j+1,k)&
                 +node_grad_cz(i-1,j,k)+node_grad_cz(i,j,k+1)+node_grad_cz(i,j+1,k+1)+&
                 node_grad_cz(i-1,j+1,k+1)+node_grad_cz(i-1,j,k+1))

         forcex(i,j,k) = -(node_grad_cx(i,j+1,k)*node_grad_cy(i,j+1,k) -&
                  node_grad_cx(i,j,k)*node_grad_cy(i,j,k))/dx(2) &
           +(0.5*(cy_local_plus*cy_local_plus+cz_local_plus*cz_local_plus-&
                 cx_local_plus*cx_local_plus) &
           -0.5*(cy_local_minus*cy_local_minus+cz_local_minus*cz_local_minus-&
                 cx_local_minus*cx_local_minus))/(dx(1)) & 
           -(node_grad_cx(i,j,k+1)*node_grad_cz(i,j,k+1) -&
                  node_grad_cx(i,j,k)*node_grad_cz(i,j,k))/dx(3)
         forcex(i,j,k) = kc*forcex(i,j,k)

         end do
        end do
       end do


      do k=lo(3),hi(3)
      do j=lo(2),hi(2)+1
      do i=lo(1),hi(1)

         cx_local_plus = (.125)*(node_grad_cx(i,j,k)+node_grad_cx(i,j+1,k)+&
                 node_grad_cx(i+1,j+1,k)+node_grad_cx(i+1,j,k)+node_grad_cx(i+1,j+1,k+1)+&
                 node_grad_cx(i,j,k+1)+node_grad_cx(i,j+1,k+1)+node_grad_cx(i+1,j,k+1))
         cy_local_plus = (.125)*(node_grad_cy(i,j,k)+node_grad_cy(i,j+1,k)+node_grad_cy(i+1,j+1,k)+&
                 node_grad_cy(i+1,j,k)+node_grad_cy(i+1,j+1,k+1)+&
                 node_grad_cy(i,j,k+1)+node_grad_cy(i,j+1,k+1)+node_grad_cy(i+1,j,k+1))
         cz_local_plus = (.125)*(node_grad_cz(i,j,k)+node_grad_cz(i,j+1,k)+&
                 node_grad_cz(i+1,j+1,k)+node_grad_cz(i+1,j,k)+node_grad_cz(i+1,j+1,k+1)+&
                 node_grad_cz(i,j,k+1)+node_grad_cz(i,j+1,k+1)+node_grad_cz(i+1,j,k+1))

         cx_local_minus = (.125)*(node_grad_cx(i,j,k)+node_grad_cx(i,j-1,k)+&
                 node_grad_cx(i+1,j-1,k)+node_grad_cx(i+1,j,k)+node_grad_cx(i,j,k+1)+&
                 node_grad_cx(i,j-1,k+1)+node_grad_cx(i+1,j-1,k+1)+node_grad_cx(i+1,j,k+1))
         cy_local_minus = (.125)*(node_grad_cy(i,j,k)+node_grad_cy(i,j-1,k)+&
                 node_grad_cy(i+1,j-1,k)+node_grad_cy(i+1,j,k)+node_grad_cy(i,j,k+1)+&
                 node_grad_cy(i,j-1,k+1)+node_grad_cy(i+1,j-1,k+1)+node_grad_cy(i+1,j,k+1))
         cz_local_minus = (.125)*(node_grad_cz(i,j,k)+node_grad_cz(i,j-1,k)+&
                 node_grad_cz(i+1,j-1,k)+node_grad_cz(i+1,j,k)+node_grad_cz(i,j,k+1)+&
                 node_grad_cz(i,j-1,k+1)+node_grad_cz(i+1,j-1,k+1)+node_grad_cz(i+1,j,k+1))

         forcey(i,j,k) = -(node_grad_cx(i+1,j,k)*node_grad_cy(i+1,j,k) - node_grad_cx(i,j,k)*node_grad_cy(i,j,k))/dx(1) &
           +(0.5*(cx_local_plus*cx_local_plus+cz_local_plus*cz_local_plus-cy_local_plus*cy_local_plus) &
           -0.5*(cx_local_minus*cx_local_minus+cz_local_minus*cz_local_minus-cy_local_minus*cy_local_minus))/(dx(2)) & 
           -(node_grad_cy(i,j,k+1)*node_grad_cz(i,j,k+1) - node_grad_cy(i,j,k)*node_grad_cz(i,j,k))/dx(3)
         forcey(i,j,k) = kc*forcey(i,j,k)

      end do
      end do
      end do

      do k=lo(3),hi(3)+1
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)

         cx_local_plus = (.125)*(node_grad_cx(i,j,k)+node_grad_cx(i,j+1,k)+&
                 node_grad_cx(i+1,j+1,k)+node_grad_cx(i+1,j,k)+node_grad_cx(i+1,j+1,k+1)+&
                 node_grad_cx(i,j,k+1)+node_grad_cx(i,j+1,k+1)+node_grad_cx(i+1,j,k+1))
         cy_local_plus = (.125)*(node_grad_cy(i,j,k)+node_grad_cy(i,j+1,k)+node_grad_cy(i+1,j+1,k)+&
                 node_grad_cy(i+1,j,k)+node_grad_cy(i+1,j+1,k+1)+node_grad_cy(i,j,k+1)+&
                 node_grad_cy(i,j+1,k+1)+node_grad_cy(i+1,j,k+1))
         cz_local_plus = (.125)*(node_grad_cz(i,j,k)+node_grad_cz(i,j+1,k)+node_grad_cz(i+1,j+1,k)+&
                 node_grad_cz(i+1,j,k)+node_grad_cz(i+1,j+1,k+1)+node_grad_cz(i,j,k+1)+&
                 node_grad_cz(i,j+1,k+1)+node_grad_cz(i+1,j,k+1))

         cx_local_minus = (.125)*(node_grad_cx(i,j,k)+node_grad_cx(i,j+1,k)+node_grad_cx(i+1,j+1,k)+& 
              node_grad_cx(i+1,j,k)+node_grad_cx(i,j,k-1)+node_grad_cx(i,j+1,k-1)& 
              +node_grad_cx(i+1,j+1,k-1)+node_grad_cx(i+1,j,k-1))
         cy_local_minus = (.125)*(node_grad_cy(i,j,k)+node_grad_cy(i,j+1,k)+node_grad_cy(i+1,j+1,k)& 
              +node_grad_cy(i+1,j,k)+node_grad_cy(i,j,k-1)+node_grad_cy(i,j+1,k-1)& 
              +node_grad_cy(i+1,j+1,k-1)+node_grad_cy(i+1,j,k-1))
         cz_local_minus = (.125)*(node_grad_cz(i,j,k)+node_grad_cz(i,j+1,k)+node_grad_cz(i+1,j+1,k)& 
              +node_grad_cz(i+1,j,k)+node_grad_cz(i,j,k-1)+node_grad_cz(i,j+1,k-1)& 
              +node_grad_cz(i+1,j+1,k-1)+node_grad_cz(i+1,j,k-1))

         forcez(i,j,k) = -(node_grad_cx(i+1,j,k)*node_grad_cz(i+1,j,k) - node_grad_cx(i,j,k)*node_grad_cz(i,j,k))/dx(1) &
           +(0.5*(cx_local_plus*cx_local_plus+cy_local_plus*cy_local_plus-cz_local_plus*cz_local_plus) &
           -0.5*(cx_local_minus*cx_local_minus+cy_local_minus*cy_local_minus-cz_local_minus*cz_local_minus))/(dx(3)) & 
           -(node_grad_cy(i,j+1,k)*node_grad_cz(i,j+1,k) - node_grad_cy(i,j,k)*node_grad_cz(i,j,k))/dx(2)
         forcez(i,j,k) = kc*forcez(i,j,k)

      end do
      end do
      end do

    end subroutine compute_div_reversible_stress_3d

  end subroutine compute_div_reversible_stress

end module reversible_stress_module
