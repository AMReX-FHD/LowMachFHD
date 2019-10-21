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
  use probin_common_module, only: molmass, k_B, rhobar, nspecies, prob_lo, prob_hi
  use probin_charged_module, only: charge_per_mass, dpdt_factor, &
                                   dielectric_const, dielectric_type, &
                                   E_ext_type, E_ext_value, zero_eps_on_wall_type, &
                                   zero_eps_on_wall_left_end, zero_eps_on_wall_right_start

  implicit none

  private

  public :: dot_with_z, dot_with_z_face, compute_charge_coef, &
            enforce_charge_neutrality, implicit_potential_coef, modify_S, &
            compute_permittivity, compute_Lorentz_force, compute_E_ext, &
            zero_eps_on_wall
  
contains

  subroutine compute_div_reversible_stress(mla,Lorentz_force,grad_Epot,permittivity,charge, &
                                   dx,the_bc_tower)


!!!need to add concentration
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: Lorentz_force(:,:)
    type(multifab) , intent(in   ) :: grad_Epot(:,:)
    type(multifab) , intent(in   ) :: permittivity(:)
    type(multifab) , intent(in   ) :: charge(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: i,n,dm,nlevs
    integer :: ng_1,ng_2,ng_3
    integer :: lo(mla%dim),hi(mla%dim)

    type(multifab) :: temp_cc(mla%nlevel)

    ! pointers into multifabs
    real(kind=dp_t), pointer :: dp1x(:,:,:,:)
    real(kind=dp_t), pointer :: dp1y(:,:,:,:)
    real(kind=dp_t), pointer :: dp1z(:,:,:,:)
    real(kind=dp_t), pointer :: dp2x(:,:,:,:)
    real(kind=dp_t), pointer :: dp2y(:,:,:,:)
    real(kind=dp_t), pointer :: dp2z(:,:,:,:)
    real(kind=dp_t), pointer :: dp3(:,:,:,:)

    logical :: use_qE_Lorentz
    
    dm = mla%dim
    nlevs = mla%nlevel

    if (dielectric_type .ne. 0) then

       ! spatially-varying permittivity
       ! add -(1/2) E^2 grad(eps)

       ng_1 = Lorentz_force(1,1)%ng
       ng_2 = grad_Epot(1,1)%ng
       ng_3 = permittivity(1)%ng

       do n=1,nlevs
          do i=1,nfabs(Lorentz_force(n,1))
             dp1x => dataptr(Lorentz_force(n,1),i)
             dp1y => dataptr(Lorentz_force(n,2),i)
             dp2x => dataptr(grad_Epot(n,1),i)
             dp2y => dataptr(grad_Epot(n,2),i)
             dp3  => dataptr(permittivity(n),i)
             lo = lwb(get_box(Lorentz_force(n,1),i))
             hi = upb(get_box(Lorentz_force(n,1),i))
             select case (dm)
             case (2)
                call compute_div_reversible_stress_2d(dp1x(:,:,1,1),dp1y(:,:,1,1),ng_1, &
                                              dp2x(:,:,1,1),dp2y(:,:,1,1),ng_2, &
                                              dp3(:,:,1,1),ng_3,lo,hi,dx(n,:))
             case (3)
                dp1z => dataptr(Lorentz_force(n,3),i)
                dp2z => dataptr(grad_Epot(n,3),i)
                call compute_div_reversible_stress_3d(dp1x(:,:,:,1),dp1y(:,:,:,1),dp1z(:,:,:,1),ng_1, &
                                              dp2x(:,:,:,1),dp2y(:,:,:,1),dp2z(:,:,:,1),ng_2, &
                                              dp3(:,:,:,1),ng_3,lo,hi,dx(n,:))
             end select
          end do
       end do

    end if

    ! set force on walls to be zero since normal velocity is zero
    do n=1,nlevs
       call zero_edgeval_walls(Lorentz_force(n,:),1,1,the_bc_tower%bc_tower_array(n))
    end do

    do n=1,nlevs
       call multifab_destroy(temp_cc(n))
    end do

  contains

    subroutine compute_div_reversible_stress_2d(forcex,forcey,ng_1,Ex,Ey,ng_2,perm,ng_3,lo,hi,dx)

      integer         :: lo(:),hi(:),ng_1,ng_2,ng_3
      real(kind=dp_t) :: forcex(lo(1)-ng_1:,lo(2)-ng_1:)
      real(kind=dp_t) :: forcey(lo(1)-ng_1:,lo(2)-ng_1:)
      real(kind=dp_t) ::     node_grad_cx(lo(1)-ng_2:,lo(2)-ng_2:)
      real(kind=dp_t) ::     node_grad_cy(lo(1)-ng_2:,lo(2)-ng_2:)
      real(kind=dp_t) ::     Ex(lo(1)-ng_2:,lo(2)-ng_2:)
      real(kind=dp_t) ::     Ey(lo(1)-ng_2:,lo(2)-ng_2:)
      real(kind=dp_t) ::   perm(lo(1)-ng_3:,lo(2)-ng_3:)
      real(kind=dp_t) :: dx(:)

      ! local variables
      integer :: i,j 
      real(kind=dp_t) :: kc
      real(kind=dp_t) :: cx_local_plus, cy_local_plus, cx_local_minus, cy_local_minus


      kc = 1.d-4 !make parameter

      !grad x & y
      do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1 !!does the plus one matter here?

              node_grad_cx(i,j) = (c(i,j)-c(i-1,j)+c(i,j-1)-c(i-1,j-1))/(2*dx(1))
              node_grad_cy(i,j) = (c(i,j)-c(i,j-1)+c(i-1,j)-c(i-1,j-1))/(2*dx(2))

       end do
      end do 

      do j=lo(2),hi(2)
      do i=lo(1),hi(1)+1
        cx_local_plus = 0.25*(node_grad_cx(i,j)+node_grad_cx(i,j+1)+node_grad_cx(i+1,j+1)+node_grad_cx(i+1,j))
        cy_local_plus = 0.25*(node_grad_cy(i,j)+node_grad_cy(i,j+1)+node_grad_cy(i+1,j+1)+node_grad_cy(i+1,j))
        cx_local_minus = 0.25*(node_grad_cx(i,j)+node_grad_cx(i,j+1)+node_grad_cx(i-1,j+1)+node_grad_cx(i-1,j))
        cy_local_minus = 0.25*(node_grad_cy(i,j)+node_grad_cy(i,j+1)+node_grad_cy(i-1,j+1)+node_grad_cy(i-1,j))\

!should first line below be ... -node_grad_cx*node_grad_cy ? first term is grad(cxcy) right? also need to check signs 
!also have not added k_c anywhere--need to fix
        forcex(i,j) = (node_grad_cx(i,j+1)*node_grad_cy(i,j+1) - node_grad_cx(i,j)*node_grad_cx(i,j))/dx(2) &
           (0.5*(cy_local_plus*cy_local_plus-cx_local_plus*cx_local_plus) &
           -0.5*(cy_local_minus*cy_local_minus-cx_local_minus*cx_local_minus))/(dx(1)) 
      end do
      end do


      do j=lo(2),hi(2)+1
      do i=lo(1),hi(1)

        cx_local_plus = 0.25*(node_grad_cx(i,j)+node_grad_cx(i,j+1)+node_grad_cx(i+1,j+1)+node_grad_cx(i+1,j))
        cy_local_plus = 0.25*(node_grad_cy(i,j)+node_grad_cy(i,j+1)+node_grad_cy(i+1,j+1)+node_grad_cy(i+1,j))
        cx_local_minus = 0.25*(node_grad_cx(i,j)+node_grad_cx(i,j-1)+node_grad_cx(i+1,j-1)+node_grad_cx(i+1,j))
        cy_local_minus = 0.25*(node_grad_cy(i,j)+node_grad_cy(i,j-1)+node_grad_cy(i+1,j-1)+node_grad_cy(i+1,j))\

        forcey(i,j) = (node_grad_cx(i+1,j)*node_grad_cy(i+1,j) - node_grad_cx(i,j)*node_grad_cx(i,j))/dx(1) &
           (0.5*(cy_local_plus*cy_local_plus-cx_local_plus*cx_local_plus) &
           -0.5*(cy_local_minus*cy_local_minus-cx_local_minus*cx_local_minus))/(dx(2)) 
      end do
      end do

    end subroutine compute_div_reversible_stress_2d

    subroutine compute_div_reversible_stress_3d(forcex,forcey,forcez,ng_1,Ex,Ey,Ez,ng_2,perm,ng_3,lo,hi,dx)

      integer         :: lo(:),hi(:),ng_1,ng_2,ng_3
      real(kind=dp_t) :: forcex(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
      real(kind=dp_t) :: forcey(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
      real(kind=dp_t) :: forcez(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
      real(kind=dp_t) :: node_grad_cx(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_1:)
      real(kind=dp_t) :: node_grad_cy(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_1:)
      real(kind=dp_t) :: node_grad_cz(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_1:)
      real(kind=dp_t) :: Ex(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_2:)
      real(kind=dp_t) :: Ey(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_2:)
      real(kind=dp_t) :: Ez(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_2:)
      real(kind=dp_t) ::   perm(lo(1)-ng_3:,lo(2)-ng_3:,lo(3)-ng_3:)
      real(kind=dp_t) :: dx(:)

      ! local variables
      integer :: i,j,k


      ! local variables
      integer :: i,j 
      real(kind=dp_t) :: kc
      real(kind=dp_t) :: cx_local_plus, cy_local_plus, cz_local_plus, cx_local_minus, cy_local_minus, cz_local_minus


      kc = 1.d-4 !make parameter

      !grad x & y & z
      do k=lo(3),hi(3)
       do j=lo(2),hi(2)
        do i=lo(1),hi(1)+1 !! again, is this correct placement for +1?

              node_grad_cx(i,j,k) = (c(i,j,k)-c(i-1,j,k)+c(i,j-1,k)-c(i-1,j-1,k)+c(i,j,k-1)-c(i-1,j,k-1)+c(i,j-1,k-1)-c(i-1,j-1,k-1))/(4*dx(1))
              node_grad_cy(i,j,k) = (c(i,j,k)-c(i,j-1,k)+c(i-1,j,k)-c(i-1,j-1,k)+c(i,j,k-1)-c(i,j-1,k-1)+c(i-1,j,k-1)-c(i-1,j-1,k-1))/(4*dx(2))
              node_grad_cz(i,j,k) = (c(i,j,k)-c(i,j,k-1)+c(i-1,j,k)-c(i-1,j,k-1)+c(i,j-1,k)-c(i,j-1,k-1)+c(i-1,j-1,k)-c(i-1,j-1,k-1))/(4*dx(3))

        end do
       end do 
      end do

      do k=lo(3),hi(3)
       do j=lo(2),hi(2)
        do i=lo(1),hi(1)+1
         cx_local_plus = (.125)*(node_grad_cx(i,j,k)+node_grad_cx(i,j+1,k)+node_grad_cx(i+1,j+1,k)+node_grad_cx(i+1,j,k)+node_grad_cx(i+1,j+1,k+1)+node_grad_cx(i,j,k+1)+node_grad_cx(i,j+1,k+1)+node_grad_cx(i+1,j,k+1))
         cy_local_plus = (.125)*(node_grad_cy(i,j,k)+node_grad_cy(i,j+1,k)+node_grad_cy(i+1,j+1,k)+node_grad_cy(i+1,j,k)+node_grad_cy(i+1,j+1,k+1)+node_grad_cy(i,j,k+1)+node_grad_cy(i,j+1,k+1)+node_grad_cy(i+1,j,k+1))
         cz_local_plus = (.125)*(node_grad_cz(i,j,k)+node_grad_cz(i,j+1,k)+node_grad_cz(i+1,j+1,k)+node_grad_cz(i+1,j,k)+node_grad_cz(i+1,j+1,k+1)+node_grad_cz(i,j,k+1)+node_grad_cz(i,j+1,k+1)+node_grad_cz(i+1,j,k+1))

         cx_local_minus = (.125)*(node_grad_cx(i,j,k)+node_grad_cx(i,j+1,k)+node_grad_cx(i-1,j+1,k)+node_grad_cx(i-1,j,k)+node_grad_cx(i,j,k+1)+node_grad_cx(i,j+1,k+1)+node_grad_cx(i-1,j+1,k+1)+node_grad_cx(i-1,j,k+1))
         cy_local_minus = (.125)*(node_grad_cy(i,j,k)+node_grad_cy(i,j+1,k)+node_grad_cy(i-1,j+1,k)+node_grad_cy(i-1,j,k)+node_grad_cy(i,j,k+1)+node_grad_cy(i,j+1,k+1)+node_grad_cy(i-1,j+1,k+1)+node_grad_cy(i-1,j,k+1))
         cz_local_minus = (.125)*(node_grad_cz(i,j,k)+node_grad_cz(i,j+1,k)+node_grad_cz(i-1,j+1,k)+node_grad_cz(i-1,j,k)+node_grad_cz(i,j,k+1)+node_grad_cz(i,j+1,k+1)+node_grad_cz(i-1,j+1,k+1)+node_grad_cz(i-1,j,k+1))


        forcex(i,j) = (node_grad_cx(i,j+1)*node_grad_cy(i,j+1) - node_grad_cx(i,j)*node_grad_cx(i,j))/dx(2) &
           (0.5*(cy_local_plus*cy_local_plus-cx_local_plus*cx_local_plus) &
           -0.5*(cy_local_minus*cy_local_minus-cx_local_minus*cx_local_minus))/(dx(1)) 

         forcex(i,j,k) =

         end do
        end do
       end do


      do k=lo(3),hi(3)
      do j=lo(2),hi(2)+1
      do i=lo(1),hi(1)

         cx_local_plus = (.125)*(node_grad_cx(i,j,k)+node_grad_cx(i,j+1,k)+node_grad_cx(i+1,j+1,k)+node_grad_cx(i+1,j,k)+node_grad_cx(i+1,j+1,k+1)+node_grad_cx(i,j,k+1)+node_grad_cx(i,j+1,k+1)+node_grad_cx(i+1,j,k+1))
         cy_local_plus = (.125)*(node_grad_cy(i,j,k)+node_grad_cy(i,j+1,k)+node_grad_cy(i+1,j+1,k)+node_grad_cy(i+1,j,k)+node_grad_cy(i+1,j+1,k+1)+node_grad_cy(i,j,k+1)+node_grad_cy(i,j+1,k+1)+node_grad_cy(i+1,j,k+1))
         cz_local_plus = (.125)*(node_grad_cz(i,j,k)+node_grad_cz(i,j+1,k)+node_grad_cz(i+1,j+1,k)+node_grad_cz(i+1,j,k)+node_grad_cz(i+1,j+1,k+1)+node_grad_cz(i,j,k+1)+node_grad_cz(i,j+1,k+1)+node_grad_cz(i+1,j,k+1))

         cx_local_minus = (.125)*(node_grad_cx(i,j,k)+node_grad_cx(i,j-1,k)+node_grad_cx(i+1,j-1,k)+node_grad_cx(i+1,j,k)+node_grad_cx(i,j,k+1)+node_grad_cx(i,j-1,k+1)+node_grad_cx(i+1,j-1,k+1)+node_grad_cx(i+1,j,k+1))
         cy_local_minus = (.125)*(node_grad_cy(i,j,k)+node_grad_cy(i,j-1,k)+node_grad_cy(i+1,j-1,k)+node_grad_cy(i+1,j,k)+node_grad_cy(i,j,k+1)+node_grad_cy(i,j-1,k+1)+node_grad_cy(i+1,j-1,k+1)+node_grad_cy(i+1,j,k+1))
         cz_local_minus = (.125)*(node_grad_cz(i,j,k)+node_grad_cz(i,j-1,k)+node_grad_cz(i+1,j-1,k)+node_grad_cz(i+1,j,k)+node_grad_cz(i,j,k+1)+node_grad_cz(i,j-1,k+1)+node_grad_cz(i+1,j-1,k+1)+node_grad_cz(i+1,j,k+1))

        forcey(i,j,k) = 

      end do
      end do
      end do

      do k=lo(3),hi(3)+1
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)

         cx_local_plus = (.125)*(node_grad_cx(i,j,k)+node_grad_cx(i,j+1,k)+node_grad_cx(i+1,j+1,k)+node_grad_cx(i+1,j,k)+node_grad_cx(i+1,j+1,k+1)+node_grad_cx(i,j,k+1)+node_grad_cx(i,j+1,k+1)+node_grad_cx(i+1,j,k+1))
         cy_local_plus = (.125)*(node_grad_cy(i,j,k)+node_grad_cy(i,j+1,k)+node_grad_cy(i+1,j+1,k)+node_grad_cy(i+1,j,k)+node_grad_cy(i+1,j+1,k+1)+node_grad_cy(i,j,k+1)+node_grad_cy(i,j+1,k+1)+node_grad_cy(i+1,j,k+1))
         cz_local_plus = (.125)*(node_grad_cz(i,j,k)+node_grad_cz(i,j+1,k)+node_grad_cz(i+1,j+1,k)+node_grad_cz(i+1,j,k)+node_grad_cz(i+1,j+1,k+1)+node_grad_cz(i,j,k+1)+node_grad_cz(i,j+1,k+1)+node_grad_cz(i+1,j,k+1))

         cx_local_minus = (.125)*(node_grad_cx(i,j,k)+node_grad_cx(i,j+1,k)+node_grad_cx(i+1,j+1,k)+node_grad_cx(i+1,j,k)+node_grad_cx(i,j,k-1)+node_grad_cx(i,j+1,k-1)+node_grad_cx(i+1,j+1,k-1)+node_grad_cx(i+1,j,k-1))
         cy_local_minus = (.125)*(node_grad_cy(i,j,k)+node_grad_cy(i,j+1,k)+node_grad_cy(i+1,j+1,k)+node_grad_cy(i+1,j,k)+node_grad_cy(i,j,k-1)+node_grad_cy(i,j+1,k-1)+node_grad_cy(i+1,j+1,k-1)+node_grad_cy(i+1,j,k-1))
         cz_local_minus = (.125)*(node_grad_cz(i,j,k)+node_grad_cz(i,j+1,k)+node_grad_cz(i+1,j+1,k)+node_grad_cz(i+1,j,k)+node_grad_cz(i,j,k-1)+node_grad_cz(i,j+1,k-1)+node_grad_cz(i+1,j+1,k-1)+node_grad_cz(i+1,j,k-1))

        forcez(i,j,k) = 

      end do
      end do
      end do

!!!!old code

      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)+1
         forcex(i,j,k) = forcex(i,j,k) - 0.5d0 * ( Ex(i,j,k)**2 &
              + (0.25d0*(Ey(i,j,k)+Ey(i,j+1,k)+Ey(i-1,j,k)+Ey(i-1,j+1,k)))*2 &
              + (0.25d0*(Ez(i,j,k)+Ez(i,j,k+1)+Ez(i-1,j,k)+Ez(i-1,j,k+1)))*2 ) &
              * (perm(i,j,k)-perm(i-1,j,k)) / dx(1)
      end do
      end do
      end do

      do k=lo(3),hi(3)
      do j=lo(2),hi(2)+1
      do i=lo(1),hi(1)
         forcey(i,j,k) = forcey(i,j,k) - 0.5d0 * ( Ey(i,j,k)**2 &
              + (0.25d0*(Ex(i,j,k)+Ex(i+1,j,k)+Ex(i,j-1,k)+Ex(i+1,j-1,k)))**2 &
              + (0.25d0*(Ez(i,j,k)+Ez(i,j,k+1)+Ez(i,j-1,k)+Ez(i,j-1,k+1)))**2 ) &
              * (perm(i,j,k)-perm(i,j-1,k)) / dx(2)
      end do
      end do
      end do

      do k=lo(3),hi(3)+1
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)
         forcez(i,j,k) = forcez(i,j,k) - 0.5d0 * ( Ez(i,j,k)**2 &
              + (0.25d0*(Ex(i,j,k)+Ex(i+1,j,k)+Ex(i,j,k-1)+Ex(i+1,j,k-1)))**2 &
              + (0.25d0*(Ey(i,j,k)+Ey(i,j+1,k)+Ey(i,j,k-1)+Ey(i,j+1,k-1)))**2 ) &
              * (perm(i,j,k)-perm(i,j,k-1)) / dx(3)
      end do
      end do
      end do

    end subroutine compute_div_reversible_stress_3d

  end subroutine compute_div_reversible_stress

  subroutine compute_E_ext(mla,E_ext)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: E_ext(:,:)

    integer :: n, nlevs, dm, i

    nlevs = mla%nlevel
    dm = mla%dim

    if (E_ext_type .eq. 1) then

       do n=1,nlevs
          do i=1,dm
             call multifab_setval(E_ext(n,i),E_ext_value(i),all=.true.)
          end do
       end do

    end if

  end subroutine compute_E_ext

  subroutine zero_eps_on_wall(mla,beta,dx)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: beta(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)

    integer :: i,dm,nlevs,n,ng_b
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: bpx(:,:,:,:)
    real(kind=dp_t), pointer :: bpy(:,:,:,:)
    real(kind=dp_t), pointer :: bpz(:,:,:,:)

    nlevs = mla%nlevel
    dm = mla%dim

    ng_b = beta(1,1)%ng

    do n=1,nlevs
       do i=1,nfabs(beta(n,1))
          lo = lwb(get_box(beta(n,1),i))
          hi = upb(get_box(beta(n,1),i))
          bpy => dataptr(beta(n,2),i)
          select case (dm)
          case (2)
             call zero_eps_on_wall_2d(bpy(:,:,1,1),ng_b,lo,hi,dx(n,:))
          case (3)
             call bl_error("zero_eps_on_wall_3d not written yet")
          end select
       end do
    end do

  end subroutine zero_eps_on_wall

  subroutine zero_eps_on_wall_2d(betay,ng_b,lo,hi,dx)
      
    integer         :: lo(:),hi(:),ng_b
    real(kind=dp_t) :: betay(lo(1)-ng_b:,lo(2)-ng_b:)
    real(kind=dp_t) :: dx(:)

    integer :: i,j

    real(kind=dp_t) :: x,y,Lx

    Lx = prob_hi(1) - prob_lo(1)    

    if (zero_eps_on_wall_type .eq. 1) then

       ! y-faces
       do j=lo(2),hi(2)+1
          y = prob_lo(2) + dx(2)*j
          do i=lo(1),hi(1)
             x = prob_lo(1) + dx(1)*(i+0.5d0)

             if (y .eq. 0.d0 .and. (x .le. zero_eps_on_wall_left_end*Lx .or. x .ge. zero_eps_on_wall_right_start*Lx) ) then
                betay(i,j) = 0.d0
             end if

          end do
       end do

    end if

  end subroutine zero_eps_on_wall_2d

end module fluid_charge_module
