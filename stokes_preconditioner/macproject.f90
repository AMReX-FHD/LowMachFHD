module macproject_module

  use ml_layout_module
  use define_bc_module
  use bl_constants_module
  use bl_error_module
  use bc_module 
  use convert_module

  implicit none

  private

  public :: macproject, subtract_weighted_gradp

contains 

  ! solve L_alpha Phi  = D x_u^* - b_p
  ! does not update any other variables
  subroutine macproject(mla,phi,umac,alpha,mac_rhs,dx,the_bc_tower)

    use probin_module            , only : mg_rel_tol, mg_max_vcycles, num_mg_vcycles
    use mac_multigrid_module     , only : mac_multigrid
    use bndry_reg_module

    type(ml_layout), intent(in   ) :: mla
    type(multifab ), intent(inout) :: phi(:)
    type(multifab ), intent(in   ) :: umac(:,:)
    type(multifab ), intent(in   ) :: alpha(:)
    type(multifab ), intent(inout) :: mac_rhs(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower ), intent(in   ) :: the_bc_tower

    ! Local  
    type(multifab)  :: zero_fab(mla%nlevel)
    type(multifab)  :: alphainv_edge(mla%nlevel,mla%dim)
    type(bndry_reg) :: fine_flx(2:mla%nlevel)
    real(kind=dp_t) :: umac_norm(mla%nlevel)
    real(kind=dp_t) :: rel_solver_eps
    real(kind=dp_t) :: abs_solver_eps
    integer         :: d,dm,i,n,nlevs,bc_comp

    nlevs = mla%nlevel
    dm = mla%dim

    bc_comp = dm + 1

    do n = 1, nlevs
       call multifab_build(zero_fab(n), mla%la(n),  1, 0)
       do d = 1,dm
          call multifab_build_edge(alphainv_edge(n,d), mla%la(n), 1, 0, d)
       end do

       call setval(zero_fab(n),ZERO,all=.true.)
       call setval(  phi(n),ZERO,all=.true.)

    end do

    ! Compute umac_norm to be used inside the MG solver as part of a stopping criterion
    umac_norm = -1.0_dp_t
    do n = 1,nlevs
       do i = 1,dm
          umac_norm(n) = max(umac_norm(n),norm_inf(umac(n,i)))
       end do
    end do

    ! compute alphainv_edge on faces by averaging and then inverting
    call average_cc_to_face_inv(nlevs,alpha,alphainv_edge,1,dm+2,1,the_bc_tower%bc_tower_array)

    ! multiply alphainv_edge by -1 so we solve L_alpha Phi = mac_rhs
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_s_c(alphainv_edge(n,i),1,-1.d0,1,0)
       end do
    end do

    ! stores (alphainv_edge/dx**2) grad phi at coarse-fine interfaces
    do n = 2,nlevs
       call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
    end do

    ! loosen the tolerance for multilevel problems
    if (nlevs .eq. 1) then
       rel_solver_eps = mg_rel_tol
    else if (nlevs .eq. 2) then
       rel_solver_eps = 1.0d1*mg_rel_tol
    else
       rel_solver_eps = 1.0d2*mg_rel_tol
    endif

    abs_solver_eps = 1.d-16

    call mac_multigrid(mla,mac_rhs,phi,fine_flx,zero_fab,alphainv_edge,dx,the_bc_tower,bc_comp, &
                       mla%mba%rr,rel_solver_eps,abs_solver_eps, &
                       mg_max_vcycles_in=mg_max_vcycles,abort_on_max_iter_in=.false.)

    num_mg_vcycles = num_mg_vcycles + mg_max_vcycles

    do n = 1, nlevs
       call multifab_destroy(zero_fab(n))
       do d = 1,dm
          call multifab_destroy(alphainv_edge(n,d))
       end do
    end do

    do n = 2,nlevs
       call bndry_reg_destroy(fine_flx(n))
    end do

  end subroutine macproject
  
  subroutine subtract_weighted_gradp(mla,x_u,alphainv_edge,phi,dx)

    type(ml_layout), intent(in   ) :: mla
    type(multifab ), intent(inout) :: x_u(:,:)
    type(multifab ), intent(in   ) :: alphainv_edge(:,:)
    type(multifab ), intent(in   ) :: phi(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)

    ! local
    integer :: i,n,dm,nlevs
    integer :: ng_u,ng_a,ng_p
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: vp(:,:,:,:)
    real(kind=dp_t), pointer :: wp(:,:,:,:)
    real(kind=dp_t), pointer :: axp(:,:,:,:)
    real(kind=dp_t), pointer :: ayp(:,:,:,:)
    real(kind=dp_t), pointer :: azp(:,:,:,:)
    real(kind=dp_t), pointer :: pp(:,:,:,:)

    nlevs = mla%nlevel
    dm = mla%dim

    ng_u = x_u(1,1)%ng
    ng_a = alphainv_edge(1,1)%ng
    ng_p = phi(1)%ng

    do n=1,nlevs
       do i=1,nfabs(phi(n))
          up => dataptr(x_u(n,1),i)
          vp => dataptr(x_u(n,2),i)
          axp => dataptr(alphainv_edge(n,1),i)
          ayp => dataptr(alphainv_edge(n,2),i)
          pp => dataptr(phi(n),i)
          lo = lwb(get_box(phi(n), i))
          hi = upb(get_box(phi(n), i))
          select case (dm)
          case (2)
             call subtract_weighted_gradp_2d(up(:,:,1,1),vp(:,:,1,1),ng_u, &
                                             axp(:,:,1,1),ayp(:,:,1,1),ng_a, &
                                             pp(:,:,1,1),ng_p,lo,hi,dx(n,:))
          case (3)
             wp => dataptr(x_u(n,3),i)
             azp => dataptr(alphainv_edge(n,3),i)
             call subtract_weighted_gradp_3d(up(:,:,:,1),vp(:,:,:,1),wp(:,:,:,1),ng_u, &
                                             axp(:,:,:,1),ayp(:,:,:,1),azp(:,:,:,1),ng_a, &
                                             pp(:,:,:,1),ng_p,lo,hi,dx(n,:))
          end select
       end do
    end do
    
  contains

    subroutine subtract_weighted_gradp_2d(xu,xv,ng_u,ainv_edgex,ainv_edgey,ng_a,phi,ng_p,lo,hi,dx)

      integer        , intent(in   ) :: lo(:),hi(:),ng_u,ng_a,ng_p
      real(kind=dp_t), intent(inout) ::         xu(lo(1)-ng_u:,lo(2)-ng_u:)
      real(kind=dp_t), intent(inout) ::         xv(lo(1)-ng_u:,lo(2)-ng_u:)
      real(kind=dp_t), intent(in   ) :: ainv_edgex(lo(1)-ng_a:,lo(2)-ng_a:)
      real(kind=dp_t), intent(in   ) :: ainv_edgey(lo(1)-ng_a:,lo(2)-ng_a:)
      real(kind=dp_t), intent(in   ) ::        phi(lo(1)-ng_p:,lo(2)-ng_p:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      ! local
      integer :: i,j

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
            xu(i,j) = xu(i,j) - ainv_edgex(i,j)*(phi(i,j)-phi(i-1,j))/dx(1)
         end do
      end do
         
      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)
            xv(i,j) = xv(i,j) - ainv_edgey(i,j)*(phi(i,j)-phi(i,j-1))/dx(2)
         end do
      end do

    end subroutine subtract_weighted_gradp_2d

    subroutine subtract_weighted_gradp_3d(xu,xv,xw,ng_u,ainv_edgex,ainv_edgey,ainv_edgez,ng_a,phi,ng_p,lo,hi,dx)

      integer        , intent(in   ) :: lo(:),hi(:),ng_u,ng_a,ng_p
      real(kind=dp_t), intent(inout) ::         xu(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
      real(kind=dp_t), intent(inout) ::         xv(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
      real(kind=dp_t), intent(inout) ::         xw(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
      real(kind=dp_t), intent(in   ) :: ainv_edgex(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
      real(kind=dp_t), intent(in   ) :: ainv_edgey(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
      real(kind=dp_t), intent(in   ) :: ainv_edgez(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
      real(kind=dp_t), intent(in   ) ::        phi(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      ! local
      integer :: i,j,k

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)+1
               xu(i,j,k) = xu(i,j,k) - ainv_edgex(i,j,k)*(phi(i,j,k)-phi(i-1,j,k))/dx(1)
            end do
         end do
      end do

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)
               xv(i,j,k) = xv(i,j,k) - ainv_edgey(i,j,k)*(phi(i,j,k)-phi(i,j-1,k))/dx(2)
            end do
         end do
      end do

      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               xw(i,j,k) = xw(i,j,k) - ainv_edgez(i,j,k)*(phi(i,j,k)-phi(i,j,k-1))/dx(3)
            end do
         end do
      end do

    end subroutine subtract_weighted_gradp_3d

  end subroutine subtract_weighted_gradp

end module macproject_module
