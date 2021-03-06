module apply_matrix_module

  use ml_layout_module
  use multifab_module
  use define_bc_module
  use div_and_grad_module
  use stag_applyop_module
  use div_and_grad_module
  use bc_module
  use multifab_physbc_module
  use multifab_physbc_stag_module
  use probin_gmres_module, only: gmres_spatial_order

  implicit none

  private

  public :: apply_matrix

contains

  ! This computes A x = b explicitly
  ! Refer to ./doc/PreconditionerNotes.tex
  subroutine apply_matrix(mla,b_u,b_p,x_u,x_p,alpha_fc,beta,beta_ed,gamma,theta_alpha, &
                          dx,the_bc_tower,vel_bc_n,vel_bc_t)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: b_u(:,:)
    type(multifab) , intent(inout) :: b_p(:)
    type(multifab) , intent(inout) :: x_u(:,:)
    type(multifab) , intent(inout) :: x_p(:)
    type(multifab) , intent(in   ) :: alpha_fc(:,:)
    type(multifab) , intent(in   ) :: beta(:)
    type(multifab) , intent(in   ) :: beta_ed(:,:)
    type(multifab) , intent(in   ) :: gamma(:)
    real(kind=dp_t), intent(in   ) :: theta_alpha
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(in), optional :: vel_bc_n(:,:)
    type(multifab) , intent(in), optional :: vel_bc_t(:,:)

    ! local
    integer :: n,nlevs,i,dm
    type(multifab) ::  gx_p(mla%nlevel,mla%dim)

    type(bl_prof_timer), save :: bpt
    
    call build(bpt, "apply_matrix")

    nlevs = mla%nlevel
    dm = mla%dim

    ! check to make sure x_u and x_p have enough ghost cells
    if (gmres_spatial_order .eq. 2) then
       if (x_u(1,1)%ng .lt. 1) then
          call bl_error("apply_matrix.f90: x_u needs at least 1 ghost cell")
       end if
       if (x_p(1)%ng .lt. 1) then
          call bl_error("apply_matrix.f90: x_p needs at least 1 ghost cell")
       end if
    else if (gmres_spatial_order .eq. 4) then
       if (x_u(1,1)%ng .lt. 2) then
          call bl_error("apply_matrix.f90: x_u needs at least 2 ghost cells")
       end if
       if (x_p(1)%ng .lt. 2) then
          call bl_error("apply_matrix.f90: x_p needs at least 2 ghost cells")
       end if
    end if

    ! fill ghost cells for x_u and x_p
    do n=1,nlevs
       call multifab_fill_boundary(x_p(n))
       if (gmres_spatial_order .eq. 2) then
          call multifab_physbc(x_p(n),1,pres_bc_comp,1,the_bc_tower%bc_tower_array(n), &
                               dx_in=dx(n,:))
       else if (gmres_spatial_order .eq. 4) then
          call multifab_physbc_ho(x_p(n),1,pres_bc_comp,1,the_bc_tower%bc_tower_array(n), &
                                  dx_in=dx(n,:))
       end if
       do i=1,dm
          ! these if present tests are a workaround for an intel compiler bug where
          ! passing optional arguments through causes crashes on edison
          if (present(vel_bc_n)) then
             call multifab_physbc_domainvel(x_u(n,i),vel_bc_comp+i-1, &
                                            the_bc_tower%bc_tower_array(n), &
                                            dx(n,:),vel_bc_n(n,:))
          else
             if (gmres_spatial_order .eq. 2) then
                call multifab_physbc_domainvel(x_u(n,i),vel_bc_comp+i-1, &
                                               the_bc_tower%bc_tower_array(n), &
                                               dx(n,:))
             else if (gmres_spatial_order .eq. 4) then
                call multifab_physbc_domainvel_ho(x_u(n,i),vel_bc_comp+i-1, &
                                                  the_bc_tower%bc_tower_array(n), &
                                                  dx(n,:))
             end if
          end if
          call multifab_fill_boundary(x_u(n,i))
          if (present(vel_bc_t)) then
             call multifab_physbc_macvel(x_u(n,i),vel_bc_comp+i-1, &
                                         the_bc_tower%bc_tower_array(n), &
                                         dx(n,:),vel_bc_t(n,:))
          else
             if (gmres_spatial_order .eq. 2) then
                call multifab_physbc_macvel(x_u(n,i),vel_bc_comp+i-1, &
                                            the_bc_tower%bc_tower_array(n), &
                                            dx(n,:))
             else if (gmres_spatial_order .eq. 4) then
                call multifab_physbc_macvel_ho(x_u(n,i),vel_bc_comp+i-1, &
                                               the_bc_tower%bc_tower_array(n), &
                                               dx(n,:))
             end if
          end if
       end do
    end do

    do n=1,nlevs
       do i=1,dm
          call multifab_build_edge(gx_p(n,i),mla%la(n),1,0,i)
       end do
    end do

    ! compute b_u = A x_u
    if (gmres_spatial_order .eq. 2) then
       call stag_applyop(mla,the_bc_tower,x_u,b_u,alpha_fc,beta,beta_ed,gamma,theta_alpha,dx)
    else if (gmres_spatial_order .eq. 4) then
       call stag_applyop_4o(mla,the_bc_tower,x_u,b_u,alpha_fc,beta,beta_ed,gamma,theta_alpha,dx)
    end if

    ! compute G x_p and add to b_u
    if (gmres_spatial_order .eq. 2) then
       call compute_grad(mla,x_p,gx_p,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)
    else if (gmres_spatial_order .eq. 4) then
       call compute_grad_4o(mla,x_p,gx_p,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)
    end if
    
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(b_u(n,i),1,gx_p(n,i),1,1,0)
       end do
    end do
        
    ! set b_p = -D x_u
    call compute_div(mla,x_u,b_p,dx,1,1,1)
    do n=1,nlevs
       call multifab_mult_mult_s_c(b_p(n),1,-1.d0,1,0)
    end do

    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(gx_p(n,i))
       end do
    end do

    call destroy(bpt)

  end subroutine apply_matrix

end module apply_matrix_module
