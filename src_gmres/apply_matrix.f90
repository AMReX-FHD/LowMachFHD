module apply_matrix_module

  use ml_layout_module
  use multifab_module
  use define_bc_module
  use probin_module
  use div_and_grad_module
  use stag_applyop_module
  use div_and_grad_module
  use multifab_physbc_module

  implicit none

  private

  public :: apply_matrix

contains

  ! This computes A x = b explicitly
  ! Refer to ./doc/PreconditionerNotes.tex
  subroutine apply_matrix(mla,b_u,b_p,x_u,x_p,alpha,beta,gamma,theta, &
                          dx,the_bc_tower,use_inhomogeneous_in)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: b_u(:,:)
    type(multifab) , intent(inout) :: b_p(:)
    type(multifab) , intent(inout) :: x_u(:,:)
    type(multifab) , intent(inout) :: x_p(:)
    type(multifab) , intent(inout) :: alpha(:)
    type(multifab) , intent(in   ) :: beta(:)
    type(multifab) , intent(in   ) :: gamma(:)
    real(kind=dp_t), intent(in   ) :: theta
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    logical,  intent(in), optional :: use_inhomogeneous_in

    ! local
    integer :: n,nlevs,i,dm
    type(multifab) ::  gx_p(mla%nlevel,mla%dim)
    logical :: use_inhomogeneous

    if (present(use_inhomogeneous_in)) then
       use_inhomogeneous = use_inhomogeneous_in
    else
       use_inhomogeneous = .false.
    end if

    nlevs = mla%nlevel
    dm = mla%dim

    if (x_u(1,1)%ng .lt. 1) then
       call bl_error("apply_matrix.f90: x_u needs at least 1 ghost cell")
    end if

    if (x_p(1)%ng .lt. 1) then
       call bl_error("apply_matrix.f90: x_p needs at least 1 ghost cell")
    end if

    ! fill ghost cells for x_u and x_p
    do n=1,nlevs
       call multifab_fill_boundary(x_p(n))
       call multifab_physbc(x_p(n),1,dm+1,1,the_bc_tower%bc_tower_array(n))
       do i=1,dm
          if (use_inhomogeneous) then
             call multifab_physbc_domainvel_inhomogeneous(x_u(n,i),1,i,1,the_bc_tower%bc_tower_array(n),dx(n,:))             
          else
             call multifab_physbc_domainvel(x_u(n,i),1,i,1,the_bc_tower%bc_tower_array(n),dx(n,:))
          end if
          call multifab_fill_boundary(x_u(n,i))
          if (use_inhomogeneous) then
             call multifab_physbc_macvel_inhomogeneous(x_u(n,i),1,i,1,the_bc_tower%bc_tower_array(n),dx(n,:))
          else
             call multifab_physbc_macvel(x_u(n,i),1,i,1,the_bc_tower%bc_tower_array(n),dx(n,:))
          end if
       end do
    end do

    do n=1,nlevs
       do i=1,dm
          call multifab_build_edge(gx_p(n,i),mla%la(n),1,0,i)
       end do
    end do

    ! compute b_u = A x_u
    call stag_applyop(mla,the_bc_tower,x_u,b_u,alpha,beta,gamma,theta,dx)

    ! compute G x_p and add to b_u
    call compute_gradp(mla,x_p,gx_p,dx)
    
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(b_u(n,i),1,gx_p(n,i),1,1,0)
       end do
    end do
        
    ! set b_p = -D x_u
    call compute_divu(mla,x_u,b_p,dx)
    do n=1,nlevs
       call multifab_mult_mult_s_c(b_p(n),1,-1.d0,1,0)
    end do

    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(gx_p(n,i))
       end do
    end do

  end subroutine apply_matrix

end module apply_matrix_module
