module initial_projection_module

  use multifab_module
  use ml_layout_module
  use convert_stag_module
  use convert_variables_module
  use define_bc_module
  use macproject_module
  use div_and_grad_module
  use mk_diffusive_fluxdiv_module
  use multifab_physbc_module

  implicit none

  private

  public :: initial_projection

contains

  subroutine initial_projection(mla,mold,umac,sold,chi,dx,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: mold(:,:)
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(inout) :: sold(:)
    type(multifab) , intent(in   ) :: chi(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: i,dm,n,nlevs

    type(multifab) ::     mac_rhs(mla%nlevel)
    type(multifab) ::         phi(mla%nlevel)
    type(multifab) ::    rho_face(mla%nlevel,mla%dim)
    type(multifab) :: rhoinv_face(mla%nlevel,mla%dim)
    type(multifab) ::    chi_face(mla%nlevel,mla%dim)

    dm = mla%dim
    nlevs = mla%nlevel

    do n=1,nlevs
       call multifab_build(mac_rhs(n),mla%la(n),1,0)
       call multifab_build(phi(n),mla%la(n),1,1)
       call setval(phi(n),0.d0)
       do i=1,dm
          call multifab_build_edge(   rho_face(n,i),mla%la(n),1,0,i)
          call multifab_build_edge(rhoinv_face(n,i),mla%la(n),1,0,i)
          call multifab_build_edge(   chi_face(n,i),mla%la(n),1,0,i)
       end do       
    end do

    ! create average rho to faces
    call average_cc_to_face(nlevs,sold,rho_face,1,dm+2,1,the_bc_tower%bc_tower_array)

    ! convert m to u
    call convert_m_to_umac(mla,rho_face,mold,umac,.true.)

    !!!!!!!!!!!!!!!!!!!!!!!!!
    ! build rhs = S - div(u)
    !!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute mac_rhs = -div(u)
    call compute_divu(mla,umac,mac_rhs,dx)
    do n=1,nlevs
       call multifab_mult_mult_s_c(mac_rhs(n),1,-1.d0,1,0)
    end do

    ! average chi to faces
    call average_cc_to_face(nlevs,chi,chi_face,1,dm+2,1,the_bc_tower%bc_tower_array)

    ! temporary way to convert rho*c to c
    do n=1,nlevs
       call multifab_div_div_c(sold(n),1,sold(n),2,1,1)
    end do

    ! add del dot rho chi grad c
    call mk_diffusive_rhoc_fluxdiv(mla,mac_rhs,1,sold,rho_face,chi_face,dx, &
                                   the_bc_tower%bc_tower_array)

    ! convert c back to rho*c
    do n=1,nlevs
       call multifab_mult_mult_c(sold(n),1,sold(n),2,1,1)
    end do

    ! project to solve for phi - use the 'full' solver
    call macproject(mla,phi,umac,sold,mac_rhs,dx,the_bc_tower,.true.)

    ! compute (1/rho)
    call average_cc_to_face_inv(nlevs,sold,rhoinv_face,1,dm+2,1,the_bc_tower%bc_tower_array)

    ! umac = umac - (1/rho) grad phi
    call subtract_weighted_gradp(mla,umac,rhoinv_face,phi,dx)

    ! fill periodic/interior and physical ghost cells for umac
    do n=1,nlevs
       do i=1,dm
          call multifab_fill_boundary(umac(n,i))
       end do
    end do

    ! convert u to m
    call convert_m_to_umac(mla,rho_face,mold,umac,.false.)

    do n=1,nlevs
       call destroy(mac_rhs(n))
       call destroy(phi(n))
       do i=1,dm
          call destroy(rho_face(n,i))
          call destroy(rhoinv_face(n,i))
          call destroy(chi_face(n,i))
       end do
    end do

  end subroutine initial_projection

end module initial_projection_module
