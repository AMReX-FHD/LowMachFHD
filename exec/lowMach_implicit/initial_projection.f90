module initial_projection_module

  use multifab_module
  use ml_layout_module
  use convert_stag_module
  use convert_variables_module
  use define_bc_module
  use macproject_module
  use div_and_grad_module
  use mk_diffusive_fluxdiv_module
  use mk_stochastic_fluxdiv_module
  use multifab_physbc_module
  use probin_lowmach_module, only: rhobar, diff_coef, nscal

  implicit none

  private

  public :: initial_projection

contains

  subroutine initial_projection(mla,mold,umac,sold,s_fc,prim,chi_fc,rhoc_d_fluxdiv, &
                                rhoc_s_fluxdiv,dx,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: mold(:,:)
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: s_fc(:,:)
    type(multifab) , intent(in   ) :: prim(:)
    type(multifab) , intent(in   ) :: chi_fc(:,:)
    type(multifab) , intent(inout) :: rhoc_d_fluxdiv(:)
    type(multifab) , intent(inout) :: rhoc_s_fluxdiv(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: i,dm,n,nlevs
    real(kind=dp_t) :: S_fac

    type(multifab) ::   mac_rhs(mla%nlevel)
    type(multifab) ::      divu(mla%nlevel)
    type(multifab) ::       phi(mla%nlevel)
    type(multifab) :: rhoinv_fc(mla%nlevel,mla%dim)

    dm = mla%dim
    nlevs = mla%nlevel

    S_fac = (1.d0/rhobar(1) - 1.d0/rhobar(2))

    do n=1,nlevs
       call multifab_build(mac_rhs(n),mla%la(n),1,0)
       call multifab_build(divu(n),mla%la(n),1,0)
       call multifab_build(phi(n),mla%la(n),1,1)
       do i=1,dm
          call multifab_build_edge(rhoinv_fc(n,i),mla%la(n),1    ,0,i)
       end do       
    end do

    do n=1,nlevs
       call setval(mac_rhs(n),0.d0)
       call setval(phi(n),0.d0)
    end do

    ! average sold to faces
    do i=1,nscal
       call average_cc_to_face(nlevs,sold,s_fc,i,dm+2,1,the_bc_tower%bc_tower_array)
    end do

    ! convert m^init to v^init in valid region
    call convert_m_to_umac(mla,s_fc,mold,umac,.true.)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! build rhs = div(v^init) - S^0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! set rhoc_d_fluxdiv = div(rho*chi grad c)^0
    call mk_diffusive_rhoc_fluxdiv(mla,rhoc_d_fluxdiv,1,prim,s_fc,chi_fc,dx, &
                                   the_bc_tower%bc_tower_array)

    ! set mac_rhs to div(rho*chi grad c)^0
    do n=1,nlevs
       call multifab_plus_plus_c(mac_rhs(n),1,rhoc_d_fluxdiv(n),1,1,0)
    end do

    ! set rhoc_s_fluxdiv = div(Psi^0)
    call mk_stochastic_s_fluxdiv(mla,the_bc_tower%bc_tower_array,rhoc_s_fluxdiv,s_fc, &
                                 chi_fc,dx,1)

    ! add div(Psi^0) to mac_rhs
    do n=1,nlevs
       call multifab_plus_plus_c(mac_rhs(n),1,rhoc_s_fluxdiv(n),1,1,0)
    end do

    ! multiply mac_rhs by -S_fac
    ! now mac_rhs = -S^0
    do n=1,nlevs
       call multifab_mult_mult_s_c(mac_rhs(n),1,-S_fac,1,0)
    end do

    ! set divu = div(v^init)
    call compute_divu(mla,umac,divu,dx)

    ! add div(v^init) to mac_rhs
    ! now mac_rhs = div(v^init) - S^0
    do n=1,nlevs
       call multifab_plus_plus_c(mac_rhs(n),1,divu(n),1,1,0)
    end do

    ! compute (1/rho^0)
    do n=1,nlevs
       do i=1,dm
          call setval(rhoinv_fc(n,i),1.d0,all=.true.)
          call multifab_div_div_c(rhoinv_fc(n,i),1,s_fc(n,i),1,1,0)
       end do
    end do

    ! solve div (1/rho^0) grad phi = div(v^init) - S^0
    ! solve to completion, i.e., use the 'full' solver
    call macproject(mla,phi,umac,rhoinv_fc,mac_rhs,dx,the_bc_tower,.true.)

    ! v^0 = v^init - (1/rho^0) grad phi
    call subtract_weighted_gradp(mla,umac,rhoinv_fc,phi,dx)

    ! fill ghost cells
    do n=1,nlevs
       do i=1,dm
          ! fill periodic and interior ghost cells
          call multifab_fill_boundary(umac(n,i))
          ! set normal velocity on physical domain boundaries
          call multifab_physbc_domainvel(umac(n,i),1,i,1, &
                                         the_bc_tower%bc_tower_array(n),dx(n,:))
          ! set transverse velocity behind physical boundaries
          call multifab_physbc_macvel(umac(n,i),1,i,1, &
                                      the_bc_tower%bc_tower_array(n),dx(n,:))
       end do
    end do

    ! compute mold by converting v^0 to m^0 in valid plus ghost region
    ! now mold has properly filled ghost cells
    call convert_m_to_umac(mla,s_fc,mold,umac,.false.)

    do n=1,nlevs
       call multifab_destroy(mac_rhs(n))
       call multifab_destroy(divu(n))
       call multifab_destroy(phi(n))
       do i=1,dm
          call multifab_destroy(rhoinv_fc(n,i))
       end do
    end do

  end subroutine initial_projection

end module initial_projection_module
