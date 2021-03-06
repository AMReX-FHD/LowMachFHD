module initial_projection_module

  use multifab_module
  use ml_layout_module
  use convert_stag_module
  use define_bc_module
  use macproject_module
  use div_and_grad_module
  use bc_module
  use multifab_physbc_stag_module
  use compute_mass_fluxdiv_energy_module
  use convert_rhoc_to_c_module
  use eos_model_wrapper_module
  use compute_rhoh_fluxdiv_module
  use probin_common_module, only: nspecies, restart, total_volume

  implicit none

  private

  public :: initial_projection

contains

  subroutine initial_projection(mla,umac,rho,rhotot,gradp_baro, &
                                diff_mass_fluxdiv, &
                                Temp,p0,dt,dx,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(inout) :: rho(:)
    type(multifab) , intent(inout) :: rhotot(:)
    type(multifab) , intent(in   ) :: gradp_baro(:,:)
    type(multifab) , intent(inout) :: diff_mass_fluxdiv(:)
    type(multifab) , intent(in   ) :: Temp(:)
    real(kind=dp_t), intent(in   ) :: p0
    real(kind=dp_t), intent(in   ) :: dt
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    ! local
    integer :: i,dm,n,nlevs

    type(multifab) ::              conc(mla%nlevel)
    type(multifab) ::         molarconc(mla%nlevel)
    type(multifab) ::               eta(mla%nlevel)
    type(multifab) ::            lambda(mla%nlevel)
    type(multifab) ::             kappa(mla%nlevel)
    type(multifab) ::               chi(mla%nlevel)
    type(multifab) ::              zeta(mla%nlevel)
    type(multifab) ::           mac_rhs(mla%nlevel)
    type(multifab) ::             theta(mla%nlevel)
    type(multifab) ::               phi(mla%nlevel)
    type(multifab) :: diff_rhoh_fluxdiv(mla%nlevel)

    type(multifab) ::      rhotot_fc(mla%nlevel,mla%dim)
    type(multifab) ::   rhototinv_fc(mla%nlevel,mla%dim)
    type(multifab) :: diff_mass_flux(mla%nlevel,mla%dim)

    real(kind=dp_t) :: Sbar

    type(bl_prof_timer), save :: bpt

    call build(bpt,"initial_projection")

    dm = mla%dim
    nlevs = mla%nlevel
 
    do n=1,nlevs
       call multifab_build(             conc(n),mla%la(n),nspecies   ,1)
       call multifab_build(        molarconc(n),mla%la(n),nspecies   ,1)
       call multifab_build(              eta(n),mla%la(n),1          ,1)
       call multifab_build(           lambda(n),mla%la(n),1          ,1)
       call multifab_build(            kappa(n),mla%la(n),1          ,1)
       call multifab_build(              chi(n),mla%la(n),nspecies**2,1)
       call multifab_build(             zeta(n),mla%la(n),nspecies   ,1)
       call multifab_build(          mac_rhs(n),mla%la(n),1          ,0)
       call multifab_build(            theta(n),mla%la(n),1          ,0)
       call multifab_build(              phi(n),mla%la(n),1          ,1)
       call multifab_build(diff_rhoh_fluxdiv(n),mla%la(n),1          ,0)
       do i=1,dm
          call multifab_build_edge(      rhotot_fc(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(   rhototinv_fc(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge( diff_mass_flux(n,i),mla%la(n),nspecies,0,i)
       end do       
    end do

    call convert_rhoc_to_c(mla,rho,rhotot,conc,.true.)
    call fill_c_ghost_cells(mla,conc,dx,the_bc_tower)
    call convert_conc_to_molarconc(mla,conc,molarconc,.true.)

    call ideal_mixture_transport_wrapper(mla,rhotot,Temp,p0,conc,molarconc, &
                                         eta,lambda,kappa,chi,zeta)

    ! compute diffusive mass fluxes
    ! this computes "-F = rho W chi [Gamma grad x... ]"
    call compute_mass_fluxdiv_energy(mla,rho,rhotot,molarconc,chi,zeta,gradp_baro,Temp, &
                                     diff_mass_fluxdiv,diff_mass_flux,dx,the_bc_tower)

    call compute_rhoh_fluxdiv(mla,lambda,Temp,diff_mass_flux,rhotot,diff_rhoh_fluxdiv, &
                              dx,0.d0,the_bc_tower)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! build rhs = div(v^init) - S^0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute -S^0, then change compute_div() to increment
    call compute_S_theta(mla,mac_rhs,theta,diff_mass_fluxdiv, &
                         diff_rhoh_fluxdiv,conc,Temp,rhotot)
    
    ! split S^n, theta^n, and Scorr into average and perturbational pieces
    do n=1,nlevs
       Sbar = multifab_sum_c(mac_rhs(n),1,1) / total_volume
       call multifab_sub_sub_s_c(mac_rhs(n),1,Sbar,1,0)
    end do

    ! change this to increment mac_rhs by div(v^init)
    call compute_div(mla,umac,mac_rhs,dx,1,1,1,increment_in=.true.)

    ! average rhotot to faces
    call average_cc_to_face(nlevs,rhotot,rhotot_fc,1,scal_bc_comp,1,the_bc_tower%bc_tower_array)

    ! compute (1/rhotot) on faces
    do n=1,nlevs
       do i=1,dm
          call setval(rhototinv_fc(n,i),1.d0,all=.true.)
          call multifab_div_div_c(rhototinv_fc(n,i),1,rhotot_fc(n,i),1,1,0)
       end do
    end do

    ! solve div (1/rhotot) grad phi = div(v^init) - S^0
    ! solve to completion, i.e., use the 'full' solver
    do n=1,nlevs
       call setval(phi(n),0.d0)
    end do
    call macproject(mla,phi,umac,rhototinv_fc,mac_rhs,dx,the_bc_tower,.true.)

    ! v^0 = v^init - (1/rho^0) grad phi
    call subtract_weighted_gradp(mla,umac,rhototinv_fc,phi,dx,the_bc_tower)

    ! fill ghost cells
    do n=1,nlevs
       do i=1,dm
          ! set normal velocity on physical domain boundaries
          call multifab_physbc_domainvel(umac(n,i),vel_bc_comp+i-1, &
               the_bc_tower%bc_tower_array(n), &
               dx(n,:))
          ! set transverse velocity behind physical boundaries
          call multifab_physbc_macvel(umac(n,i),vel_bc_comp+i-1, &
               the_bc_tower%bc_tower_array(n), &
               dx(n,:))
          ! fill periodic and interior ghost cells
          call multifab_fill_boundary(umac(n,i))
       end do
    end do

    do n=1,nlevs
       call multifab_destroy(conc(n))
       call multifab_destroy(molarconc(n))
       call multifab_destroy(eta(n))
       call multifab_destroy(lambda(n))
       call multifab_destroy(kappa(n))
       call multifab_destroy(chi(n))
       call multifab_destroy(zeta(n))
       call multifab_destroy(mac_rhs(n))
       call multifab_destroy(theta(n))
       call multifab_destroy(phi(n))
       call multifab_destroy(diff_rhoh_fluxdiv(n))
       do i=1,dm
          call multifab_destroy(rhotot_fc(n,i))
          call multifab_destroy(rhototinv_fc(n,i))
          call multifab_destroy(diff_mass_flux(n,i))
       end do
    end do

    call destroy(bpt)

  end subroutine initial_projection

end module initial_projection_module
