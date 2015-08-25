module initial_projection_module

  use ml_layout_module
  use define_bc_module
  use bc_module
  use convert_stag_module
  use convert_rhoc_to_c_module
  use eos_model_wrapper_module
  use mass_fluxdiv_energy_module
  use rhoh_fluxdiv_energy_module
  use div_and_grad_module
  use macproject_module
  use multifab_physbc_stag_module
  use probin_common_module, only: total_volume
  use probin_multispecies_module, only: nspecies

  implicit none

  private

  public :: initial_projection

contains

  subroutine initial_projection(mla,umac,rho,rhotot,rhoh,p0, &
                                gradp_baro,Temp,dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(inout) :: rho(:)
    type(multifab) , intent(in   ) :: rhotot(:)
    type(multifab) , intent(in   ) :: rhoh(:)
    real(kind=dp_t), intent(in   ) :: p0
    type(multifab) , intent(in   ) :: gradp_baro(:,:)
    type(multifab) , intent(in   ) :: Temp(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: n,nlevs,i,dm

    ! temporary storage for concentrations and mole fractions
    type(multifab) :: conc(mla%nlevel)
    type(multifab) :: molefrac(mla%nlevel)

    ! viscosity
    type(multifab) :: eta(mla%nlevel)

    ! bulk viscosity
    type(multifab) :: kappa(mla%nlevel)

    ! thermal diffusivity
    type(multifab) :: lambda(mla%nlevel)

    ! diffusion matrix
    type(multifab) :: chi(mla%nlevel)

    ! thermodiffusion coefficients
    type(multifab) :: zeta(mla%nlevel)

    ! store F^0, div(F^0), and div(Q)^0 + sum(div(hk*Fk))^0 + rho*Hext^0
    type(multifab) :: mass_flux(mla%nlevel,mla%dim)
    type(multifab) :: mass_fluxdiv(mla%nlevel)
    type(multifab) :: rhoh_fluxdiv(mla%nlevel)

    ! delta_S^0 and delta_theta^0
    type(multifab) :: delta_S(mla%nlevel)
    type(multifab) :: delta_theta(mla%nlevel)

    ! S_proj
    type(multifab) :: Sproj(mla%nlevel)

    ! rhotot^0 and 1/rhotot^0 on faces
    type(multifab) :: rhotot_fc(mla%nlevel,mla%dim)
    type(multifab) :: rhototinv_fc(mla%nlevel,mla%dim)

    ! solution of the pressure-projection solve
    type(multifab) :: phi(mla%nlevel)

    real(kind=dp_t) :: Sbar
    real(kind=dp_t) :: thetabar

    nlevs = mla%nlevel
    dm = mla%dim

    do n=1,nlevs
       call multifab_build(conc(n)    ,mla%la(n),nspecies   ,1)
       call multifab_build(molefrac(n),mla%la(n),nspecies   ,1)
       call multifab_build(eta(n)     ,mla%la(n),1          ,1)
       call multifab_build(kappa(n)   ,mla%la(n),1          ,1)
       call multifab_build(lambda(n)  ,mla%la(n),1          ,1)
       call multifab_build(chi(n)     ,mla%la(n),nspecies**2,1)
       call multifab_build(zeta(n)    ,mla%la(n),nspecies   ,1)
       do i=1,dm
          call multifab_build_edge(mass_flux(n,i),mla%la(n),nspecies,0,i)
       end do
       call multifab_build(mass_fluxdiv(n),mla%la(n),nspecies,0)
       call multifab_build(rhoh_fluxdiv(n),mla%la(n),1       ,0)
       call multifab_build(delta_S(n)     ,mla%la(n),1       ,0)
       call multifab_build(delta_theta(n) ,mla%la(n),1       ,0)
       call multifab_build(Sproj(n)       ,mla%la(n),1       ,0)
       do i=1,dm
          call multifab_build_edge(rhotot_fc(n,i)   ,mla%la(n),1,0,i)
          call multifab_build_edge(rhototinv_fc(n,i),mla%la(n),1,0,i)
       end do
       call multifab_build(phi(n),mla%la(n),1,1)
    end do

    ! compute rhotot^0 on faces
    call average_cc_to_face(nlevs,rhotot,rhotot_fc,1,scal_bc_comp,1, &
                            the_bc_tower%bc_tower_array)

    ! compute (1/rhotot)^0 on faces
    do n=1,nlevs
       do i=1,dm
          call setval(rhototinv_fc(n,i),1.d0,all=.true.)
          call multifab_div_div_c(rhototinv_fc(n,i),1,rhotot_fc(n,i),1,1,0)
       end do
    end do

    ! compute mass fractions in valid region and then fill ghost cells
    call convert_rhoc_to_c(mla,rho,rhotot,conc,.true.)
    call fill_c_ghost_cells(mla,conc,dx,the_bc_tower)

    ! compute mole fractions in VALID + GHOST regions
    call convert_conc_to_molefrac(mla,conc,molefrac,.true.)

    ! compute t^0 transport properties
    call ideal_mixture_transport_wrapper(mla,rhotot,Temp,p0,conc,molefrac, &
                                         eta,lambda,kappa,chi,zeta)

    ! compute mass_flux = F^0 and mass_fluxdiv = div(F^0)
    call mass_fluxdiv_energy(mla,rho,rhotot,molefrac,chi,zeta, &
                             gradp_baro,Temp,mass_fluxdiv, &
                             mass_flux,dx,the_bc_tower)

    ! compute rhoh_fluxdiv = div(Q)^0 + sum(div(hk*Fk))^0 + rho*Hext^0
    call rhoh_fluxdiv_energy(mla,lambda,Temp,mass_flux,rhotot, &
                             rhoh_fluxdiv,dx,0.d0,the_bc_tower)

    ! compute S^0 and theta^0 (store them in delta_S and delta_theta)
    call compute_S_theta(mla,delta_S,delta_theta,mass_fluxdiv, &
                         rhoh_fluxdiv,conc,Temp,rhotot)

    ! split S and theta into average and perturbational pieces
    ! S = Sbar + delta_S
    ! theta = thetabar + delta_theta
    do n=1,nlevs
       Sbar     = multifab_sum_c(delta_S(n)    ,1,1) / total_volume
       thetabar = multifab_sum_c(delta_theta(n),1,1) / total_volume
       call multifab_sub_sub_s_c(delta_S(n)    ,1,Sbar    ,1,0)
       call multifab_sub_sub_s_c(delta_theta(n),1,thetabar,1,0)
    end do

    ! compute Sproj = delta_S - delta_theta * Sbar /thetabar
    do n=1,nlevs
       call multifab_copy_c(Sproj(n),1,delta_theta(n),1,1,0)
       call multifab_mult_mult_s_c(Sproj(n),1,-Sbar/thetabar,1,0)
       call multifab_plus_plus_c(Sproj(n),1,delta_S(n),1,1,0)
    end do

    ! build rhs for projection, div(v^init) - Sproj
    ! first multiply Sproj by -1
    do n=1,nlevs
       call multifab_mult_mult_s_c(Sproj(n),1,-1.d0,1,0)
    end do

    ! add div(v^init) to Sproj
    call compute_div(mla,umac,Sproj,dx,1,1,1,increment_in=.true.)

    ! solve div (1/rhotot^0) grad phi = div(v^init) - Sproj^0
    ! solve to completion, i.e., use the 'full' solver
    call macproject(mla,phi,umac,rhototinv_fc,Sproj,dx,the_bc_tower,.true.)

    ! v^0 = v^init - (1/rho^0) grad phi
    call subtract_weighted_gradp(mla,umac,rhototinv_fc,phi,dx,the_bc_tower)
    
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
       call multifab_destroy(molefrac(n))
       call multifab_destroy(eta(n))
       call multifab_destroy(kappa(n))
       call multifab_destroy(lambda(n))
       call multifab_destroy(chi(n))
       call multifab_destroy(zeta(n))
       do i=1,dm
          call multifab_destroy(mass_flux(n,i))
       end do
       call multifab_destroy(mass_fluxdiv(n))
       call multifab_destroy(rhoh_fluxdiv(n))
       call multifab_destroy(delta_S(n))
       call multifab_destroy(delta_theta(n))
       call multifab_destroy(Sproj(n))
       do i=1,dm
          call multifab_destroy(rhotot_fc(n,i))
          call multifab_destroy(rhototinv_fc(n,i))
       end do
       call multifab_destroy(phi(n))
    end do

  end subroutine initial_projection

end module initial_projection_module
