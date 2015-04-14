module advance_timestep_module

  use ml_layout_module
  use define_bc_module
  use bc_module
  use convert_stag_module
  use convert_rhoc_to_c_module
  use energy_eos_wrapper_module
  use mass_fluxdiv_energy_module
  use rhoh_fluxdiv_energy_module
  use div_and_grad_module
  use macproject_module
  use probin_common_module, only: n_cells
  use probin_multispecies_module, only: nspecies

  implicit none

  private

  public :: advance_timestep

contains

  subroutine advance_timestep(mla,umac_old,rho_old,rhotot_old,rhoh_old,p0_old, &
                              gradp_baro,Temp_old,dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: umac_old(:,:)
    type(multifab) , intent(inout) :: rho_old(:)
    type(multifab) , intent(in   ) :: rhotot_old(:)
    type(multifab) , intent(in   ) :: rhoh_old(:)
    real(kind=dp_t), intent(in   ) :: p0_old
    type(multifab) , intent(in   ) :: gradp_baro(:,:)
    type(multifab) , intent(in   ) :: Temp_old(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: n,nlevs,i,dm
    integer :: n_cell

    ! concentrations and mole fractions
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

    ! store F^n, div(F^n), and div(Q)^n + sum(div(hk*Fk))^n + rho_old*Hext^n
    type(multifab) :: mass_flux_old(mla%nlevel,mla%dim)
    type(multifab) :: mass_fluxdiv_old(mla%nlevel)
    type(multifab) :: rhoh_fluxdiv_old(mla%nlevel)

    ! (delta_S,delta_alpha)^n
    type(multifab) :: delta_S_old(mla%nlevel)
    type(multifab) :: delta_alpha_old(mla%nlevel)

    ! (rhotot,rho,rhoh)^n on faces
    type(multifab) :: rhotot_old_fc(mla%nlevel,mla%dim)
    type(multifab) :: rho_old_fc(mla%nlevel,mla%dim)
    type(multifab) :: rhoh_old_fc(mla%nlevel,mla%dim)

    ! Scorr and delta_Scorr
    type(multifab) :: Scorr(mla%nlevel)
    type(multifab) :: delta_Scorr(mla%nlevel)

    real(kind=dp_t) :: Sbar_old, Sbar_new
    real(kind=dp_t) :: alphabar_old, alphabar_new
    real(kind=dp_t) :: Scorrbar

    nlevs = mla%nlevel
    dm = mla%dim

    if (dm .eq. 2) then
       n_cell = n_cells(1)*n_cells(2)
    else
       n_cell = n_cells(1)*n_cells(2)*n_cells(3)
    end if

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
       call multifab_build(mass_fluxdiv_old(n),mla%la(n),nspecies,0)
       call multifab_build(rhoh_fluxdiv_old(n),mla%la(n),1       ,0)
       call multifab_build(delta_S_old(n)     ,mla%la(n),1       ,0)
       call multifab_build(delta_alpha_old(n) ,mla%la(n),1       ,0)
       do i=1,dm
          call multifab_build_edge(rhotot_old_fc(n,i)   ,mla%la(n),1,0,i)
       end do
    end do

    ! compute rhotot_old^n on faces
    call average_cc_to_face(nlevs,rhotot_old,rhotot_old_fc,1,scal_bc_comp,1, &
                            the_bc_tower%bc_tower_array)

    ! compute rhoh^n on faces


    ! compute mass fractions in valid region and then fill ghost cells
    call convert_rho_oldc_to_c(mla,rho_old,rhotot_old,conc,.true.)
    call fill_c_ghost_cells(mla,conc,dx,the_bc_tower)

    ! compute mole fractions in VALID + GHOST regions
    call convert_conc_to_molefrac(mla,conc,molefrac,.true.)

    ! compute t^n transport properties
    call ideal_mixture_transport_wrapper(mla,rhotot_old,Temp_old,p0_old,conc,molefrac, &
                                         eta,lambda,kappa,chi,zeta)

    ! compute mass_flux = F^n and mass_fluxdiv_old = div(F^n)
    call mass_fluxdiv_old_energy(mla,rho_old,rhotot_old,molefrac,chi,zeta, &
                             gradp_bar_oldo,Temp_old,mass_fluxdiv_old, &
                             mass_flux_old,dx,the_bc_tower)

    ! compute rhoh_fluxdiv_old = div(Q)^n + sum(div(hk*Fk))^n + rho_old*Hext^n
    call rhoh_fluxdiv_old_energy(mla,lambda,Temp_old,mass_flux_old,rhotot_old, &
                             rhoh_fluxdiv_old,dx,0.d0,the_bc_tower)

    ! compute S^n and alpha^n (store them in delta_S_old and delta_alpha_old)
    call compute_S_alpha(mla,delta_S_old,delta_alpha_old,mass_fluxdiv_old, &
                         rhoh_fluxdiv_old,conc,Temp_old,rhotot_old,p0_old)

    ! split S and alpha into average and perturbational pieces
    ! S = Sbar_old + delta_S_old
    ! alpha = alphabar_old + delta_alpha_old
    do n=1,nlevs
       Sbar_old     = multifab_sum_c(delta_S_old(n)    ,1,1) / dble(n_cell)
       alphabar_old = multifab_sum_c(delta_alpha_old(n),1,1) / dble(n_cell)
       call multifab_sub_sub_s_c(delta_S_old(n)    ,1,Sbar_old    ,1,0)
       call multifab_sub_sub_s_c(delta_alpha_old(n),1,alphabar_old,1,0)
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
       call multifab_destroy(mass_fluxdiv_old(n))
       call multifab_destroy(rhoh_fluxdiv_old(n))
       call multifab_destroy(delta_S_old(n))
       call multifab_destroy(delta_alpha_old(n))
       do i=1,dm
          call multifab_destroy(rhotot_old_fc(n,i))
       end do
    end do

  end subroutine advance_timestep

end module advance_timestep_module
