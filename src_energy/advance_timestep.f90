module advance_timestep_module

  use ml_layout_module
  use define_bc_module
  use bc_module
  use convert_stag_module
  use convert_rhoc_to_c_module
  use convert_rhoh_to_h_module
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

  subroutine advance_timestep(mla,umac_old,umac_new,rho_old,rho_new, &
                              rhotot_old,rhotot_new,rhoh_old,rhoh_new, &
                              p0_old,p0_new,gradp_baro,Temp_old,Temp_new, &
                              dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: umac_old(:,:)
    type(multifab) , intent(inout) :: umac_new(:,:)
    type(multifab) , intent(inout) :: rho_old(:)
    type(multifab) , intent(inout) :: rho_new(:)
    type(multifab) , intent(in   ) :: rhotot_old(:)
    type(multifab) , intent(in   ) :: rhotot_new(:)
    type(multifab) , intent(inout) :: rhoh_old(:)
    type(multifab) , intent(inout) :: rhoh_new(:)
    real(kind=dp_t), intent(in   ) :: p0_old
    real(kind=dp_t), intent(in   ) :: p0_new
    type(multifab) , intent(in   ) :: gradp_baro(:,:)
    type(multifab) , intent(in   ) :: Temp_old(:)
    type(multifab) , intent(in   ) :: Temp_new(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: n,nlevs,i,dm,comp
    integer :: n_cell

    ! concentrations and mole fractions
    type(multifab) :: conc_old(mla%nlevel)
    type(multifab) :: molefrac_old(mla%nlevel)

    ! viscosity
    type(multifab) :: eta_old(mla%nlevel)

    ! bulk viscosity
    type(multifab) :: kappa_old(mla%nlevel)

    ! thermal diffusivity
    type(multifab) :: lambda_old(mla%nlevel)

    ! diffusion matrix
    type(multifab) :: chi_old(mla%nlevel)

    ! thermodiffusion coefficients
    type(multifab) :: zeta_old(mla%nlevel)

    ! enthalpy
    type(multifab) :: enth_old(mla%nlevel)

    ! store F^n, div(F^n), and div(Q)^n + sum(div(hk*Fk))^n + rho_old*Hext^n
    type(multifab) :: mass_flux_old(mla%nlevel,mla%dim)
    type(multifab) :: mass_fluxdiv_old(mla%nlevel)
    type(multifab) :: rhoh_fluxdiv_old(mla%nlevel)

    ! (delta_S,delta_alpha)^n
    type(multifab) :: delta_S_old(mla%nlevel)
    type(multifab) :: delta_alpha_old(mla%nlevel)

    ! (rhotot,rho,rhoh)^n on faces
    type(multifab) :: rhotot_fc_old(mla%nlevel,mla%dim)
    type(multifab) :: rho_fc_old(mla%nlevel,mla%dim)
    type(multifab) :: rhoh_fc_old(mla%nlevel,mla%dim)

    ! Scorr and delta_Scorr
    type(multifab) :: Scorr(mla%nlevel)
    type(multifab) :: delta_Scorr(mla%nlevel)

    ! pressure from the EOS
    type(multifab) :: Peos(mla%nlevel)

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
       call multifab_build(conc_old(n)    ,mla%la(n),nspecies   ,1)
       call multifab_build(molefrac_old(n),mla%la(n),nspecies   ,1)
       call multifab_build(eta_old(n)     ,mla%la(n),1          ,1)
       call multifab_build(kappa_old(n)   ,mla%la(n),1          ,1)
       call multifab_build(lambda_old(n)  ,mla%la(n),1          ,1)
       call multifab_build(chi_old(n)     ,mla%la(n),nspecies**2,1)
       call multifab_build(zeta_old(n)    ,mla%la(n),nspecies   ,1)
       call multifab_build(enth_old(n)    ,mla%la(n),1          ,1)
       do i=1,dm
          call multifab_build_edge(mass_flux_old(n,i),mla%la(n),nspecies,0,i)
       end do
       call multifab_build(mass_fluxdiv_old(n),mla%la(n),nspecies,0)
       call multifab_build(rhoh_fluxdiv_old(n),mla%la(n),1       ,0)
       call multifab_build(delta_S_old(n)     ,mla%la(n),1       ,0)
       call multifab_build(delta_alpha_old(n) ,mla%la(n),1       ,0)
       do i=1,dm
          call multifab_build_edge(rhotot_fc_old(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(rho_fc_old(n,i)   ,mla%la(n),nspecies,0,i)
          call multifab_build_edge(rhoh_fc_old(n,i)  ,mla%la(n),1       ,0,i)
       end do
       call multifab_build(Scorr(n)      ,mla%la(n),1,0)
       call multifab_build(delta_Scorr(n),mla%la(n),1,0)
       call multifab_build(Peos(n)       ,mla%la(n),1,0)
    end do

    ! compute rhotot^n on faces
    call average_cc_to_face(nlevs,rhotot_old,rhotot_fc_old,1,scal_bc_comp,1, &
                            the_bc_tower%bc_tower_array)

    ! compute rho_i^n on faces
    ! first compute c^n and fill ghost cells
    call convert_rhoc_to_c(mla,rho_old,rhotot_old,conc_old,.true.)
    call fill_c_ghost_cells(mla,conc_old,dx,the_bc_tower)

    ! average c^n to faces
    call average_cc_to_face(nlevs,conc_old,rho_fc_old,1,c_bc_comp,nspecies, &
                            the_bc_tower%bc_tower_array)

    ! multiply c^n on faces by rhotot^n on faces
    do n=1,nlevs
       do i=1,dm
          do comp=1,nspecies
             call multifab_mult_mult_c(rho_fc_old(n,i),comp,rhotot_fc_old(n,i),1,1,0)
          end do
       end do
    end do

    ! compute rhoh^n on faces
    ! first compute h^n and fill ghost cells
    call convert_rhoh_to_h(mla,rhoh_old,rhotot_old,enth_old,.true.)
    call fill_h_ghost_cells(mla,enth_old,dx,the_bc_tower)

    ! average h^n to faces
    call average_cc_to_face(nlevs,enth_old,rhoh_fc_old,1,h_bc_comp,1, &
                            the_bc_tower%bc_tower_array)

    ! multiply h^n on faces by rhotot^n on faces
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_c(rhoh_fc_old(n,i),1,rhotot_fc_old(n,i),1,1,0)
       end do
    end do

    ! compute mole fractions in VALID + GHOST regions
    call convert_conc_to_molefrac(mla,conc_old,molefrac_old,.true.)

    ! compute t^n transport properties
    call ideal_mixture_transport_wrapper(mla,rhotot_old,Temp_old,p0_old,conc_old, &
                                         molefrac_old,eta_old,lambda_old,kappa_old, &
                                         chi_old,zeta_old)

    ! compute mass_flux_old = F^n and mass_fluxdiv_old = div(F^n)
    call mass_fluxdiv_energy(mla,rho_old,rhotot_old,molefrac_old,chi_old,zeta_old, &
                             gradp_baro,Temp_old,mass_fluxdiv_old, &
                             mass_flux_old,dx,the_bc_tower)

    ! compute rhoh_fluxdiv_old = div(Q)^n + sum(div(hk*Fk))^n + rho_old*Hext^n
    call rhoh_fluxdiv_energy(mla,lambda_old,Temp_old,mass_flux_old,rhotot_old, &
                             rhoh_fluxdiv_old,dx,0.d0,the_bc_tower)

    ! compute S^n and alpha^n (store them in delta_S_old and delta_alpha_old)
    call compute_S_alpha(mla,delta_S_old,delta_alpha_old,mass_fluxdiv_old, &
                         rhoh_fluxdiv_old,conc_old,Temp_old,rhotot_old,p0_old)

    ! compute P_eos^n
    call compute_p(mla,rhotot_old,Temp_old,conc_old,Peos)

    ! Scorr = Scorr + alpha^n * (Peos^n - P0^n)/dt
    do n=1,nlevs
       call multifab_sub_sub_s_c(Peos(n),1,p0_old,1,0)
       call multifab_mult_mult_s_c(Peos(n),1,1.d0/dt,1,0)
       call multifab_mult_mult_c(Peos(n),1,delta_alpha_old(n),1,1,0)
       call multifab_copy_c(Scorr(n),1,Peos(n),1,1,0)
    end do

    ! split S^n, alpha^n, and Scorr into average and perturbational pieces
    do n=1,nlevs
       Sbar_old     = multifab_sum_c(delta_S_old(n)    ,1,1) / dble(n_cell)
       alphabar_old = multifab_sum_c(delta_alpha_old(n),1,1) / dble(n_cell)
       Scorrbar     = multifab_sum_c(Scorr(n),1,1) / dble(n_cell)
       call multifab_sub_sub_s_c(delta_S_old(n)    ,1,Sbar_old    ,1,0)
       call multifab_sub_sub_s_c(delta_alpha_old(n),1,alphabar_old,1,0)
       call multifab_copy_c(delta_Scorr(n),1,Scorr(n),1,1,0)
       call multifab_sub_sub_s_c(delta_Scorr(n),1,Scorrbar,1,0)
    end do

    ! new state begins as a copy of old state
    do n=1,nlevs
       call multifab_copy_c(rho_new(n)   ,1,rho_old(n)   ,1,nspecies,rho_new(n)%ng)
       call multifab_copy_c(rhotot_new(n),1,rhotot_old(n),1,1       ,rhotot_new(n)%ng)
       call multifab_copy_c(rhoh_new(n)  ,1,rhoh_old(n)  ,1,1       ,rhoh_new(n)%ng)
       call multifab_copy_c(Temp_new(n)  ,1,Temp_old(n)  ,1,1       ,Temp_new(n)%ng)

       do i=1,dm
          call multifab_copy_c(umac_new(n,i),1,umac_old(n,i),1,1,umac_new(n,i)%ng)
       end do
    end do




    do n=1,nlevs
       call multifab_destroy(conc_old(n))
       call multifab_destroy(molefrac_old(n))
       call multifab_destroy(eta_old(n))
       call multifab_destroy(kappa_old(n))
       call multifab_destroy(lambda_old(n))
       call multifab_destroy(chi_old(n))
       call multifab_destroy(zeta_old(n))
       call multifab_destroy(enth_old(n))
       do i=1,dm
          call multifab_destroy(mass_flux_old(n,i))
       end do
       call multifab_destroy(mass_fluxdiv_old(n))
       call multifab_destroy(rhoh_fluxdiv_old(n))
       call multifab_destroy(delta_S_old(n))
       call multifab_destroy(delta_alpha_old(n))
       do i=1,dm
          call multifab_destroy(rhotot_fc_old(n,i))
          call multifab_destroy(rho_fc_old(n,i))
          call multifab_destroy(rhoh_fc_old(n,i))
       end do
       call multifab_destroy(Scorr(n))
       call multifab_destroy(delta_Scorr(n))
    end do

  end subroutine advance_timestep

end module advance_timestep_module
