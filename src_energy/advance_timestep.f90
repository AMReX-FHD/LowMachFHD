module advance_timestep_module

  use ml_layout_module
  use bndry_reg_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use convert_stag_module
  use convert_rhoc_to_c_module
  use convert_rhoh_to_h_module
  use energy_eos_wrapper_module
  use mass_flux_utilities_module
  use mass_fluxdiv_energy_module
  use rhoh_fluxdiv_energy_module
  use div_and_grad_module
  use macproject_module
  use mk_advective_s_fluxdiv_module
  use ml_solve_module
  use probin_common_module, only: n_cells
  use probin_multispecies_module, only: nspecies
  use energy_EOS_module, only: dpdt_iters, deltaT_iters

  use fabio_module

  implicit none

  private

  public :: advance_timestep

contains

  subroutine advance_timestep(mla,umac_old,umac_new,rho_old,rho_new, &
                              rhotot_old,rhotot_new,rhoh_old,rhoh_new, &
                              p0_old,p0_new,gradp_baro,Temp_old,Temp_new, &
                              dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: umac_old(:,:)
    type(multifab) , intent(inout) :: umac_new(:,:)
    type(multifab) , intent(inout) :: rho_old(:)
    type(multifab) , intent(inout) :: rho_new(:)
    type(multifab) , intent(inout) :: rhotot_old(:)
    type(multifab) , intent(inout) :: rhotot_new(:)
    type(multifab) , intent(inout) :: rhoh_old(:)
    type(multifab) , intent(inout) :: rhoh_new(:)
    real(kind=dp_t), intent(in   ) :: p0_old
    real(kind=dp_t), intent(inout) :: p0_new
    type(multifab) , intent(in   ) :: gradp_baro(:,:)
    type(multifab) , intent(in   ) :: Temp_old(:)
    type(multifab) , intent(inout) :: Temp_new(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: n,nlevs,i,dm,comp
    integer :: n_cell
    integer :: k,l

    ! concentration
    type(multifab) :: conc_old(mla%nlevel)
    type(multifab) :: conc_new(mla%nlevel)

    ! mole fractions
    type(multifab) :: molefrac_old(mla%nlevel)
    type(multifab) :: molefrac_new(mla%nlevel)

    ! viscosity
    type(multifab) :: eta_old(mla%nlevel)
    type(multifab) :: eta_new(mla%nlevel)

    ! bulk viscosity
    type(multifab) :: kappa_old(mla%nlevel)
    type(multifab) :: kappa_new(mla%nlevel)

    ! thermal diffusivity
    type(multifab) :: lambda_old(mla%nlevel)
    type(multifab) :: lambda_new(mla%nlevel)

    ! diffusion matrix
    type(multifab) :: chi_old(mla%nlevel)
    type(multifab) :: chi_new(mla%nlevel)

    ! thermodiffusion coefficients
    type(multifab) :: zeta_old(mla%nlevel)
    type(multifab) :: zeta_new(mla%nlevel)

    ! enthalpy
    type(multifab) :: enth_old(mla%nlevel)
    type(multifab) :: enth_new(mla%nlevel)

    ! store F
    type(multifab) :: mass_flux_old(mla%nlevel,mla%dim)
    type(multifab) :: mass_flux_new(mla%nlevel,mla%dim)

    ! div(F)
    type(multifab) :: mass_fluxdiv_old(mla%nlevel)
    type(multifab) :: mass_fluxdiv_new(mla%nlevel)

    ! -div(rho_i*v) + div(F)
    type(multifab) :: mass_update_old(mla%nlevel)
    type(multifab) :: mass_update_new(mla%nlevel)

    ! div(Q) + sum(div(hk*Fk)) + rho*Hext
    type(multifab) :: rhoh_fluxdiv_old(mla%nlevel)
    type(multifab) :: rhoh_fluxdiv_new(mla%nlevel)

    ! -div(rhoh*v) + p0_update + div(Q) + sum(div(hk*Fk)) + rho*Hext
    type(multifab) :: rhoh_update_old(mla%nlevel)
    type(multifab) :: rhoh_update_new(mla%nlevel)

    ! delta_S
    type(multifab) :: delta_S_old(mla%nlevel)
    type(multifab) :: delta_S_new(mla%nlevel)

    ! delta_alpha
    type(multifab) :: delta_alpha_old(mla%nlevel)
    type(multifab) :: delta_alpha_new(mla%nlevel)

    ! rhotot on faces
    type(multifab) :: rhotot_fc_old(mla%nlevel,mla%dim)
    type(multifab) :: rhotot_fc_new(mla%nlevel,mla%dim)

    ! rho on faces
    type(multifab) :: rho_fc_old(mla%nlevel,mla%dim)
    type(multifab) :: rho_fc_new(mla%nlevel,mla%dim)

    ! rhoh on faces
    type(multifab) :: rhoh_fc_old(mla%nlevel,mla%dim)
    type(multifab) :: rhoh_fc_new(mla%nlevel,mla%dim)

    ! Scorr and delta_Scorr
    type(multifab) :: Scorr(mla%nlevel)
    type(multifab) :: delta_Scorr(mla%nlevel)

    ! pressure from the EOS
    type(multifab) :: Peos(mla%nlevel)

    ! solution and RHS for deltaT solve
    type(multifab) :: deltaT(mla%nlevel)
    type(multifab) :: deltaT_rhs(mla%nlevel)

    ! coefficients for deltaT (energy) solve
    type(multifab) :: cc_solver_alpha(mla%nlevel)
    type(multifab) :: cc_solver_beta(mla%nlevel,mla%dim)

    ! for energy implicit solve
    ! doesn't actually do anything for single-level solves
    type(bndry_reg) :: fine_flx(2:mla%nlevel)

    real(kind=dp_t) :: Sbar_old, Sbar_new
    real(kind=dp_t) :: alphabar_old, alphabar_new
    real(kind=dp_t) :: Scorrbar

    real(kind=dp_t) :: p0_update_old
    real(kind=dp_t) :: p0_update_new

    nlevs = mla%nlevel
    dm = mla%dim

    if (dm .eq. 2) then
       n_cell = n_cells(1)*n_cells(2)
    else
       n_cell = n_cells(1)*n_cells(2)*n_cells(3)
    end if

    do n=1,nlevs
       call multifab_build(conc_old(n)    ,mla%la(n),nspecies   ,1)
       call multifab_build(conc_new(n)    ,mla%la(n),nspecies   ,1)
       call multifab_build(molefrac_old(n),mla%la(n),nspecies   ,1)
       call multifab_build(molefrac_new(n),mla%la(n),nspecies   ,1)
       call multifab_build(eta_old(n)     ,mla%la(n),1          ,1)
       call multifab_build(eta_new(n)     ,mla%la(n),1          ,1)
       call multifab_build(kappa_old(n)   ,mla%la(n),1          ,1)
       call multifab_build(kappa_new(n)   ,mla%la(n),1          ,1)
       call multifab_build(lambda_old(n)  ,mla%la(n),1          ,1)
       call multifab_build(lambda_new(n)  ,mla%la(n),1          ,1)
       call multifab_build(chi_old(n)     ,mla%la(n),nspecies**2,1)
       call multifab_build(chi_new(n)     ,mla%la(n),nspecies**2,1)
       call multifab_build(zeta_old(n)    ,mla%la(n),nspecies   ,1)
       call multifab_build(zeta_new(n)    ,mla%la(n),nspecies   ,1)
       call multifab_build(enth_old(n)    ,mla%la(n),1          ,1)
       call multifab_build(enth_new(n)    ,mla%la(n),1          ,1)
       do i=1,dm
          call multifab_build_edge(mass_flux_old(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge(mass_flux_new(n,i),mla%la(n),nspecies,0,i)
       end do
       call multifab_build(mass_fluxdiv_old(n),mla%la(n),nspecies,0)
       call multifab_build(mass_fluxdiv_new(n),mla%la(n),nspecies,0)
       call multifab_build(mass_update_old(n) ,mla%la(n),nspecies,0)
       call multifab_build(mass_update_new(n) ,mla%la(n),nspecies,0)
       call multifab_build(rhoh_fluxdiv_old(n),mla%la(n),1       ,0)
       call multifab_build(rhoh_fluxdiv_new(n),mla%la(n),1       ,0)
       call multifab_build(rhoh_update_old(n) ,mla%la(n),1       ,0)
       call multifab_build(rhoh_update_new(n) ,mla%la(n),1       ,0)
       call multifab_build(delta_S_old(n)     ,mla%la(n),1       ,0)
       call multifab_build(delta_S_new(n)     ,mla%la(n),1       ,0)
       call multifab_build(delta_alpha_old(n) ,mla%la(n),1       ,0)
       call multifab_build(delta_alpha_new(n) ,mla%la(n),1       ,0)
       do i=1,dm
          call multifab_build_edge(rhotot_fc_old(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(rhotot_fc_new(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(rho_fc_old(n,i)   ,mla%la(n),nspecies,0,i)
          call multifab_build_edge(rho_fc_new(n,i)   ,mla%la(n),nspecies,0,i)
          call multifab_build_edge(rhoh_fc_old(n,i)  ,mla%la(n),1       ,0,i)
          call multifab_build_edge(rhoh_fc_new(n,i)  ,mla%la(n),1       ,0,i)
       end do
       call multifab_build(Scorr(n)          ,mla%la(n),1,0)
       call multifab_build(delta_Scorr(n)    ,mla%la(n),1,0)
       call multifab_build(Peos(n)           ,mla%la(n),1,0)
       call multifab_build(deltaT(n)         ,mla%la(n),1,1)
       call multifab_build(deltaT_rhs(n)     ,mla%la(n),1,0)
       call multifab_build(cc_solver_alpha(n),mla%la(n),1,0)
       do i=1,dm
          call multifab_build_edge(cc_solver_beta(n,i),mla%la(n),1,0,i)
       end do
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

    p0_update_old = (Sbar_old + Scorrbar)/alphabar_old

    ! mass_update_old = [-div(rho_i*v) + div(F)]^n
    do n=1,nlevs
       call multifab_copy_c(mass_update_old(n),1,mass_fluxdiv_old(n),1,nspecies,0)
    end do
    call mk_advective_s_fluxdiv(mla,umac_old,rho_fc_old,mass_update_old,dx,1,nspecies)

    ! rhoh_update_old = [-div(rhoh*v) + p0_update + div(Q) + sum(div(hk*Fk)) + rho*Hext]^n
    do n=1,nlevs
       call multifab_copy_c(rhoh_update_old(n),1,rhoh_fluxdiv_old(n),1,1,0)
       call multifab_plus_plus_s_c(rhoh_update_old(n),1,p0_update_old,1,0)
    end do
    call mk_advective_s_fluxdiv(mla,umac_old,rhoh_fc_old,rhoh_update_old,dx,1,1)
    
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

    p0_new = p0_old

    do k=1,dpdt_iters
       
       if (k .ne. 1) then

          ! compute t^{n+1} transport properties
          call ideal_mixture_transport_wrapper(mla,rhotot_new,Temp_new,p0_new,conc_new, &
                                               molefrac_new,eta_new,lambda_new,kappa_new, &
                                               chi_new,zeta_new)

          ! compute mass_flux_new = F^{n+1} and mass_fluxdiv_new = div(F^{n+1})
          call mass_fluxdiv_energy(mla,rho_new,rhotot_new,molefrac_new,chi_new,zeta_new, &
                                   gradp_baro,Temp_new,mass_fluxdiv_new, &
                                   mass_flux_new,dx,the_bc_tower)

          ! compute rhoh_fluxdiv_new = div(Q)^{n+1} + sum(div(hk*Fk))^{n+1} + rho_new*Hext^{n+1}
          call rhoh_fluxdiv_energy(mla,lambda_new,Temp_new,mass_flux_new,rhotot_new, &
                                   rhoh_fluxdiv_new,dx,0.d0,the_bc_tower)

          ! compute S^{n+1} and alpha^{n+1} (store them in delta_S_new and delta_alpha_new)
          call compute_S_alpha(mla,delta_S_new,delta_alpha_new,mass_fluxdiv_new, &
                               rhoh_fluxdiv_new,conc_new,Temp_new,rhotot_new,p0_new)

          ! compute P_eos^{n+1}
          call compute_p(mla,rhotot_new,Temp_new,conc_new,Peos)

          ! Scorr = Scorr + 2 * alpha^{n+1} * (Peos^{n+1} - P0^{n+1})/dt
          do n=1,nlevs
             call multifab_sub_sub_s_c(Peos(n),1,p0_new,1,0)
             call multifab_mult_mult_s_c(Peos(n),1,2.d0/dt,1,0)
             call multifab_mult_mult_c(Peos(n),1,delta_alpha_new(n),1,1,0)
             call multifab_plus_plus_c(Scorr(n),1,Peos(n),1,1,0)
          end do

          ! split S^{n+1}, alpha^{n+1}, and Scorr into average and perturbational pieces
          do n=1,nlevs
             Sbar_new     = multifab_sum_c(delta_S_new(n)    ,1,1) / dble(n_cell)
             alphabar_new = multifab_sum_c(delta_alpha_new(n),1,1) / dble(n_cell)
             Scorrbar     = multifab_sum_c(Scorr(n),1,1) / dble(n_cell)
             call multifab_sub_sub_s_c(delta_S_new(n)    ,1,Sbar_new    ,1,0)
             call multifab_sub_sub_s_c(delta_alpha_new(n),1,alphabar_new,1,0)
             call multifab_copy_c(delta_Scorr(n),1,Scorr(n),1,1,0)
             call multifab_sub_sub_s_c(delta_Scorr(n),1,Scorrbar,1,0)
          end do
          
       else

          ! if this is the first iteration, save a few FLOPs by copying old values
          ! and also we do not need to compute the volume discrepancy correction

          do n=1,nlevs
             call multifab_copy_c(eta_new(n)   ,1,eta_old(n)   ,1,1,1)
             call multifab_copy_c(lambda_new(n),1,lambda_old(n),1,1,1)
             call multifab_copy_c(kappa_new(n) ,1,kappa_old(n) ,1,1,1)
             call multifab_copy_c(chi_new(n)   ,1,chi_old(n)   ,1,1,1)
             call multifab_copy_c(zeta_new(n)  ,1,zeta_old(n)  ,1,1,1)

             do i=1,dm
                call multifab_copy_c(mass_flux_new(n,i),1,mass_flux_old(n,i),1,nspecies,0)
             end do
             call multifab_copy_c(mass_fluxdiv_new(n),1,mass_fluxdiv_old(n),1,nspecies,0)
             call multifab_copy_c(rhoh_fluxdiv_new(n),1,rhoh_fluxdiv_old(n),1,1       ,0)

             call multifab_copy_c(delta_S_new(n)    ,1,delta_S_old(n)    ,1,1,0)
             call multifab_copy_c(delta_alpha_new(n),1,delta_alpha_old(n),1,1,0)
          end do

          Sbar_new = Sbar_old
          alphabar_new = alphabar_old

       end if

       p0_update_new = (Sbar_new + Scorrbar)/alphabar_new

       ! update pressure
       p0_new = p0_old + 0.5d0*dt*(p0_update_old + p0_update_new)

       ! mass_update_new = [-div(rho_i*v) + div(F)]^{n+1}
       do n=1,nlevs
          call multifab_copy_c(mass_update_new(n),1,mass_fluxdiv_new(n),1,nspecies,0)
       end do
       call mk_advective_s_fluxdiv(mla,umac_new,rho_fc_new,mass_update_new,dx,1,nspecies)

       ! rhoh_update_new = [-div(rhoh*v) + p0_update + div(Q) + sum(div(hk*Fk)) + rho*Hext]^{n+1}
       do n=1,nlevs
          call multifab_copy_c(rhoh_update_new(n),1,rhoh_fluxdiv_new(n),1,1,0)
          call multifab_plus_plus_s_c(rhoh_update_new(n),1,p0_update_new,1,0)
       end do
       call mk_advective_s_fluxdiv(mla,umac_new,rhoh_fc_new,rhoh_update_new,dx,1,1)

       ! update densities
       do n=1,nlevs
          call multifab_saxpy_5(rho_new(n),1.d0,rho_old(n),0.5d0*dt,mass_update_old(n))
          call multifab_saxpy_3(rho_new(n),0.5d0*dt,mass_update_new(n))
       end do
       ! compute rhotot_new = sum(rho_new) in VALID REGION ONLY
       call compute_rhotot(mla,rho_new,rhotot_new)

       ! fill ghost cells
       do n=1,nlevs
          call multifab_fill_boundary(rhotot_new(n))
          call multifab_physbc(rhotot_new(n),1,scal_bc_comp,1, &
                               the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
       end do

       ! compute mass fractions in valid region and then fill ghost cells
       ! then convert back to densities to fill ghost cells
       call convert_rhoc_to_c(mla,rho_new,rhotot_new,conc_new,.true.)
       call fill_c_ghost_cells(mla,conc_new,dx,the_bc_tower)
       call convert_rhoc_to_c(mla,rho_new,rhotot_new,conc_new,.false.)

       ! compute mole fractions in VALID + GHOST regions
       call convert_conc_to_molefrac(mla,conc_new,molefrac_new,.true.)

       do l=1,deltaT_iters

          ! build RHS for deltaT solve
          do n=1,nlevs
             call multifab_copy_c(deltaT_rhs(n),1,rhoh_old(n),1,1,0)
             call multifab_sub_sub_c(deltaT_rhs(n),1,rhoh_new(n),1,1,0)
             call multifab_mult_mult_s_c(deltaT_rhs(n),1,1.d0/dt,1,0)
             call multifab_saxpy_3(deltaT_rhs(n),0.5d0,rhoh_update_old(n))
             call multifab_saxpy_3(deltaT_rhs(n),0.5d0,rhoh_update_new(n))
          end do

          ! cc_solver_alpha = rho^{n+1} c_p^{n+1,l} / dt
          call compute_cp(mla,cc_solver_alpha,conc_new,Temp_new)
          do n=1,nlevs
             call multifab_mult_mult_c(cc_solver_alpha(n),1,rhotot_new(n),1,1,0)
             call multifab_mult_mult_s_c(cc_solver_alpha(n),1,1.d0/dt,1,0)
          end do

          ! cc_solver_beta = (1/2) lambda^{n+1,l}
          call average_cc_to_face(nlevs,lambda_new,cc_solver_beta,1,tran_bc_comp,1, &
                                  the_bc_tower%bc_tower_array)          
          do n=1,nlevs
             do i=1,dm
                call multifab_mult_mult_s_c(cc_solver_beta(n,i),1,0.5d0,1,0)
             end do
          end do

          ! solve for deltaT
          do n = 2,nlevs
             call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
          end do

          do n=1,nlevs
             call setval(deltaT(n),0.d0,all=.true.)
          end do

          call ml_cc_solve(mla,deltaT_rhs,deltaT,fine_flx,cc_solver_alpha,cc_solver_beta,dx, &
                           the_bc_tower,temp_bc_comp)

          ! debugging statements
          if (.false.) then

             print*,'deltaT norm',multifab_norm_inf_c(deltaT(1),1,1)
             if (l .eq. 1) then
                call fabio_ml_multifab_write_d(deltaT,mla%mba%rr(:,1),"a_deltaT1")
             else if (l .eq. 2) then
                call fabio_ml_multifab_write_d(deltaT,mla%mba%rr(:,1),"a_deltaT2")
             else if (l .eq. 3) then
                call fabio_ml_multifab_write_d(deltaT,mla%mba%rr(:,1),"a_deltaT3")
             else if (l .eq. 4) then
                call fabio_ml_multifab_write_d(deltaT,mla%mba%rr(:,1),"a_deltaT4")
             else if (l .eq. 5) then
                call fabio_ml_multifab_write_d(deltaT,mla%mba%rr(:,1),"a_deltaT5")
             else if (l .eq. 6) then
                call fabio_ml_multifab_write_d(deltaT,mla%mba%rr(:,1),"a_deltaT6")
             else if (l .eq. 7) then
                call fabio_ml_multifab_write_d(deltaT,mla%mba%rr(:,1),"a_deltaT7")
             else if (l .eq. 8) then
                call fabio_ml_multifab_write_d(deltaT,mla%mba%rr(:,1),"a_deltaT8")
             else if (l .eq. 9) then
                call fabio_ml_multifab_write_d(deltaT,mla%mba%rr(:,1),"a_deltaT9")
             end if

          end if

          do n = 2,nlevs
             call bndry_reg_destroy(fine_flx(n))
          end do

          ! T^{n+1,l+1} = T^{n+1,l} + deltaT
          do n=1,nlevs
             call multifab_plus_plus_c(Temp_new(n),1,deltaT(n),1,1,0)
          end do

          ! fill T ghost cells
          do n=1,nlevs
             call multifab_fill_boundary(Temp_new(n))
             call multifab_physbc(Temp_new(n),1,temp_bc_comp,1,the_bc_tower%bc_tower_array(n), &
                                  dx_in=dx(n,:))
          end do

          ! h^{n+1,l+1} = h(rho^{n+1},w^{n+1},T^{n+1,l+1})
          call compute_h(mla,Temp_new,enth_new,conc_new)
          call convert_rhoh_to_h(mla,rhoh_new,rhotot_new,enth_new,.false.)

          ! compute t^{n+1,l} transport properties
          call ideal_mixture_transport_wrapper(mla,rhotot_new,Temp_new,p0_new,conc_new, &
                                               molefrac_new,eta_new,lambda_new,kappa_new, &
                                               chi_new,zeta_new)

          ! rhoh_fluxdiv_new = div(Q)^{n+1,l} + sum(div(hk*Fk))^{n+1,l} + rho_new*Hext^{n+1,l}
          call rhoh_fluxdiv_energy(mla,lambda_new,Temp_new,mass_flux_new,rhotot_new, &
                                   rhoh_fluxdiv_new,dx,0.d0,the_bc_tower)

          ! rhoh_update_new = [-div(rhoh*v) + p0_update + div(Q) + sum(div(hk*Fk)) + rho*Hext]^{n+1}
          do n=1,nlevs
             call multifab_copy_c(rhoh_update_new(n),1,rhoh_fluxdiv_new(n),1,1,0)
             call multifab_plus_plus_s_c(rhoh_update_new(n),1,p0_update_new,1,0)
          end do
          ! FIXME, shouldn't have to recompute this
          call mk_advective_s_fluxdiv(mla,umac_new,rhoh_fc_new,rhoh_update_new,dx,1,1)



       end do  ! end loop l over deltaT_iters


       stop

    end do  ! end loop k over dpdt_iters

    do n=1,nlevs
       call multifab_destroy(conc_old(n))
       call multifab_destroy(conc_new(n))
       call multifab_destroy(molefrac_old(n))
       call multifab_destroy(molefrac_new(n))
       call multifab_destroy(eta_old(n))
       call multifab_destroy(eta_new(n))
       call multifab_destroy(kappa_old(n))
       call multifab_destroy(kappa_new(n))
       call multifab_destroy(lambda_old(n))
       call multifab_destroy(lambda_new(n))
       call multifab_destroy(chi_old(n))
       call multifab_destroy(chi_new(n))
       call multifab_destroy(zeta_old(n))
       call multifab_destroy(zeta_new(n))
       call multifab_destroy(enth_old(n))
       call multifab_destroy(enth_new(n))
       do i=1,dm
          call multifab_destroy(mass_flux_old(n,i))
          call multifab_destroy(mass_flux_new(n,i))
       end do
       call multifab_destroy(mass_fluxdiv_old(n))
       call multifab_destroy(mass_fluxdiv_new(n))
       call multifab_destroy(mass_update_old(n))
       call multifab_destroy(mass_update_new(n))
       call multifab_destroy(rhoh_fluxdiv_old(n))
       call multifab_destroy(rhoh_fluxdiv_new(n))
       call multifab_destroy(rhoh_update_old(n))
       call multifab_destroy(rhoh_update_new(n))
       call multifab_destroy(delta_S_old(n))
       call multifab_destroy(delta_alpha_old(n))
       do i=1,dm
          call multifab_destroy(rhotot_fc_old(n,i))
          call multifab_destroy(rhotot_fc_new(n,i))
          call multifab_destroy(rho_fc_old(n,i))
          call multifab_destroy(rho_fc_new(n,i))
          call multifab_destroy(rhoh_fc_old(n,i))
          call multifab_destroy(rhoh_fc_new(n,i))
       end do
       call multifab_destroy(Scorr(n))
       call multifab_destroy(delta_Scorr(n))
       call multifab_destroy(Peos(n))
       call multifab_destroy(deltaT(n))
       call multifab_destroy(deltaT_rhs(n))
       call multifab_destroy(cc_solver_alpha(n))
       do i=1,dm
          call multifab_destroy(cc_solver_beta(n,i))
       end do
    end do

  end subroutine advance_timestep

end module advance_timestep_module
