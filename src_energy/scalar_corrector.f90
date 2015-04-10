module scalar_corrector_module

  use ml_layout_module
  use bndry_reg_module
  use ml_solve_module
  use multifab_physbc_module
  use bc_module
  use define_bc_module
  use energy_eos_module
  use energy_eos_wrapper_module
  use convert_rhoc_to_c_module
  use convert_rhoh_to_h_module
  use mass_fluxdiv_energy_module
  use rhoh_fluxdiv_energy_module
  use macproject_module
  use convert_stag_module
  use convert_m_to_umac_module
  use div_and_grad_module
  use mk_advective_s_fluxdiv_module
  use mk_advective_m_fluxdiv_module
  use mk_grav_force_module
  use mass_flux_utilities_module
  use multifab_physbc_stag_module
  use gmres_module
  use diffusive_m_fluxdiv_module
  use probin_multispecies_module, only: nspecies
  use probin_common_module, only: n_cells, grav, barodiffusion_type

  use fabio_module

  implicit none

  private

  public :: scalar_corrector

contains

  ! this routine performs "Step 2" of the algorithm (see exec/energy/doc/)
  ! -Compute P_0^{n+1}
  ! -Compute v^{*,n+1} with a Stokes solver
  ! -Compute rho_i^{n+1} explicitly
  ! -Compute (rho h)^{n+1} implicitly
  ! -If necessary, compute volume discrepancy correction and return to beginning of step
  subroutine scalar_corrector(mla,umac_old,umac_new,rho_old,rho_new,rhotot_old,rhotot_new, &
                              rhoh_old,rhoh_new,p0_old,p0_new,pi, &
                              gradp_baro,Temp_old,Temp_new,eta_old,eta_old_ed, &
                              eta_new,eta_new_ed, &
                              mass_update_old,rhoh_update_old,pres_update_old, &
                              dx,dt,time,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: umac_old(:,:)
    type(multifab) , intent(inout) :: umac_new(:,:)
    type(multifab) , intent(inout) :: rho_old(:)
    type(multifab) , intent(inout) :: rho_new(:)
    type(multifab) , intent(inout) :: rhotot_old(:)
    type(multifab) , intent(inout) :: rhotot_new(:)
    type(multifab) , intent(inout) :: rhoh_old(:)
    type(multifab) , intent(inout) :: rhoh_new(:)
    real(kind=dp_t), intent(in   ) :: p0_old
    real(kind=dp_t), intent(inout) :: p0_new
    type(multifab) , intent(inout) :: pi(:)
    type(multifab) , intent(inout) :: gradp_baro(:,:)
    type(multifab) , intent(inout) :: Temp_old(:)
    type(multifab) , intent(inout) :: Temp_new(:)
    ! enters with eta^n
    type(multifab) , intent(inout) :: eta_old(:)
    type(multifab) , intent(inout) :: eta_old_ed(:,:) ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(inout) :: eta_new(:)
    type(multifab) , intent(inout) :: eta_new_ed(:,:) ! nodal (2d); edge-centered (3d)
    ! enters with div(F^n) - div(rho*v)^n
    type(multifab) , intent(inout) :: mass_update_old(:)
    ! enters with [-div(rhoh*v) + (Sbar+Scorrbar)/alphabar + div(Q) + div(h*F) + (rhoHext)]^n
    type(multifab) , intent(inout) :: rhoh_update_old(:)
    ! enters with (Sbar^n + Scorrbar^n) / alphabar^n
    real(kind=dp_t), intent(in   ) :: pres_update_old
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt,time
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables

    ! holds (Sbar^{*,n+1} + Scorrbar^{*,n+1}) / alphabar^{*,n+1}
    real(kind=dp_t) :: pres_update_new

    ! holds the average P0 update at the old and new times
    real(kind=dp_t) :: pres_update_avg

    ! temporary storage for F_k
    type(multifab) :: mass_flux_new(mla%nlevel,mla%dim)

    ! this will hold div(F_k)
    type(multifab) :: mass_fluxdiv_new(mla%nlevel)

    ! this will hold div(Q) + sum(div(hk*Fk)) + rho*Hext
    type(multifab) :: rhoh_fluxdiv_new(mla%nlevel)

    ! holds -div(rho*v)^{*,n+1}
    type(multifab) :: mass_update_new(mla%nlevel)

    ! holds -(1/2)[div(rhoh*v) + (Sbar+Scorrbar)/alphabar]^{*,n+1}
    type(multifab) :: rhoh_update_new(mla%nlevel)

    ! temporary storage for enthalpy
    type(multifab) :: enth(mla%nlevel)

    ! this will hold rho and rhoh on faces (for computing advective fluxes)
    type(multifab) :: rho_fc_new(mla%nlevel,mla%dim)
    type(multifab) :: rhoh_fc_new(mla%nlevel,mla%dim)

    ! temporary storage for enthalpy and temprature
    type(multifab) :: enth_tmp(mla%nlevel)
    type(multifab) :: Temp_tmp(mla%nlevel)

    ! temporary storage for momentum
    type(multifab) :: mtemp(mla%nlevel,mla%dim)

    ! This will hold (rhoh)^n/dt - (1/2)div(rhoh*v)^n - (1/2)div(rhoh*v)^{*,n+1}
    !                + (1/2)(Sbar^n+Sbarcorr^n)/alphabar^n 
    !                + (1/2)(Sbar^{*,n+1}+Sbarcorr^{*,n+1})/alphabar^{*,n+1}
    !                + (1/2)(div(Q^n) + sum(div(h_k^n F_k^n)) + (rho Hext)^n)
    ! for the RHS of the temperature diffusion solve.
    ! Each of these terms stays fixed over all l iterations.
    type(multifab) :: deltaT_rhs1(mla%nlevel)

    ! This will hold -(rho^{*,n+1}h^{*,n+1,l})/dt
    !                + (1/2)(div(Q^{*,n+1,l}) + sum(div(h_k^{*,n+1,l}F_k^{*,n+1,l}))
    !                + (1/2)(rho Hext)^(*,n+1))
    ! for the RHS of the temperature diffusion solve.
    ! Each of these terms may change for each l iteration.
    type(multifab) :: deltaT_rhs2(mla%nlevel)

    ! temporary storage for concentrations and mole fractions
    type(multifab) :: conc    (mla%nlevel)
    type(multifab) :: molefrac(mla%nlevel)

    ! the implicit energy solve computes deltaT
    type(multifab) :: deltaT(mla%nlevel)

    ! coefficients for deltaT (energy) solve
    type(multifab) :: cc_solver_alpha(mla%nlevel)
    type(multifab) :: cc_solver_beta (mla%nlevel,mla%dim)

    ! temporary storage for viscosity
    type(multifab) :: eta(mla%nlevel)

    ! bulk viscosity
    type(multifab) :: kappa(mla%nlevel)

    ! thermal diffusivity
    type(multifab) :: lambda(mla%nlevel)

    ! diffusion matrix
    type(multifab) :: chi(mla%nlevel)

    ! thermodiffusion coefficients
    type(multifab) :: zeta(mla%nlevel)

    ! div(u) + alpha dP_0/dt = S_new, where
    ! S_new = Sbar_new + deltaS_new
    real(kind=dp_t) :: Sbar_new
    type(multifab)  :: deltaS_new(mla%nlevel)

    ! volume discrepancy correction
    ! Scorr_new = Scorrbar_new + deltaScorr_new
    type(multifab)  :: Scorr_new(mla%nlevel)
    real(kind=dp_t) :: Scorrbar_new
    type(multifab)  :: deltaScorr_new(mla%nlevel)

    ! coefficient multiplying dP_0/dt in constraint at new-time
    ! alpha_new = alphabar_new + deltaalpha_new
    type(multifab)  :: alpha_new(mla%nlevel)
    real(kind=dp_t) :: alphabar_new
    type(multifab)  :: deltaalpha_new(mla%nlevel)

    ! this holds the thermodynamic pressure
    type(multifab) :: Peos(mla%nlevel)

    ! coefficient for projection, also used to average rho*h to faces
    type(multifab) :: rhotot_fc_old(mla%nlevel,mla%dim)
    type(multifab) :: rhotot_fc_new(mla%nlevel,mla%dim)

    ! for energy implicit solve
    ! doesn't actually do anything for single-level solves
    type(bndry_reg) :: fine_flx(2:mla%nlevel)

    ! RHS for the Stokes solver
    type(multifab) :: gmres_rhs_p(mla%nlevel)
    type(multifab) :: gmres_rhs_v(mla%nlevel,mla%dim)

    ! temporary copy of incomping pi
    type(multifab) :: pi_tmp(mla%nlevel)
    
    ! stores grad(pi)
    type(multifab) :: gradpi(mla%nlevel,mla%dim)

    ! temporary storage for A_0(u)
    type(multifab) :: m_d_fluxdiv(mla%nlevel,mla%dim)

    ! for stokes solver - will/can be passed in
    type(multifab) :: dumac(mla%nlevel,mla%dim)
    type(multifab) :: dpi(mla%nlevel)

    integer :: n,nlevs,i,dm,n_cell,k,l
    real(kind=dp_t) :: theta_alpha, norm_pre_rhs

    nlevs = mla%nlevel
    dm = mla%dim

    theta_alpha = 1.d0/dt

    if (dm .eq. 2) then
       n_cell = n_cells(1)*n_cells(2)
    else
       n_cell = n_cells(1)*n_cells(2)*n_cells(3)
    end if

    do n=1,nlevs
       do i=1,dm
          call multifab_build_edge(mass_flux_new(n,i),mla%la(n),nspecies,0,i)
       end do
       call multifab_build(mass_fluxdiv_new(n),mla%la(n),nspecies,0)
       call multifab_build(rhoh_fluxdiv_new(n),mla%la(n),1,0)

       call multifab_build(mass_update_new(n),mla%la(n),nspecies,0)
       call multifab_build(rhoh_update_new(n),mla%la(n),1       ,0)

       call multifab_build(enth(n),mla%la(n),1,rho_old(n)%ng)

       do i=1,dm
          call multifab_build_edge(rho_fc_new(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge(rhoh_fc_new(n,i),mla%la(n),1,0,i)

          call multifab_build_edge(mtemp(n,i),mla%la(n),1,1,i)
       end do

       call multifab_build(enth_tmp(n),mla%la(n),1,rhoh_new(n)%ng)
       call multifab_build(Temp_tmp(n),mla%la(n),1,Temp_new(n)%ng)

       call multifab_build(deltaT_rhs1(n),mla%la(n),1,0)
       call multifab_build(deltaT_rhs2(n),mla%la(n),1,0)

       call multifab_build(    conc(n),mla%la(n),nspecies,rho_old(n)%ng)
       call multifab_build(molefrac(n),mla%la(n),nspecies,rho_old(n)%ng)

       call multifab_build(deltaT(n),mla%la(n),1,1)

       call multifab_build(cc_solver_alpha(n),mla%la(n),1,0)
       do i=1,dm
          call multifab_build_edge(cc_solver_beta(n,i),mla%la(n),1,0,i)
       end do

       call multifab_build(eta(n),mla%la(n),1,1)

       call multifab_build(kappa(n),mla%la(n),1,1)

       call multifab_build(lambda(n),mla%la(n),1,2)

       call multifab_build(chi(n),mla%la(n),nspecies**2,1)

       call multifab_build(zeta(n),mla%la(n),nspecies,1)

       call multifab_build(deltaS_new(n),mla%la(n),1,0)

       call multifab_build(     Scorr_new(n),mla%la(n),1,0)
       call multifab_build(deltaScorr_new(n),mla%la(n),1,0)

       call multifab_build(     alpha_new(n),mla%la(n),1,0)
       call multifab_build(deltaalpha_new(n),mla%la(n),1,0)

       call multifab_build(Peos(n),mla%la(n),1,0)

       do i=1,dm
          call multifab_build_edge(rhotot_fc_old(n,i),mla%la(n),1,1,i)
          call multifab_build_edge(rhotot_fc_new(n,i),mla%la(n),1,1,i)

          call multifab_build_edge(gmres_rhs_v(n,i),mla%la(n),1,0,i)
       end do

       call multifab_build(pi_tmp(n),mla%la(n),1,1)

       call multifab_build(gmres_rhs_p(n),mla%la(n),1,0)
       do i=1,dm
          call multifab_build_edge(gradpi(n,i),mla%la(n),1,0,i)
          call multifab_build_edge(m_d_fluxdiv(n,i),mla%la(n),1,0,i)
       end do

       call multifab_build(dpi(n),mla%la(n),1,1)
       do i=1,dm
          call multifab_build_edge(dumac(n,i),mla%la(n),1,1,i)
       end do

    end do

    ! temporary copies
    do n=1,nlevs
       call multifab_copy_c(pi_tmp(n),1,pi(n),1,1,1)
       call multifab_copy_c(Temp_tmp(n),1,Temp_new(n),1,1,Temp_tmp(n)%ng)
    end do

    ! compute mass fractions in valid region and then fill ghost cells
    call convert_rhoc_to_c(mla,rho_new,rhotot_new,conc,.true.)
    call fill_c_ghost_cells(mla,conc,dx,the_bc_tower)

    ! compute mole fractions in VALID + GHOST regions
    call convert_conc_to_molefrac(mla,conc,molefrac,.true.)

    ! compute t^{*,n+1} transport properties
    call ideal_mixture_transport_wrapper(mla,rhotot_new,Temp_new,p0_new,conc,molefrac, &
                                         eta_new,lambda,kappa,chi,zeta)

    ! eta^{*,n+1} on nodes (2d) or edges (3d)
    if (dm .eq. 2) then
       call average_cc_to_node(nlevs,eta_new,eta_new_ed(:,1),1,tran_bc_comp,1, &
                               the_bc_tower%bc_tower_array)
    else if (dm .eq. 3) then
       call average_cc_to_edge(nlevs,eta_new,eta_new_ed,1,tran_bc_comp,1, &
                               the_bc_tower%bc_tower_array)
    end if

    ! compute mass_fluxdiv_new = div(F^{*,n+1})
    call mass_fluxdiv_energy(mla,rho_new,rhotot_new,molefrac,chi,zeta, &
                             gradp_baro,Temp_new,mass_fluxdiv_new,mass_flux_new,dx,the_bc_tower)

    ! compute rhoh_fluxdiv_new = div(Q)^{*,n+1} + sum(div(hk*Fk))^{*,n+1} + rho*Hext^{*,n+1}
    call rhoh_fluxdiv_energy(mla,lambda,Temp_new,mass_flux_new,rhotot_new,rhoh_fluxdiv_new, &
                             dx,time,the_bc_tower)

    ! compute S_new and alpha_new (store them in deltaS_new and deltaalpha_new)
    call compute_S_alpha(mla,deltaS_new,deltaalpha_new,mass_fluxdiv_new,rhoh_fluxdiv_new,conc, &
                         Temp_new,rhotot_new,p0_new)

    ! split S_new and alpha_new into average and perturbational pieces
    ! S_new = Sbar_new + deltaS_new
    ! alpha_new = alphabar_new + deltaalpha_new
    do n=1,nlevs
       Sbar_new     = multifab_sum_c(deltaS_new(n)    ,1,1) / dble(n_cell)
       alphabar_new = multifab_sum_c(deltaalpha_new(n),1,1) / dble(n_cell)
       call multifab_sub_sub_s_c(deltaS_new(n)    ,1,Sbar_new    ,1,0)
       call multifab_sub_sub_s_c(deltaalpha_new(n),1,alphabar_new,1,0)
    end do

    ! zero out volume discrepancy correction and its decomposition
    Scorrbar_new = 0.d0
    do n=1,nlevs
       call multifab_setval(Scorr_new(n)     ,0.d0,all=.true.)
       call multifab_setval(deltaScorr_new(n),0.d0,all=.true.)
    end do

    ! set rho_fc_new to rho^{*,n+1} on faces
    call average_cc_to_face(nlevs,rho_new,rho_fc_new,1,c_bc_comp,nspecies, &
                            the_bc_tower%bc_tower_array)

    ! set rhotot_fc_old to rhotot^n on faces
    call average_cc_to_face(nlevs,rhotot_old,rhotot_fc_old,1,scal_bc_comp,1, &
                            the_bc_tower%bc_tower_array)

    ! set rhotot_fc_new to rhotot^{*,n+1} on faces
    call average_cc_to_face(nlevs,rhotot_new,rhotot_fc_new,1,scal_bc_comp,1, &
                            the_bc_tower%bc_tower_array)

    ! compute h^{*,n+1} and fill ghost cells
    call convert_rhoh_to_h(mla,rhoh_new,rhotot_new,enth,.true.)
    call fill_h_ghost_cells(mla,enth,dx,the_bc_tower)

    ! average h^{*,n+1} to faces
    call average_cc_to_face(nlevs,enth,rhoh_fc_new,1,h_bc_comp,1, &
                            the_bc_tower%bc_tower_array)

    ! multiply h^{*,n+1} on faces by rhotot^{*,n+1} on faces
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_c(rhoh_fc_new(n,i),1,rhotot_fc_new(n,i),1,1,0)
       end do
    end do

    ! begin loop here over Steps 1a-1e
    do k=1,dpdt_iters

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Step 2a: Compute a pressure update
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! update pressure
       pres_update_new = (Sbar_new + Scorrbar_new)/alphabar_new
       pres_update_avg = 0.5d0*(pres_update_old+pres_update_new)
       p0_new = p0_old + 0.5d0*dt*pres_update_avg

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Step 2b: Compute the velocity and dynamic pressure using a Stokes solver
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! compute gmres_rhs_p = deltaS_new + deltaScorr_new 
       !                       - deltaalpha_new * (Sbar_new + Scorrbar_new)/alphabar_new
       do n=1,nlevs
          call multifab_copy_c(gmres_rhs_p(n),1,deltaalpha_new(n),1,1,0)
          call multifab_mult_mult_s_c(gmres_rhs_p(n),1,-pres_update_new,1,0)
          call multifab_plus_plus_c(gmres_rhs_p(n),1,deltaS_new(n),1,1,0)
          call multifab_plus_plus_c(gmres_rhs_p(n),1,deltaScorr_new(n),1,1,0)
       end do

       ! construct gmres_rhs_p = div(vbar^n) - gmres_rhs_p
       ! first multiply gmres_rhs_p by -1
       do n=1,nlevs
          call multifab_mult_mult_s_c(gmres_rhs_p(n),1,-1.d0,1,0)
       end do

       ! FIXME: umac_old will need vbar boundary conditions
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(umac_new(n,i),1,umac_old(n,i),1,1,umac_new(n,i)%ng)
          end do
       end do

       ! add div(vbar) to gmres_rhs_p
       call compute_div(mla,umac_new,gmres_rhs_p,dx,1,1,1,increment_in=.true.)

       ! construct gmres_rhs_v
       ! set gmres_rhs_v to rho^n v^n
       call convert_m_to_umac(mla,rhotot_fc_old,mtemp,umac_old,.false.)
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(gmres_rhs_v(n,i),1,mtemp(n,i),1,1,0)
          end do
       end do

       ! subtract mtemp = rho^{*,n+1} vbar^n from gmres_rhs_v and then multiply by 1/dt
       call convert_m_to_umac(mla,rhotot_fc_new,mtemp,umac_new,.false.)
       do n=1,nlevs
          do i=1,dm
             call multifab_sub_sub_c(gmres_rhs_v(n,i),1,mtemp(n,i),1,1,0)
             call multifab_mult_mult_s_c(gmres_rhs_v(n,i),1,1.d0/dt,1,0)
          end do
       end do

       ! reset pi back to pi^n
       do n=1,nlevs
          call multifab_copy_c(pi(n),1,pi_tmp(n),1,1,1)
       end do

       ! subtract grad(pi^n) from gmres_rhs_v
       call compute_grad(mla,pi,gradpi,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)
       do n=1,nlevs
          do i=1,dm
             call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradpi(n,i),1,1,0)
          end do
       end do

       if (barodiffusion_type .ne. 0) then
          call bl_error("scalar_corrector: barodiffusion not implemented yet")
       end if

       ! add -div(rho*v*v)^n to gmres_rhs_v
       call convert_m_to_umac(mla,rhotot_fc_old,mtemp,umac_old,.false.)
       call mk_advective_m_fluxdiv(mla,umac_old,mtemp,gmres_rhs_v,dx, &
                                the_bc_tower%bc_tower_array)

       ! add (1/2) A_0^n v^n to gmres_rhs_v
       do n=1,nlevs
          do i=1,dm
             call setval(m_d_fluxdiv(n,i),0.d0,all=.true.)
          end do
       end do
       call diffusive_m_fluxdiv(mla,m_d_fluxdiv,umac_old,eta_old,eta_old_ed,kappa,dx, &
                                the_bc_tower%bc_tower_array)
       do n=1,nlevs
          do i=1,dm
             call multifab_saxpy_3(gmres_rhs_v(n,i),0.5d0,m_d_fluxdiv(n,i))
          end do
       end do

       ! add (1/2) A_0^{*,n+1} vbar^n to gmres_rhs_v
       do n=1,nlevs
          do i=1,dm
             call setval(m_d_fluxdiv(n,i),0.d0,all=.true.)
          end do
       end do
       call diffusive_m_fluxdiv(mla,m_d_fluxdiv,umac_new,eta_new,eta_new_ed,kappa,dx, &
                                the_bc_tower%bc_tower_array)
       do n=1,nlevs
          do i=1,dm
             call multifab_saxpy_3(gmres_rhs_v(n,i),0.5d0,m_d_fluxdiv(n,i))
          end do
       end do

       ! add (1/2)(rho^n + rho^{*,n+1})g to gmres_rhs_v
       if (any(grav(1:dm) .ne. 0.d0)) then
          call mk_grav_force(mla,gmres_rhs_v,rhotot_fc_old,rhotot_fc_new,the_bc_tower)
       end if

       ! set the initial guess to zero
       do n=1,nlevs
          do i=1,dm
             call multifab_setval(dumac(n,i),0.d0,all=.true.)
          end do
          call multifab_setval(dpi(n),0.d0,all=.true.)
       end do

       ! multiply eta^{*,n+1} and kappa by 1/2 to put in proper form for gmres solve
       do n=1,nlevs
          call multifab_mult_mult_s_c(eta_new(n),1,1.d0/2.d0,1,eta_new(n)%ng)
          call multifab_mult_mult_s_c(kappa(n),1,1.d0/2.d0,1,kappa(n)%ng)
          do i=1,size(eta_new_ed,dim=2)
             call multifab_mult_mult_s_c(eta_new_ed(n,i),1,1.d0/2.d0,1,eta_new_ed(n,i)%ng)
          end do
       end do

       ! call the Stokes solver
       call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dpi,rhotot_fc_new, &
                  eta_new,eta_new_ed,kappa,theta_alpha,norm_pre_rhs)

       ! restore eta and kappa
       do n=1,nlevs
          call multifab_mult_mult_s_c(eta_new(n),1,2.d0,1,eta_new(n)%ng)
          call multifab_mult_mult_s_c(kappa(n),1,2.d0,1,kappa(n)%ng)
          do i=1,size(eta_new_ed,dim=2)
             call multifab_mult_mult_s_c(eta_new_ed(n,i),1,2.d0,1,eta_new_ed(n,i)%ng)
          end do
       end do

       ! increment velocity and dynamic pressure
       ! compute v^{*,n+1} = v^n + dumac
       ! compute pi^{*,n+1}= pi^n + dpi
       do n=1,nlevs
          do i=1,dm
             call multifab_plus_plus_c(umac_new(n,i),1,dumac(n,i),1,1,0)
          end do
          call multifab_plus_plus_c(pi(n),1,dpi(n),1,1,0)
       end do

       do n=1,nlevs
          ! presure ghost cells
          call multifab_fill_boundary(pi(n))
          call multifab_physbc(pi(n),1,pres_bc_comp,1,the_bc_tower%bc_tower_array(n), &
                               dx_in=dx(n,:))
          do i=1,dm
             ! set normal velocity on physical domain boundaries
             call multifab_physbc_domainvel(umac_new(n,i),vel_bc_comp+i-1, &
                                            the_bc_tower%bc_tower_array(n), &
                                            dx(n,:))
             ! set transverse velocity behind physical boundaries
             call multifab_physbc_macvel(umac_new(n,i),vel_bc_comp+i-1, &
                                         the_bc_tower%bc_tower_array(n), &
                                         dx(n,:))
             ! fill periodic and interior ghost cells
             call multifab_fill_boundary(umac_new(n,i))
          end do
       end do

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Step 2c: Advance the densities using trapezoidal advective and diffusive fluxes
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! zero out mass update so we can increment it
       do n=1,nlevs
          call multifab_setval(mass_update_new(n),0.d0,all=.true.)
       end do

       ! compute -div(rho*v)^{*,n+1} and store it in mass_update_new
       call mk_advective_s_fluxdiv(mla,umac_new,rho_fc_new,mass_update_new,dx,1,nspecies)

       ! add div(F)^{*,n+1} to mass_update_new
       do n=1,nlevs
          call multifab_plus_plus_c(mass_update_new(n),1,mass_fluxdiv_new(n),1,nspecies,0)
       end do

       ! rho_new = rho_old + (dt/2)[-div(rho*v)^n - div(rho*v)^{*,n+1} + div(F^n) + div(F^{*,n+1})]
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
       call convert_rhoc_to_c(mla,rho_new,rhotot_new,conc,.true.)
       call fill_c_ghost_cells(mla,conc,dx,the_bc_tower)
       call convert_rhoc_to_c(mla,rho_new,rhotot_new,conc,.false.)

       ! compute mole fractions in VALID + GHOST regions
       call convert_conc_to_molefrac(mla,conc,molefrac,.true.)

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Step 2d: Advance the enthalpy by iteratively looping over an energy solve
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! The portion of the RHS that stays fixed over all l iterations is:
       !   (rhoh)^n/dt - (1/2)div(rhoh*v)^n - (1/2)div(rhoh*v)^{*,n+1}
       !               + (1/2)(Sbar^n+Sbarcorr^n)/alphabar^n 
       !               + (1/2)(Sbar^{*,n+1}+Sbarcorr^{*,n+1})/alphabar^{*,n+1}
       !               + (1/2)(div(Q^n) + sum(div(h_k^n F_k^n)) + (rho Hext)^n)
       ! store this in deltaT_rhs1

       ! rhoh_update_old already contains
       !   [-div(rhoh*v) + (Sbar+Scorrbar)/alphabar + div(Q) + div(h*F) + (rhoHext)]^n

       ! rhoh_update_new is a temporary that will hold 
       !   -(1/2)[div(rhoh*v) + (Sbar+Scorrbar)/alphabar]^{*,n+1}

       ! set rhoh_update_new to [Sbar+Scorrbar)/alphabar]^{*,n+1}
       do n=1,nlevs
          call multifab_setval(rhoh_update_new(n),pres_update_new,all=.true.)
       end do

       ! add -div(rhoh*v)^{*,n+1} to rhoh_update_new
       call mk_advective_s_fluxdiv(mla,umac_new,rhoh_fc_new,rhoh_update_new,dx,1,1)

       ! deltaT_rhs1 = (rhoh)^n/dt
       do n=1,nlevs
          call multifab_copy_c(deltaT_rhs1(n),1,rhoh_old(n),1,1,0)
          call multifab_mult_mult_s_c(deltaT_rhs1(n),1,1.d0/dt,1,0)
       end do

       ! add (1/2)[-div(rhoh*v) + [Sbar+Scorrbar)/alphabar]^{*,n+1} to deltaT_rhs1
       do n=1,nlevs
          call multifab_saxpy_3(deltaT_rhs1(n),0.5d0,rhoh_update_new(n))
       end do

       ! add (1/2)[-div(rhoh*v) + (Sbar+Scorrbar)/alphabar + div(Q) + div(h*F) + (rhoHext)]^n
       ! to deltaT_rhs1
       do n=1,nlevs
          call multifab_saxpy_3(deltaT_rhs1(n),0.5d0,rhoh_update_old(n))
       end do

       ! reset Temp_new back to T^{*,n+1}
       do n=1,nlevs
          call multifab_copy_c(Temp_new(n),1,Temp_tmp(n),1,1,Temp_tmp(n)%ng)
       end do

       ! set rhoh_new = rho^{n+1} h^{n+1,l} = rho^{n+1} h^{*,n+1}
       do n=1,nlevs
          call multifab_copy_c(rhoh_new(n),1,enth(n),1,1,rhoh_new(n)%ng)
          call multifab_mult_mult_c(rhoh_new(n),1,rhotot_new(n),1,1,rhoh_new(n)%ng)
       end do

       do l=1,deltaT_iters

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Step 2d-1: Compute (lambda,cp,F)^{n+1,l}
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! compute t^{n+1,l} transport properties
          ! don't overwrite eta_new since that is needed at the beginning of the next dpdt_iter
          call ideal_mixture_transport_wrapper(mla,rhotot_new,Temp_new,p0_new,conc,molefrac, &
                                               eta,lambda,kappa,chi,zeta)

          ! compute mass_flux_new = F^{n+1,l}
          ! compute mass_fluxdiv_new = div(F^{n+1,l}) (not actually needed)
          call mass_fluxdiv_energy(mla,rho_new,rhotot_new,molefrac,chi,zeta, &
                                   gradp_baro,Temp_new,mass_fluxdiv_new,mass_flux_new,dx,the_bc_tower)

          ! compute rhoh_fluxdiv_new = div(Q)^{n+1,l} + sum(div(hk*Fk))^{n+1,l} + rho*Hext^{n+1,l}
          call rhoh_fluxdiv_energy(mla,lambda,Temp_new,mass_flux_new,rhotot_new,rhoh_fluxdiv_new, &
                                   dx,time,the_bc_tower)

          ! The portion of the RHS that changes for each l iteration is:
          !  -(rho^{n+1}h^{n+1,l})/dt 
          !               + (1/2)(div(Q^{n+1,l}) + sum(div(h_k^{n+1,l}F_k^{n+1,l})) + (rho Hext)^(n+1))
          ! store this in deltaT_rhs2

          ! set deltaT_rhs2 = -(rho^{n+1} h^{n+1,l})/dt 
          do n=1,nlevs
             call multifab_copy_c(deltaT_rhs2(n),1,rhoh_new(n),1,1,0)
             call multifab_mult_mult_s_c(deltaT_rhs2(n),1,-1.d0/dt,1,0)
          end do

          ! add (1/2)[div(Q) + sum(div(h_k F_k)) + (rho Hext)]^{n+1,l} to deltaT_rhs2
          do n=1,nlevs
             call multifab_saxpy_3(deltaT_rhs2(n),0.5d0,rhoh_fluxdiv_new(n))
          end do

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Step 2d-2: Solve for deltaT implicitly
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! cc_solver_alpha = rho^{n+1} c_p^{n+1,l} / dt
          call compute_cp(mla,cc_solver_alpha,conc,Temp_new)
          do n=1,nlevs
             call multifab_mult_mult_c(cc_solver_alpha(n),1,rhotot_new(n),1,1,0)
             call multifab_mult_mult_s_c(cc_solver_alpha(n),1,1.d0/dt,1,0)
          end do

          ! cc_solver_beta = (1/2) lambda^{n+1,l}
          call average_cc_to_face(nlevs,lambda,cc_solver_beta,1,tran_bc_comp,1, &
                                  the_bc_tower%bc_tower_array)
          
          do n=1,nlevs
             do i=1,dm
                call multifab_mult_mult_s_c(cc_solver_beta(n,i),1,0.5d0,1,0)
             end do
          end do

          ! cc_solver_rhs = deltaT_rhs1 + deltaT_rhs2 (store this in deltaT_rhs2)
          do n=1,nlevs
             call multifab_plus_plus_c(deltaT_rhs2(n),1,deltaT_rhs1(n),1,1,0)
          end do

          ! solve for deltaT
          do n = 2,nlevs
             call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
          end do

          do n=1,nlevs
             call setval(deltaT(n),0.d0,all=.true.)
          end do

          call ml_cc_solve(mla,deltaT_rhs2,deltaT,fine_flx,cc_solver_alpha,cc_solver_beta,dx, &
                           the_bc_tower,temp_bc_comp)

          do n = 2,nlevs
             call bndry_reg_destroy(fine_flx(n))
          end do

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Step 2d-3: Update the temperature and enthalpy
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
          call compute_h(mla,Temp_new,enth_tmp,conc)
          call convert_rhoh_to_h(mla,rhoh_new,rhotot_new,enth_tmp,.false.)

       end do ! end loop l over deltaT iterations

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Step 2e: If the thermodynamic drift is unacceptable, update the volume
       !          discrepancy correction
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! compute thermodynamic pressure
       call compute_p(mla,rhotot_new,Temp_new,conc,Peos)

       ! compute alpha^{*,n+1}
       call compute_alpha(mla,alpha_new,conc,Temp_new,p0_new)

       ! Scorr = Scorr + alpha^{n+1} * 2*[(Peos^{n+1} - P0^{n+1})/dt]
       do n=1,nlevs
          call multifab_sub_sub_s_c(Peos(n),1,p0_new,1,0)

!          if (k .eq. 1) then
!             call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift1")
!          else if (k .eq. 2) then
!             call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift2")
!          else if (k .eq. 3) then
!             call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift3")
!          else if (k .eq. 4) then
!             call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift4")
!          else if (k .eq. 5) then
!             call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift5")
!          else if (k .eq. 6) then
!             call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift6")
!          else if (k .eq. 7) then
!             call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift7")
!          else if (k .eq. 8) then
!             call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift8")
!          else if (k .eq. 9) then
!             call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift9")
!          end if

          call multifab_mult_mult_s_c(Peos(n),1,2.d0/dt,1,0)
          call multifab_mult_mult_c(Peos(n),1,alpha_new(n),1,1,0)
          call multifab_plus_plus_c(Scorr_new(n),1,Peos(n),1,1,0)
       end do

       ! split Scorr into average and perturbational components
       do n=1,nlevs
          Scorrbar_new = multifab_sum_c(Scorr_new(n),1,1) / dble(n_cell)
          call multifab_copy_c(deltaScorr_new(n),1,Scorr_new(n),1,1,0)
          call multifab_sub_sub_s_c(deltaScorr_new(n),1,Scorrbar_new,1,0)
       end do

    end do  ! end loop k over dpdt iterations

    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(mass_flux_new(n,i))
       end do
       call multifab_destroy(mass_fluxdiv_new(n))
       call multifab_destroy(rhoh_fluxdiv_new(n))
       call multifab_destroy(mass_update_new(n))
       call multifab_destroy(rhoh_update_new(n))
       call multifab_destroy(enth(n))
       do i=1,dm
          call multifab_destroy(rho_fc_new(n,i))
          call multifab_destroy(rhoh_fc_new(n,i))
          call multifab_destroy(mtemp(n,i))
       end do
       call multifab_destroy(enth_tmp(n))
       call multifab_destroy(deltaT_rhs1(n))
       call multifab_destroy(deltaT_rhs2(n))
       call multifab_destroy(conc(n))
       call multifab_destroy(molefrac(n))
       call multifab_destroy(deltaT(n))
       call multifab_destroy(cc_solver_alpha(n))
       do i=1,dm
          call multifab_destroy(cc_solver_beta(n,i))
       end do
       call multifab_destroy(eta(n))
       call multifab_destroy(kappa(n))
       call multifab_destroy(lambda(n))
       call multifab_destroy(chi(n))
       call multifab_destroy(zeta(n))
       call multifab_destroy(deltaS_new(n))
       call multifab_destroy(Scorr_new(n))
       call multifab_destroy(deltaScorr_new(n))
       call multifab_destroy(alpha_new(n))
       call multifab_destroy(deltaalpha_new(n))
       call multifab_destroy(Peos(n))
       do i=1,dm
          call multifab_destroy(rhotot_fc_old(n,i))
          call multifab_destroy(rhotot_fc_new(n,i))
          call multifab_destroy(gmres_rhs_v(n,i))
       end do
       call multifab_destroy(pi_tmp(n))
       call multifab_destroy(gmres_rhs_p(n))
       do i=1,dm
          call multifab_destroy(gradpi(n,i))
          call multifab_destroy(m_d_fluxdiv(n,i))
       end do
       call multifab_destroy(dpi(n))
       do i=1,dm
          call multifab_destroy(dumac(n,i))
       end do
    end do

  end subroutine scalar_corrector

end module scalar_corrector_module
