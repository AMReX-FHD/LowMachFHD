module advance_timestep_module

  use ml_layout_module
  use bndry_reg_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use convert_stag_module
  use convert_rhoc_to_c_module
  use convert_rhoh_to_h_module
  use eos_model_wrapper_module
  use mass_flux_utilities_module
  use mass_fluxdiv_energy_module
  use rhoh_fluxdiv_energy_module
  use div_and_grad_module
  use macproject_module
  use mk_advective_s_fluxdiv_module
  use mk_advective_m_fluxdiv_module
  use diffusive_m_fluxdiv_module
  use mk_grav_force_module
  use ml_solve_module
  use gmres_module
  use convert_m_to_umac_module
  use multifab_physbc_stag_module
  use probin_common_module, only: n_cells, grav
  use probin_multispecies_module, only: nspecies
  use probin_energy_module, only: dpdt_iters, dpdt_factor, deltaT_iters

  use fabio_module

  implicit none

  private

  public :: advance_timestep

contains

  subroutine advance_timestep(mla,umac_old,umac_new,rho_old,rho_new, &
                              rhotot_old,rhotot_new,rhoh_old,rhoh_new, &
                              p0_old,p0_new,gradp_baro,Temp_old,Temp_new, &
                              pi,dx,dt,the_bc_tower)

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
    type(multifab) , intent(in   ) :: gradp_baro(:,:)
    type(multifab) , intent(in   ) :: Temp_old(:)
    type(multifab) , intent(inout) :: Temp_new(:)
    type(multifab) , intent(inout) :: pi(:)
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

    ! RHS for the Stokes solver
    type(multifab) :: gmres_rhs_v(mla%nlevel,mla%dim)
    type(multifab) :: gmres_rhs_p(mla%nlevel)
    
    ! stores grad(pi)
    type(multifab) :: gradpi(mla%nlevel,mla%dim)

    ! temporary storage for momentum-like quantities
    type(multifab) :: mtemp(mla%nlevel,mla%dim)
    type(multifab) :: mtemp2(mla%nlevel,mla%dim)

    ! for stokes solver - will/can be passed in
    type(multifab) :: dumac(mla%nlevel,mla%dim)
    type(multifab) :: dpi(mla%nlevel)

    type(multifab), allocatable :: eta_ed_old(:,:)
    type(multifab), allocatable :: eta_ed_new(:,:)

    ! for energy implicit solve
    ! doesn't actually do anything for single-level solves
    type(bndry_reg) :: fine_flx(2:mla%nlevel)

    real(kind=dp_t) :: Sbar_old, Sbar_new
    real(kind=dp_t) :: alphabar_old, alphabar_new
    real(kind=dp_t) :: Scorrbar

    real(kind=dp_t) :: p0_update_old
    real(kind=dp_t) :: p0_update_new

    real(kind=dp_t) :: theta_alpha, norm_pre_rhs, norm

    logical :: nodal_temp(3)

    nlevs = mla%nlevel
    dm = mla%dim

    if (dm .eq. 2) then
       n_cell = n_cells(1)*n_cells(2)
    else
       n_cell = n_cells(1)*n_cells(2)*n_cells(3)
    end if

    if (dm .eq. 2) then
       allocate(eta_ed_old(nlevs,1))
       allocate(eta_ed_new(nlevs,1))
    else if (dm .eq. 3) then
       allocate(eta_ed_old(nlevs,3))
       allocate(eta_ed_new(nlevs,3))
    end if

    theta_alpha = 1.d0/dt

    do n=1,nlevs
       call multifab_build(conc_old(n)    ,mla%la(n),nspecies   ,rho_old(n)%ng)
       call multifab_build(conc_new(n)    ,mla%la(n),nspecies   ,rho_old(n)%ng)
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
       call multifab_build(enth_old(n)    ,mla%la(n),1          ,rho_old(n)%ng)
       call multifab_build(enth_new(n)    ,mla%la(n),1          ,rho_old(n)%ng)
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
          call multifab_build_edge(rhotot_fc_old(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(rhotot_fc_new(n,i),mla%la(n),1       ,1,i)
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
          call multifab_build_edge(gmres_rhs_v(n,i)   ,mla%la(n),1,0,i)
       end do
       call multifab_build(gmres_rhs_p(n),mla%la(n),1,0)
       do i=1,dm
          call multifab_build_edge(gradpi(n,i),mla%la(n),1,0,i)
          call multifab_build_edge(mtemp(n,i) ,mla%la(n),1,1,i)
          call multifab_build_edge(mtemp2(n,i),mla%la(n),1,1,i)
          call multifab_build_edge(dumac(n,i) ,mla%la(n),1,1,i)
       end do
       call multifab_build(dpi(n),mla%la(n),1,1)

       ! eta_old on nodes (2d) or edges (3d)
       if (dm .eq. 2) then
          call multifab_build_nodal(eta_ed_old(n,1),mla%la(n),1,0)
          call multifab_build_nodal(eta_ed_new(n,1),mla%la(n),1,0)
       else
          nodal_temp(1) = .true.
          nodal_temp(2) = .true.
          nodal_temp(3) = .false.
          call multifab_build(eta_ed_old(n,1),mla%la(n),1,0,nodal_temp)
          call multifab_build(eta_ed_new(n,1),mla%la(n),1,0,nodal_temp)
          nodal_temp(1) = .true.
          nodal_temp(2) = .false.
          nodal_temp(3) = .true.
          call multifab_build(eta_ed_old(n,2),mla%la(n),1,0,nodal_temp)
          call multifab_build(eta_ed_new(n,2),mla%la(n),1,0,nodal_temp)
          nodal_temp(1) = .false.
          nodal_temp(2) = .true.
          nodal_temp(3) = .true.
          call multifab_build(eta_ed_old(n,3),mla%la(n),1,0,nodal_temp)
          call multifab_build(eta_ed_new(n,3),mla%la(n),1,0,nodal_temp)
       end if

    end do

    ! temporary boundary register needed for deltaT solve
    do n=2,nlevs
       call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
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

    ! compute mole fractions, x^n, in VALID + GHOST regions
    call convert_conc_to_molefrac(mla,conc_old,molefrac_old,.true.)

    ! compute t^n transport properties (eta,lambda,kappa,chi,zeta)
    call ideal_mixture_transport_wrapper(mla,rhotot_old,Temp_old,p0_old,conc_old, &
                                         molefrac_old,eta_old,lambda_old,kappa_old, &
                                         chi_old,zeta_old)

    ! eta^n on nodes (2d) or edges (3d)
    if (dm .eq. 2) then
       call average_cc_to_node(nlevs,eta_old,eta_ed_old(:,1),1,tran_bc_comp,1, &
                               the_bc_tower%bc_tower_array)
    else if (dm .eq. 3) then
       call average_cc_to_edge(nlevs,eta_old,eta_ed_old,1,tran_bc_comp,1, &
                               the_bc_tower%bc_tower_array)
    end if

    ! compute mass_flux_old = F^n and mass_fluxdiv_old = div(F^n)
    call mass_fluxdiv_energy(mla,rho_old,rhotot_old,molefrac_old,chi_old,zeta_old, &
                             gradp_baro,Temp_old,mass_fluxdiv_old, &
                             mass_flux_old,dx,the_bc_tower)

    ! compute rhoh_fluxdiv_old = div(Q)^n + sum(div(hk*Fk))^n + (rho*Hext)^n
    call rhoh_fluxdiv_energy(mla,lambda_old,Temp_old,mass_flux_old,rhotot_old, &
                             rhoh_fluxdiv_old,dx,0.d0,the_bc_tower)

    ! compute S^n and alpha^n (store them in delta_S_old and delta_alpha_old)
    call compute_S_alpha(mla,delta_S_old,delta_alpha_old,mass_fluxdiv_old, &
                         rhoh_fluxdiv_old,conc_old,Temp_old,rhotot_old,p0_old)

    ! compute P_eos^n
    call compute_p(mla,rhotot_old,Temp_old,conc_old,Peos)

    ! Scorr = Scorr + (1 / p0^n) * (Peos^n - P0^n)/dt
    do n=1,nlevs
       call multifab_sub_sub_s_c(Peos(n),1,p0_old,1,0)
       call multifab_mult_mult_s_c(Peos(n),1,1.d0/(p0_old*dt),1,0)
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
          call multifab_copy_c(umac_new(n,i)     ,1,umac_old(n,i)     ,1,1       ,umac_new(n,i)%ng)
          call multifab_copy_c(rhotot_fc_new(n,i),1,rhotot_fc_old(n,i),1,1       ,rhotot_fc_new(n,i)%ng)
          call multifab_copy_c(rho_fc_new(n,i)   ,1,rho_fc_old(n,i)   ,1,nspecies,rho_fc_new(n,i)%ng)
          call multifab_copy_c(rhoh_fc_new(n,i)  ,1,rhoh_fc_old(n,i)  ,1,1       ,rhoh_fc_new(n,i)%ng)
          call multifab_copy_c(mass_flux_new(n,i),1,mass_flux_old(n,i),1,nspecies,0)
       end do
       call multifab_copy_c(eta_new(n)         ,1,eta_old(n)         ,1,1          ,1)
       call multifab_copy_c(lambda_new(n)      ,1,lambda_old(n)      ,1,1          ,1)
       call multifab_copy_c(kappa_new(n)       ,1,kappa_old(n)       ,1,1          ,1)
       call multifab_copy_c(chi_new(n)         ,1,chi_old(n)         ,1,nspecies**2,1)
       call multifab_copy_c(zeta_new(n)        ,1,zeta_old(n)        ,1,nspecies   ,1)
       call multifab_copy_c(mass_fluxdiv_new(n),1,mass_fluxdiv_old(n),1,nspecies   ,0)
       call multifab_copy_c(rhoh_fluxdiv_new(n),1,rhoh_fluxdiv_old(n),1,1          ,0)
       call multifab_copy_c(delta_S_new(n)     ,1,delta_S_old(n)     ,1,1          ,0)
       call multifab_copy_c(delta_alpha_new(n) ,1,delta_alpha_old(n) ,1,1          ,0)
    end do

    p0_new = p0_old
    Sbar_new = Sbar_old
    alphabar_new = alphabar_old

    ! compute grad(pi^n)
    call compute_grad(mla,pi,gradpi,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

    do k=1,dpdt_iters

       ! mtemp2 will hold -div(rho*v*v)^n - div(rho*v*v)^{n+1,m}
       do n=1,nlevs
          do i=1,dm
             call multifab_setval(mtemp2(n,i),0.d0,all=.true.)
          end do
       end do

       ! add -div(rho*v*v)^{n+1,m} to mtemp2
       call convert_m_to_umac(mla,rhotot_fc_new,mtemp,umac_new,.false.)
       call mk_advective_m_fluxdiv(mla,umac_new,mtemp,mtemp2,dx, &
                                   the_bc_tower%bc_tower_array)

       p0_update_new = (Sbar_new + Scorrbar)/alphabar_new

       ! update pressure
       p0_new = p0_old + 0.5d0*dt*(p0_update_old + p0_update_new)

       if (parallel_IOProcessor()) then
          print*,'p0_old,new',p0_old,p0_new
       end if

       ! mass_update_new = [-div(rho_i*v) + div(F)]^{n+1,m}
       do n=1,nlevs
          call multifab_copy_c(mass_update_new(n),1,mass_fluxdiv_new(n),1,nspecies,0)
       end do
       call mk_advective_s_fluxdiv(mla,umac_new,rho_fc_new,mass_update_new,dx,1,nspecies)

       ! rhoh_update_new = [-div(rhoh*v) + p0_update + div(Q) + sum(div(hk*Fk)) + rho*Hext]^{n+1,m}
       do n=1,nlevs
          call multifab_copy_c(rhoh_update_new(n),1,rhoh_fluxdiv_new(n),1,1,0)
          call multifab_plus_plus_s_c(rhoh_update_new(n),1,p0_update_new,1,0)
       end do
       call mk_advective_s_fluxdiv(mla,umac_new,rhoh_fc_new,rhoh_update_new,dx,1,1)

       ! compute rho_i^{n+1,m+1}
       do n=1,nlevs
          call multifab_saxpy_5(rho_new(n),1.d0,rho_old(n),0.5d0*dt,mass_update_old(n))
          call multifab_saxpy_3(rho_new(n),0.5d0*dt,mass_update_new(n))
       end do

       ! compute rhotot^{n+1,m+1} and fill ghost cells
       call compute_rhotot(mla,rho_new,rhotot_new)
       do n=1,nlevs
          call multifab_fill_boundary(rhotot_new(n))
          call multifab_physbc(rhotot_new(n),1,scal_bc_comp,1, &
                               the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
       end do

       ! compute w^{n+1,m+1} in valid region and then fill ghost cells
       ! then convert back to rho_i^{n+1,m+1} to fill ghost cells
       call convert_rhoc_to_c(mla,rho_new,rhotot_new,conc_new,.true.)
       call fill_c_ghost_cells(mla,conc_new,dx,the_bc_tower)
       call convert_rhoc_to_c(mla,rho_new,rhotot_new,conc_new,.false.)

       ! compute x^{n+1,m+1} in VALID + GHOST regions
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

          ! cc_solver_alpha = rhotot^{n+1,m+1} c_p^{n+1,m+1,l} / dt
          call compute_cp(mla,cc_solver_alpha,conc_new,Temp_new)
          do n=1,nlevs
             call multifab_mult_mult_c(cc_solver_alpha(n),1,rhotot_new(n),1,1,0)
             call multifab_mult_mult_s_c(cc_solver_alpha(n),1,1.d0/dt,1,0)
          end do

          ! cc_solver_beta = (1/2) lambda^{n+1,m+1,l}
          call average_cc_to_face(nlevs,lambda_new,cc_solver_beta,1,tran_bc_comp,1, &
                                  the_bc_tower%bc_tower_array)          
          do n=1,nlevs
             do i=1,dm
                call multifab_mult_mult_s_c(cc_solver_beta(n,i),1,0.5d0,1,0)
             end do
          end do

          ! initialize deltaT to zero
          do n=1,nlevs
             call setval(deltaT(n),0.d0,all=.true.)
          end do

          ! solve for deltaT
          call ml_cc_solve(mla,deltaT_rhs,deltaT,fine_flx,cc_solver_alpha,cc_solver_beta,dx, &
                           the_bc_tower,temp_bc_comp)

          norm = multifab_norm_l1_c(deltaT(1),1,1)/n_cell
          if (parallel_IOProcessor()) then
             print*,'deltaT_norm',norm
          end if

          ! T^{n+1,m+1,l+1} = T^{n+1,m+1,l} + deltaT
          do n=1,nlevs
             call multifab_plus_plus_c(Temp_new(n),1,deltaT(n),1,1,0)
          end do

          ! fill T ghost cells
          do n=1,nlevs
             call multifab_fill_boundary(Temp_new(n))
             call multifab_physbc(Temp_new(n),1,temp_bc_comp,1,the_bc_tower%bc_tower_array(n), &
                                  dx_in=dx(n,:))
          end do

          ! h^{n+1,m+1,l+1} = h(rhotot^{n+1,m+1},w^{n+1,m+1},T^{n+1,m+1,l+1})
          call compute_h(mla,Temp_new,enth_new,conc_new)
          call convert_rhoh_to_h(mla,rhoh_new,rhotot_new,enth_new,.false.)

          ! compute t^{n+1,m+1,l+1} transport properties (eta,lambda,kappa,chi,zeta)
          call ideal_mixture_transport_wrapper(mla,rhotot_new,Temp_new,p0_new,conc_new, &
                                               molefrac_new,eta_new,lambda_new,kappa_new, &
                                               chi_new,zeta_new)

          ! rhoh_fluxdiv_new = div(Q)^{n+1,m+1,l} + sum(div(hk^{n+1,m+1,l}*Fk^{n+1,m})) + (rho*Hext)^{n+1,m+1}
          call rhoh_fluxdiv_energy(mla,lambda_new,Temp_new,mass_flux_new,rhotot_new, &
                                   rhoh_fluxdiv_new,dx,0.d0,the_bc_tower)

          ! rhoh_update_new = [-div(rhoh*v) + p0_update + div(Q) + sum(div(hk*Fk)) + rho*Hext]^{n+1}
          do n=1,nlevs
             call multifab_copy_c(rhoh_update_new(n),1,rhoh_fluxdiv_new(n),1,1,0)
             call multifab_plus_plus_s_c(rhoh_update_new(n),1,p0_update_new,1,0)
          end do
          call mk_advective_s_fluxdiv(mla,umac_new,rhoh_fc_new,rhoh_update_new,dx,1,1)

       end do  ! end loop l over deltaT_iters

       ! eta^{n+1,m+1,l+1} on nodes (2d) or edges (3d)
       if (dm .eq. 2) then
          call average_cc_to_node(nlevs,eta_new,eta_ed_new(:,1),1,tran_bc_comp,1, &
                                  the_bc_tower%bc_tower_array)
       else if (dm .eq. 3) then
          call average_cc_to_edge(nlevs,eta_new,eta_ed_new,1,tran_bc_comp,1, &
                                  the_bc_tower%bc_tower_array)
       end if

       ! compute rhotot^{n+1,m+1} on faces
       call average_cc_to_face(nlevs,rhotot_new,rhotot_fc_new,1,scal_bc_comp,1, &
                               the_bc_tower%bc_tower_array)

       ! compute rho_i^{n+1,m+1} on faces
       ! first, average c^{n+1,m+1} to faces
       call average_cc_to_face(nlevs,conc_new,rho_fc_new,1,c_bc_comp,nspecies, &
                               the_bc_tower%bc_tower_array)

       ! multiply c^{n+1,m+1} on faces by rhotot^{n+1,m+1} on faces
       do n=1,nlevs
          do i=1,dm
             do comp=1,nspecies
                call multifab_mult_mult_c(rho_fc_new(n,i),comp,rhotot_fc_new(n,i),1,1,0)
             end do
          end do
       end do

       ! compute rhoh^{n+1,m+1} on faces
       ! average h^{n+1,m+1} to faces
       call average_cc_to_face(nlevs,enth_new,rhoh_fc_new,1,h_bc_comp,1, &
                               the_bc_tower%bc_tower_array)

       ! multiply h^{n+1,m+1} on faces by rhotot^{n+1,m+1} on faces
       do n=1,nlevs
          do i=1,dm
             call multifab_mult_mult_c(rhoh_fc_new(n,i),1,rhotot_fc_new(n,i),1,1,0)
          end do
       end do
       
       ! compute mass_flux_new = F^{n+1,m+1} and mass_fluxdiv_new = div(F^{n+1,m+1})
       call mass_fluxdiv_energy(mla,rho_new,rhotot_new,molefrac_new,chi_new,zeta_new, &
                                gradp_baro,Temp_new,mass_fluxdiv_new, &
                                mass_flux_new,dx,the_bc_tower)

       ! compute rhoh_fluxdiv_new = div(Q)^{n+1,m+1} + sum(div(hk*Fk))^{n+1,m+1} + rho_new*Hext^{n+1,m+1}
       call rhoh_fluxdiv_energy(mla,lambda_new,Temp_new,mass_flux_new,rhotot_new, &
                                rhoh_fluxdiv_new,dx,0.d0,the_bc_tower)

       ! compute S^{n+1,m+1} and alpha^{n+1,m+1} (store them in delta_S_new and delta_alpha_new)
       call compute_S_alpha(mla,delta_S_new,delta_alpha_new,mass_fluxdiv_new, &
                            rhoh_fluxdiv_new,conc_new,Temp_new,rhotot_new,p0_new)

       ! compute P_eos^{n+1,m+1}
       call compute_p(mla,rhotot_new,Temp_new,conc_new,Peos)

       ! Scorr = Scorr + (dpdt_factor / p0^{n+1,m+1}) * (Peos^{n+1,m+1} - P0^{n+1,m+1})/dt
       do n=1,nlevs
          call multifab_sub_sub_s_c(Peos(n),1,p0_new,1,0)

          ! debugging statements
          if (.false.) then
             if (k .eq. 1) then
                call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift1")
             else if (k .eq. 2) then
                call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift2")
             else if (k .eq. 3) then
                call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift3")
             else if (k .eq. 4) then
                call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift4")
             else if (k .eq. 5) then
                call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift5")
             else if (k .eq. 6) then
                call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift6")
             else if (k .eq. 7) then
                call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift7")
             else if (k .eq. 8) then
                call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift8")
             else if (k .eq. 9) then
                call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift9")
             end if
          end if

          norm = multifab_norm_l1_c(Peos(n),1,1)/n_cell
          if (parallel_IOProcessor()) then
             print*,'drift_norm',norm
          end if

          call multifab_mult_mult_s_c(Peos(n),1,dpdt_factor/(p0_new*dt),1,0)

          if (k .ge. 2) then
             call multifab_plus_plus_c(Scorr(n),1,Peos(n),1,1,0)
          end if
       end do

       ! split S^{n+1,m+1}, alpha^{n+1,m+1}, and Scorr into average and perturbational pieces
       do n=1,nlevs
          Sbar_new     = multifab_sum_c(delta_S_new(n)    ,1,1) / dble(n_cell)
          alphabar_new = multifab_sum_c(delta_alpha_new(n),1,1) / dble(n_cell)
          Scorrbar     = multifab_sum_c(Scorr(n),1,1) / dble(n_cell)
          call multifab_sub_sub_s_c(delta_S_new(n)    ,1,Sbar_new    ,1,0)
          call multifab_sub_sub_s_c(delta_alpha_new(n),1,alphabar_new,1,0)
          call multifab_copy_c(delta_Scorr(n),1,Scorr(n),1,1,0)
          call multifab_sub_sub_s_c(delta_Scorr(n),1,Scorrbar,1,0)
       end do

       ! compute gmres_rhs_p = delta_S_new + delta_Scorr
       !                       - delta_alpha_new * (Sbar_new + Scorrbar)/alphabar_new
       do n=1,nlevs
          call multifab_copy_c(gmres_rhs_p(n),1,delta_alpha_new(n),1,1,0)
          call multifab_mult_mult_s_c(gmres_rhs_p(n),1,-p0_update_new,1,0)
          call multifab_plus_plus_c(gmres_rhs_p(n),1,delta_S_new(n),1,1,0)
          call multifab_plus_plus_c(gmres_rhs_p(n),1,delta_Scorr(n),1,1,0)
       end do

       ! construct gmres_rhs_p = div(vbar^n) - gmres_rhs_p
       ! first multiply gmres_rhs_p by -1
       do n=1,nlevs
          call multifab_mult_mult_s_c(gmres_rhs_p(n),1,-1.d0,1,0)
       end do

       ! add -div(rho*v*v)^n to mtemp2
       call convert_m_to_umac(mla,rhotot_fc_old,mtemp,umac_old,.false.)
       call mk_advective_m_fluxdiv(mla,umac_old,mtemp,mtemp2,dx, &
                                   the_bc_tower%bc_tower_array)

       ! overwrite umac_new with vbar^n
       ! FIXME: vbar will need boundary conditions at t^{n+1} but for now this is fine
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(umac_new(n,i),1,umac_old(n,i),1,1,umac_new(n,i)%ng)
          end do
       end do

       ! add div(vbar) to gmres_rhs_p
       call compute_div(mla,umac_new,gmres_rhs_p,dx,1,1,1,increment_in=.true.)

       ! construct gmres_rhs_v
       ! set gmres_rhs_v to rhotot^n v^n
       call convert_m_to_umac(mla,rhotot_fc_old,mtemp,umac_old,.false.)
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(gmres_rhs_v(n,i),1,mtemp(n,i),1,1,0)
          end do
       end do

       ! subtract mtemp = rhotot^{n+1,m+1} vbar^n from gmres_rhs_v and then multiply by 1/dt
       call convert_m_to_umac(mla,rhotot_fc_new,mtemp,umac_new,.false.)
       do n=1,nlevs
          do i=1,dm
             call multifab_sub_sub_c(gmres_rhs_v(n,i),1,mtemp(n,i),1,1,0)
             call multifab_mult_mult_s_c(gmres_rhs_v(n,i),1,1.d0/dt,1,0)
          end do
       end do

       ! subtract grad(pi^n) from gmres_rhs_v
       do n=1,nlevs
          do i=1,dm
             call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradpi(n,i),1,1,0)
          end do
       end do

       ! add (1/2)[-div(rho*v*v)^n - div(rho*v*v)^{n+1,m}] to gmres_rhs_v
       do n=1,nlevs
          do i=1,dm
             call multifab_saxpy_3(gmres_rhs_v(n,i),0.5d0,mtemp2(n,i))
          end do
       end do

       ! add (1/2) A_0^n v^n to gmres_rhs_v
       do n=1,nlevs
          do i=1,dm
             call setval(mtemp(n,i),0.d0,all=.true.)
          end do
       end do
       call diffusive_m_fluxdiv(mla,mtemp,umac_old,eta_old,eta_ed_old,kappa_old,dx, &
                                the_bc_tower%bc_tower_array)
       do n=1,nlevs
          do i=1,dm
             call multifab_saxpy_3(gmres_rhs_v(n,i),0.5d0,mtemp(n,i))
          end do
       end do

       ! add (1/2) A_0^{n+1,m+1} vbar^n to gmres_rhs_v
       do n=1,nlevs
          do i=1,dm
             call setval(mtemp(n,i),0.d0,all=.true.)
          end do
       end do
       call diffusive_m_fluxdiv(mla,mtemp,umac_new,eta_new,eta_ed_new,kappa_new,dx, &
                                the_bc_tower%bc_tower_array)
       do n=1,nlevs
          do i=1,dm
             call multifab_saxpy_3(gmres_rhs_v(n,i),0.5d0,mtemp(n,i))
          end do
       end do

       ! add (1/2)(rhotot^n + rhotot^{n+1,m+1})g to gmres_rhs_v
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

       ! multiply (eta,kappa)^{n+1,m+1} by 1/2 to put in proper form for gmres solve
       do n=1,nlevs
          call multifab_mult_mult_s_c(eta_new(n),1,1.d0/2.d0,1,eta_new(n)%ng)
          call multifab_mult_mult_s_c(kappa_new(n),1,1.d0/2.d0,1,kappa_new(n)%ng)
          do i=1,size(eta_ed_new,dim=2)
             call multifab_mult_mult_s_c(eta_ed_new(n,i),1,1.d0/2.d0,1,eta_ed_new(n,i)%ng)
          end do
       end do

       ! call the Stokes solver
       call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dpi,rhotot_fc_new, &
                  eta_new,eta_ed_new,kappa_new,theta_alpha,norm_pre_rhs)

       ! restore (eta,kappa)^{n+1,m+1}
       do n=1,nlevs
          call multifab_mult_mult_s_c(eta_new(n),1,2.d0,1,eta_new(n)%ng)
          call multifab_mult_mult_s_c(kappa_new(n),1,2.d0,1,kappa_new(n)%ng)
          do i=1,size(eta_ed_new,dim=2)
             call multifab_mult_mult_s_c(eta_ed_new(n,i),1,2.d0,1,eta_ed_new(n,i)%ng)
          end do
       end do

       ! increment velocity
       ! compute v^{n+1,m+1} = vbar^n + dumac
       ! note: no need to compute pi^{n+1,m+1} until the final iteration
       do n=1,nlevs
          do i=1,dm
             call multifab_plus_plus_c(umac_new(n,i),1,dumac(n,i),1,1,0)
          end do
       end do

       do n=1,nlevs
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

    end do  ! end loop k over dpdt_iters

    do n=1,nlevs
       ! compute pi^{n+1}= pi^n + dpi and update ghost cells
       call multifab_plus_plus_c(pi(n),1,dpi(n),1,1,0)
       call multifab_fill_boundary(pi(n))
       call multifab_physbc(pi(n),1,pres_bc_comp,1,the_bc_tower%bc_tower_array(n), &
                            dx_in=dx(n,:))
    end do

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
       call multifab_destroy(delta_S_new(n))
       call multifab_destroy(delta_alpha_old(n))
       call multifab_destroy(delta_alpha_new(n))
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
          call multifab_destroy(gmres_rhs_v(n,i))
       end do
       call multifab_destroy(gmres_rhs_p(n))
       do i=1,dm
          call multifab_destroy(gradpi(n,i))
          call multifab_destroy(mtemp(n,i))
          call multifab_destroy(mtemp2(n,i))
          call multifab_destroy(dumac(n,i))
       end do
       call multifab_destroy(dpi(n))
       do i=1,size(eta_ed_old,dim=2)
          call multifab_destroy(eta_ed_old(n,i))
          call multifab_destroy(eta_ed_new(n,i))
       end do
    end do

    do n=2,nlevs
       call bndry_reg_destroy(fine_flx(n))
    end do

    deallocate(eta_ed_old,eta_ed_new)

  end subroutine advance_timestep

end module advance_timestep_module
