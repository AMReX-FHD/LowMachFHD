module advance_timestep_inertial_module

  use ml_layout_module
  use define_bc_module
  use bc_module
  use convert_stag_module
  use mk_advective_s_fluxdiv_module
  use diffusive_m_fluxdiv_module
  use compute_mass_fluxdiv_module
  use project_onto_eos_module
  use compute_HSE_pres_module
  use convert_m_to_umac_module
  use convert_rhoc_to_c_module
  use mk_advective_m_fluxdiv_module
  use bds_module
  use gmres_module
  use div_and_grad_module
  use mk_grav_force_module
  use compute_mixture_properties_module
  use mass_flux_utilities_module
  use multifab_physbc_module
  use multifab_physbc_stag_module
  use zero_edgeval_module
  use fill_rho_ghost_cells_module
  use bl_rng_module
  use bl_random_module
  use probin_common_module, only: advection_type, grav, variance_coef_mass, &
                                  variance_coef_mom, barodiffusion_type, &
                                  use_bl_rng, nspecies
  use probin_gmres_module, only: gmres_abs_tol, gmres_rel_tol
  use probin_charged_module, only: use_charged_fluid
  use probin_chemistry_module, only: nreactions

  implicit none

  private

  public :: advance_timestep_inertial

contains

  subroutine advance_timestep_inertial(mla,umac,rho_old,rho_new,rhotot_old,rhotot_new, &
                                       rhoh_old,rhoh_new,p0_old,p0_new, &
                                       gradp_baro,pi,eta,eta_ed,kappa,Temp,Temp_ed, &
                                       diff_mass_fluxdiv,stoch_mass_fluxdiv,stoch_mass_flux, &
                                       dx,dt,time,the_bc_tower,istep)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(inout) :: rho_old(:)
    type(multifab) , intent(inout) :: rho_new(:)
    type(multifab) , intent(inout) :: rhotot_old(:)
    type(multifab) , intent(inout) :: rhotot_new(:)
    type(multifab) , intent(inout) :: rhoh_old(:)
    type(multifab) , intent(inout) :: rhoh_new(:)
    type(multifab) , intent(inout) :: gradp_baro(:,:)
    type(multifab) , intent(inout) :: pi(:)
    ! eta and kappa need to enter consistent with old and leave consistent with new
    type(multifab) , intent(inout) :: eta(:)
    type(multifab) , intent(inout) :: eta_ed(:,:) ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(inout) :: kappa(:)
    type(multifab) , intent(inout) :: Temp(:)
    type(multifab) , intent(inout) :: Temp_ed(:,:) ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(inout) :: diff_mass_fluxdiv(:)
    type(multifab) , intent(inout) :: stoch_mass_fluxdiv(:)
    type(multifab) , intent(inout) :: stoch_mass_flux(:,:)
    real(kind=dp_t), intent(inout) :: p0_old,p0_new
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt,time
    type(bc_tower) , intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: istep

    ! local
    type(multifab) ::  rho_update(mla%nlevel)
    type(multifab) ::   bds_force(mla%nlevel)
    type(multifab) :: gmres_rhs_p(mla%nlevel)
    type(multifab) ::         dpi(mla%nlevel)
    type(multifab) ::  rho_nd_old(mla%nlevel)
    type(multifab) ::     rho_tmp(mla%nlevel)
    type(multifab) ::      p_baro(mla%nlevel)

    type(multifab) ::              mold(mla%nlevel,mla%dim)
    type(multifab) ::             mtemp(mla%nlevel,mla%dim)
    type(multifab) ::   adv_mom_fluxdiv(mla%nlevel,mla%dim)
    type(multifab) ::  diff_mom_fluxdiv(mla%nlevel,mla%dim)
    type(multifab) ::       gmres_rhs_v(mla%nlevel,mla%dim)
    type(multifab) ::             dumac(mla%nlevel,mla%dim)
    type(multifab) ::          umac_tmp(mla%nlevel,mla%dim)
    type(multifab) ::     rhotot_fc_old(mla%nlevel,mla%dim)
    type(multifab) ::     rhotot_fc_new(mla%nlevel,mla%dim)
    type(multifab) ::            gradpi(mla%nlevel,mla%dim)
    type(multifab) ::            rho_fc(mla%nlevel,mla%dim)
    type(multifab) ::    diff_mass_flux(mla%nlevel,mla%dim)
    type(multifab) ::   total_mass_flux(mla%nlevel,mla%dim)

    ! energy specific multifabs
    type(multifab) :: rhoh_fc(mla%nlevel,mla%dim)
    
    integer :: i,dm,n,nlevs

    real(kind=dp_t) :: theta_alpha, norm_pre_rhs, gmres_abs_tol_in

    real(kind=dp_t) :: weights(1)

    type(bl_prof_timer), save :: bpt

    if (use_charged_fluid) then
       call bl_error("energy advance_timestep_inertial does not support charges")
    end if

    if (nreactions .gt. 0) then
       call bl_error("energy advance_timestep_inertial does not support reactions")
    end if

    if (variance_coef_mass .ne. 0.d0 .or. variance_coef_mom .ne. 0.d0) then
       call bl_error("energy advance_timestep_inertial does not support fluctuations")
    end if

    if (advection_type .ge. 1) then
       call bl_error("energy advance_timestep_inertial does not support bds advection")
    end if

    if (barodiffusion_type .ge. 1) then
       call bl_error("energy advance_timestep_inertial does not support barodiffusion")
    end if

    if (is_nonisothermal) then
       call bl_error("energy advance_timestep_inertial does not support is_nonisothermal")
    end if

    call build(bpt, "advance_timestep_inertial")

    weights(1) = 1.d0

    nlevs = mla%nlevel
    dm = mla%dim

    theta_alpha = 1.d0/dt
    
    do n=1,nlevs
       call multifab_build( rho_update(n),mla%la(n),nspecies,0)
       call multifab_build(  bds_force(n),mla%la(n),nspecies,1)
       call multifab_build(gmres_rhs_p(n),mla%la(n),1       ,0)
       call multifab_build(        dpi(n),mla%la(n),1       ,1)
       call multifab_build(     p_baro(n),mla%la(n),1       ,1)
       do i=1,dm
          call multifab_build_edge(            mold(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(           mtemp(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge( adv_mom_fluxdiv(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(diff_mom_fluxdiv(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(     gmres_rhs_v(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(           dumac(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(        umac_tmp(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(          gradpi(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(   rhotot_fc_old(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(   rhotot_fc_new(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(          rho_fc(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge(  diff_mass_flux(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge( total_mass_flux(n,i),mla%la(n),nspecies,0,i)
       end do
    end do

    ! energy specific
    do n=1,nlevs
       do i=1,dm
          call multifab_build_edge(rhoh_fc(n,i),mla%la(n),1,0,i)
       end do
    end do

    ! make copies of old quantities
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(umac_tmp(n,i),1,umac(n,i),1,1,1)
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 1 - Calculate Predictor Diffusive and Stochastic Fluxes
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! diff/stoch_mass_fluxdiv already contain F_i
    ! this was already done in Step 0 (initialization) or Step 6 from the previous time step

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 2 - Predictor Euler Step
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! average rho_old and rhotot_old to faces
    call average_cc_to_face(nlevs,   rho_old,   rho_fc    ,1,   c_bc_comp,nspecies,the_bc_tower%bc_tower_array)
    call average_cc_to_face(nlevs,rhotot_old,rhotot_fc_old,1,scal_bc_comp,       1,the_bc_tower%bc_tower_array)

    ! energy
    call average_cc_to_face(nlevs,rhoh_old  ,rhoh_fc      ,1,h_bc_comp   ,1       ,the_bc_tower%bc_tower_array)

    ! add D^n and St^n to rho_update
    do n=1,nlevs
       call multifab_copy_c(rho_update(n),1,diff_mass_fluxdiv(n),1,nspecies,0)
    end do

    ! add A^n to rho_update
    ! compute A^n for scalars using centered advection
    call mk_advective_s_fluxdiv(mla,umac,rho_fc,rho_update,.true.,dx,1,nspecies)

    ! energy - FIXME



    ! set rho_new = rho_old + dt * (A^n + D^n + St^n)
    do n=1,nlevs
       call multifab_saxpy_4(rho_new(n),rho_old(n),dt,rho_update(n))
    end do

    ! compute rhotot from rho in VALID REGION
    call compute_rhotot(mla,rho_new,rhotot_new)

    ! fill rho and rhotot ghost cells
    call fill_rho_rhotot_ghost(mla,rho_new,rhotot_new,dx,the_bc_tower)

    ! average rho_new and rhotot_new to faces
    call average_cc_to_face(nlevs,   rho_new,   rho_fc    ,1,   c_bc_comp,nspecies,the_bc_tower%bc_tower_array)
    call average_cc_to_face(nlevs,rhotot_new,rhotot_fc_new,1,scal_bc_comp,       1,the_bc_tower%bc_tower_array)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 3 - Calculate Corrector Diffusive and Stochastic Fluxes
    ! Step 4 - Predictor Crank-Nicolson Step
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute mold
    call convert_m_to_umac(mla,rhotot_fc_old,mold,umac,.false.)

    ! build up rhs_v for gmres solve: first set gmres_rhs_v to mold/dt
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_5(gmres_rhs_v(n,i),1.d0/dt,mold(n,i),0.d0,mold(n,i))
       end do
    end do

    ! compute grad pi^n
    call compute_grad(mla,pi,gradpi,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

    ! subtract grad pi^n from gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradpi(n,i),1,1,0)
       end do
    end do

    ! compute adv_mom_fluxdiv = A^n for momentum
    call mk_advective_m_fluxdiv(mla,umac,mold,adv_mom_fluxdiv,.false.,dx, &
                                the_bc_tower%bc_tower_array)

    ! add A^n for momentum to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,adv_mom_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! compute diff_mom_fluxdiv = A_0^n v^n
    call diffusive_m_fluxdiv(mla,diff_mom_fluxdiv,.false.,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    ! add (1/2) A_0^n v^n to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3(gmres_rhs_v(n,i),0.5d0,diff_mom_fluxdiv(n,i))
       end do
    end do

    ! add rho^n*g to gmres_rhs_v
    if (any(grav(1:dm) .ne. 0.d0)) then
       call mk_grav_force(mla,gmres_rhs_v,.true.,rhotot_fc_old,rhotot_fc_old,the_bc_tower)
    end if

    ! compute diffusive and stochastic mass fluxes
    ! this computes "-F = rho W chi [Gamma grad x... ]"
    call compute_mass_fluxdiv(mla,rho_new,rhotot_new,gradp_baro,Temp, &
                              diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                              diff_mass_flux,stoch_mass_flux, &
                              dt,time,dx,weights,the_bc_tower)

    ! energy - FIXME call diffusive_rhoh_fluxdiv


    ! compute gmres_rhs_p (we will later add divu)
    do n=1,nlevs
       call setval(gmres_rhs_p(n),0.d0,all=.true.)
       ! FIXME


    end do

    ! modify umac to respect the boundary conditions we want after the next gmres solve
    ! thus when we add A_0^n vbar^n to gmres_rhs_v and add div vbar^n to gmres_rhs_p we
    ! are automatically putting the system in delta form WITH homogeneous boundary conditions
    do n=1,nlevs
       do i=1,dm
          ! set normal velocity on physical domain boundaries
          call multifab_physbc_domainvel(umac(n,i),vel_bc_comp+i-1, &
                                         the_bc_tower%bc_tower_array(n), &
                                         dx(n,:),vel_bc_n(n,:))
          ! set transverse velocity behind physical boundaries
          call multifab_physbc_macvel(umac(n,i),vel_bc_comp+i-1, &
                                      the_bc_tower%bc_tower_array(n), &
                                      dx(n,:),vel_bc_t(n,:))
          ! fill periodic and interior ghost cells
          call multifab_fill_boundary(umac(n,i))
       end do
    end do

    ! compute mtemp = rho^{*,n+1} * vbar^n
    call convert_m_to_umac(mla,rhotot_fc_new,mtemp,umac,.false.)

    ! subtract rho^{*,n+1} * vbar^n / dt from gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3(gmres_rhs_v(n,i),-1.d0/dt,mtemp(n,i))
       end do
    end do

    ! compute mtemp = A_0^n vbar^n
    call diffusive_m_fluxdiv(mla,mtemp,.false.,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    ! add (1/2) A_0^n vbar^n to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3(gmres_rhs_v(n,i),0.5d0,mtemp(n,i))
       end do
    end do

    ! set physical boundary values to zero
    do n=1,nlevs
       call zero_edgeval_physical(gmres_rhs_v(n,:),1,1,the_bc_tower%bc_tower_array(n))
    end do

    ! set rho_update to F^{*,n+1} = div(rho*chi grad c)^{*,n+1} + div(Psi^n)
    ! it is used in Step 5 below
    do n=1,nlevs
       call multifab_copy_c(rho_update(n),1, diff_mass_fluxdiv(n),1,nspecies)
    end do

    ! compute div vbar^n and add to gmres_rhs_p
    ! now gmres_rhs_p = div vbar^n - S^{*,n+1}
    ! the sign convention is correct since we solve -div(delta v) = gmres_rhs_p
    call compute_div(mla,umac,gmres_rhs_p,dx,1,1,1,increment_in=.true.)

    ! multiply eta and kappa by 1/2 to put in proper form for gmres solve
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,1.d0/2.d0,1,eta(n)%ng)
       call multifab_mult_mult_s_c(kappa(n),1,1.d0/2.d0,1,kappa(n)%ng)
       do i=1,size(eta_ed,dim=2)
          call multifab_mult_mult_s_c(eta_ed(n,i),1,1.d0/2.d0,1,eta_ed(n,i)%ng)
       end do
    end do

    ! set the initial guess to zero
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(dumac(n,i),0.d0,all=.true.)
       end do
       call multifab_setval(dpi(n),0.d0,all=.true.)
    end do

    gmres_abs_tol_in = gmres_abs_tol ! Save this 

    ! This relies entirely on relative tolerance and can fail if the rhs is roundoff error only:
    ! gmres_abs_tol = 0.d0 ! It is better to set gmres_abs_tol in namelist to a sensible value

    ! call gmres to compute delta v and delta pi
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dpi,rhotot_fc_new, &
               eta,eta_ed,kappa,theta_alpha,norm_pre_rhs)

    ! for the corrector gmres solve we want the stopping criteria based on the
    ! norm of the preconditioned rhs from the predictor gmres solve.  otherwise
    ! for cases where du in the corrector should be small the gmres stalls
    gmres_abs_tol = max(gmres_abs_tol_in, norm_pre_rhs*gmres_rel_tol)

    ! restore eta and kappa
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,2.d0,1,eta(n)%ng)
       call multifab_mult_mult_s_c(kappa(n),1,2.d0,1,kappa(n)%ng)
       do i=1,size(eta_ed,dim=2)
          call multifab_mult_mult_s_c(eta_ed(n,i),1,2.d0,1,eta_ed(n,i)%ng)
       end do
    end do

    ! compute v^{*,n+1} = v^n + dumac
    ! compute pi^{*,n+1}= pi^n + dpi
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(umac(n,i),1,dumac(n,i),1,1,0)
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
          call multifab_physbc_domainvel(umac(n,i),vel_bc_comp+i-1, &
                                         the_bc_tower%bc_tower_array(n), &
                                         dx(n,:),vel_bc_n(n,:))
          ! set transverse velocity behind physical boundaries
          call multifab_physbc_macvel(umac(n,i),vel_bc_comp+i-1, &
                                      the_bc_tower%bc_tower_array(n), &
                                      dx(n,:),vel_bc_t(n,:))
          ! fill periodic and interior ghost cells
          call multifab_fill_boundary(umac(n,i))
       end do
    end do

    ! convert v^{*,n+1} to rho^{*,n+1}v^{*,n+1} in valid and ghost region
    ! now mnew has properly filled ghost cells
    call convert_m_to_umac(mla,rhotot_fc_new,mtemp,umac,.false.)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 5 - Trapezoidal Scalar Corrector
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! rho_update already contains D^{*,n+1} + St^{*,n+1} for rho from above
    ! add A^{*,n+1} for rho to rho_update
    ! compute A^{*,n+1} for scalars using centered advection
    call mk_advective_s_fluxdiv(mla,umac,rho_fc,rho_update,.true.,dx,1,nspecies)

    ! snew = s^{n+1} 
    !      = (1/2)*(s^n + s^{*,n+1} + dt*(A^{*,n+1} + D^{*,n+1} + St^{*,n+1}))
    do n=1,nlevs
       call multifab_plus_plus_c(rho_new(n),1,rho_old(n),1,nspecies,0)
       call multifab_saxpy_3(rho_new(n),dt,rho_update(n))
       call multifab_mult_mult_s_c(rho_new(n),1,0.5d0,nspecies,0)
    end do

    ! compute rhotot from rho in VALID REGION
    call compute_rhotot(mla,rho_new,rhotot_new)

    ! fill rho and rhotot ghost cells
    call fill_rho_rhotot_ghost(mla,rho_new,rhotot_new,dx,the_bc_tower)

    ! average rho_new and rhotot_new to faces
    call average_cc_to_face(nlevs,   rho_new,   rho_fc    ,1,   c_bc_comp,nspecies,the_bc_tower%bc_tower_array)
    call average_cc_to_face(nlevs,rhotot_new,rhotot_fc_new,1,scal_bc_comp,       1,the_bc_tower%bc_tower_array)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 6 - Calculate Diffusive and Stochastic Fluxes
    ! Step 7 - Corrector Crank-Nicolson Step
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! build up rhs_v for gmres solve: first set gmres_rhs_v to mold/dt
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_5(gmres_rhs_v(n,i),1.d0/dt,mold(n,i),0.d0,mold(n,i))
       end do
    end do

    ! compute grad pi^{*,n+1}
    call compute_grad(mla,pi,gradpi,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

    ! subtract grad pi^{*,n+1} from gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradpi(n,i),1,1,0)
       end do
    end do

    ! adv_mom_fluxdiv already contains A^n for momentum
    ! add A^{*,n+1} = -rho^{*,n+1} v^{*,n+1} v^{*,n+1} for momentum to adv_mom_fluxdiv
    call mk_advective_m_fluxdiv(mla,umac,mtemp,adv_mom_fluxdiv,.true.,dx, &
                                the_bc_tower%bc_tower_array)

    ! add (1/2) adv_mom_fluxdiv to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3(gmres_rhs_v(n,i),0.5d0,adv_mom_fluxdiv(n,i))
       end do
    end do

    ! add (1/2) A_0^n v^n to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3(gmres_rhs_v(n,i),0.5d0,diff_mom_fluxdiv(n,i))
       end do
    end do

    ! add gravity term
    if (any(grav(1:dm) .ne. 0.d0)) then
       call mk_grav_force(mla,gmres_rhs_v,.true.,rhotot_fc_old,rhotot_fc_new,the_bc_tower)
    end if

    ! compute diffusive and stochastic mass fluxes
    ! this computes "-F = rho W chi [Gamma grad x... ]"
    call compute_mass_fluxdiv(mla,rho_new,rhotot_new,gradp_baro,Temp, &
                              diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                              diff_mass_flux,stoch_mass_flux, &
                              dt,time,dx,weights,the_bc_tower)

    ! compute gmres_rhs_p (we will later add divu)
    do n=1,nlevs
       call setval(gmres_rhs_p(n),0.d0,all=.true.)
       ! FIXME


    end do

    ! modify umac to respect the boundary conditions we want after the next gmres solve
    ! thus when we add A_0^{n+1} vbar^{*,n+1} to gmres_rhs_v and add div vbar^{*,n+1} to gmres_rhs_p we
    ! are automatically putting the system in delta form WITH homogeneous boundary conditions
    do n=1,nlevs
       do i=1,dm
          ! set normal velocity on physical domain boundaries
          call multifab_physbc_domainvel(umac(n,i),vel_bc_comp+i-1, &
                                         the_bc_tower%bc_tower_array(n), &
                                         dx(n,:),vel_bc_n(n,:))
          ! set transverse velocity behind physical boundaries
          call multifab_physbc_macvel(umac(n,i),vel_bc_comp+i-1, &
                                      the_bc_tower%bc_tower_array(n), &
                                      dx(n,:),vel_bc_t(n,:))
          ! fill periodic and interior ghost cells
          call multifab_fill_boundary(umac(n,i))
       end do
    end do

   ! compute mtemp = rho^{n+1} * vbar^{*,n+1}
   call convert_m_to_umac(mla,rhotot_fc_new,mtemp,umac,.false.)

   ! subtract rho^{n+1} * vbar^{*,n+1} / dt from gmres_rhs_v
   do n=1,nlevs
      do i=1,dm
         call multifab_saxpy_3(gmres_rhs_v(n,i),-1.d0/dt,mtemp(n,i))
      end do
   end do

    ! set diff_mom_fluxdiv = A_0^{n+1} vbar^{n+1,*}
    call diffusive_m_fluxdiv(mla,diff_mom_fluxdiv,.false.,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    ! add (1/2) A_0^{n+1} vbar^{n+1,*} to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3(gmres_rhs_v(n,i),0.5d0,diff_mom_fluxdiv(n,i))
       end do
    end do

    ! set physical boundary values to zero
    do n=1,nlevs
       call zero_edgeval_physical(gmres_rhs_v(n,:),1,1,the_bc_tower%bc_tower_array(n))
    end do

    ! compute div(vbar^{n+1,*}) and add to gmres_rhs_p
    ! now gmres_rhs_p = div(vbar^{n+1,*}) - S^{n+1}
    ! the sign convention is correct since we solve -div(delta v) = gmres_rhs_p
    call compute_div(mla,umac,gmres_rhs_p,dx,1,1,1,increment_in=.true.)

    ! multiply eta and kappa by 1/2 to put in proper form for gmres solve
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,1.d0/2.d0,1,eta(n)%ng)
       call multifab_mult_mult_s_c(kappa(n),1,1.d0/2.d0,1,eta(n)%ng)
       do i=1,size(eta_ed,dim=2)
          call multifab_mult_mult_s_c(eta_ed(n,i),1,1.d0/2.d0,1,eta_ed(n,i)%ng)
       end do
    end do

    ! set the initial guess to zero
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(dumac(n,i),0.d0,all=.true.)
       end do
          call multifab_setval(dpi(n),0.d0,all=.true.)
    end do

    ! call gmres to compute delta v and delta pi
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dpi,rhotot_fc_new, &
               eta,eta_ed,kappa,theta_alpha)
                              
    gmres_abs_tol = gmres_abs_tol_in ! Restore the desired tolerance   

    ! restore eta and kappa
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,2.d0,1,eta(n)%ng)
       call multifab_mult_mult_s_c(kappa(n),1,2.d0,1,kappa(n)%ng)
       do i=1,size(eta_ed,dim=2)
          call multifab_mult_mult_s_c(eta_ed(n,i),1,2.d0,1,eta_ed(n,i)%ng)
       end do
    end do

    ! compute v^{n+1} = v^{n+1,*} + dumac
    ! compute pi^{n+1} = pi^{n+1,*} + dpi
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(umac(n,i),1,dumac(n,i),1,1,0)
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
          call multifab_physbc_domainvel(umac(n,i),vel_bc_comp+i-1, &
                                         the_bc_tower%bc_tower_array(n), &
                                         dx(n,:),vel_bc_n(n,:))
          ! set transverse velocity behind physical boundaries
          call multifab_physbc_macvel(umac(n,i),vel_bc_comp+i-1, &
                                      the_bc_tower%bc_tower_array(n), &
                                      dx(n,:),vel_bc_t(n,:))
          ! fill periodic and interior ghost cells
          call multifab_fill_boundary(umac(n,i))
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End Time-Advancement
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do n=1,nlevs
       call multifab_destroy(rho_update(n))
       call multifab_destroy(bds_force(n))
       call multifab_destroy(gmres_rhs_p(n))
       call multifab_destroy(dpi(n))
       call multifab_destroy(p_baro(n))
       do i=1,dm
          call multifab_destroy(mold(n,i))
          call multifab_destroy(mtemp(n,i))
          call multifab_destroy(adv_mom_fluxdiv(n,i))
          call multifab_destroy(diff_mom_fluxdiv(n,i))
          call multifab_destroy(gmres_rhs_v(n,i))
          call multifab_destroy(dumac(n,i))
          call multifab_destroy(umac_tmp(n,i))
          call multifab_destroy(rhotot_fc_old(n,i))
          call multifab_destroy(rhotot_fc_new(n,i))
          call multifab_destroy(gradpi(n,i))
          call multifab_destroy(rho_fc(n,i))
          call multifab_destroy(diff_mass_flux(n,i))
          call multifab_destroy(total_mass_flux(n,i))
       end do
    end do

    call destroy(bpt)

  end subroutine advance_timestep_inertial

end module advance_timestep_inertial_module
