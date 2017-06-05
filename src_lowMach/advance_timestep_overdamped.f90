module advance_timestep_overdamped_module

  use ml_layout_module
  use define_bc_module
  use bc_module
  use convert_stag_module
  use mk_advective_s_fluxdiv_module
  use diffusive_m_fluxdiv_module
  use stochastic_m_fluxdiv_module
  use stochastic_mass_fluxdiv_module
  use compute_mass_fluxdiv_module
  use project_onto_eos_module
  use compute_HSE_pres_module
  use convert_rhoc_to_c_module
  use reservoir_bc_fill_module
  use bds_module
  use gmres_module
  use div_and_grad_module
  use mk_grav_force_module
  use compute_mixture_properties_module
  use mass_flux_utilities_module
  use multifab_physbc_module
  use multifab_physbc_stag_module
  use fill_rho_ghost_cells_module
  use bl_rng_module
  use bl_random_module
  use chemical_rates_module
  use probin_common_module, only: advection_type, grav, rhobar, variance_coef_mass, &
                                  variance_coef_mom, barodiffusion_type, restart, &
                                  molmass, nspecies, project_eos_int, use_bl_rng
  use probin_gmres_module, only: gmres_abs_tol, gmres_rel_tol
  use probin_charged_module, only: use_charged_fluid
  use probin_chemistry_module, only: nreactions, use_Poisson_rng
  use probin_multispecies_module, only: midpoint_stoch_mass_flux_type

  implicit none

  private

  public :: advance_timestep_overdamped

contains

  ! eta and kappa can be local temps inside advance_timestep_overdamped
  ! This is consistent with what is done for mass diffusion coefficients
  ! They are local to the wrapper and not really needed outside
  ! Note for future: In general Temp can depend on time so here one should pass
  ! both temperature at the beginning and at the end of the timestep
  subroutine advance_timestep_overdamped(mla,umac,rho_old,rho_new,rhotot_old,rhotot_new, &
                                         gradp_baro,pi,eta,eta_ed,kappa,Temp,Temp_ed, &
                                         diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                         stoch_mass_flux,chem_rate, &
                                         dx,dt,time,the_bc_tower,istep)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(inout) :: rho_old(:)
    type(multifab) , intent(inout) :: rho_new(:)
    type(multifab) , intent(inout) :: rhotot_old(:)
    type(multifab) , intent(inout) :: rhotot_new(:)
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
    type(multifab) , intent(inout) :: chem_rate(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt,time
    type(bc_tower) , intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: istep

    ! local
    type(multifab) ::     rho_update(mla%nlevel)
    type(multifab) ::      bds_force(mla%nlevel)
    type(multifab) ::    gmres_rhs_p(mla%nlevel)
    type(multifab) ::            dpi(mla%nlevel)
    type(multifab) ::     rho_nd_old(mla%nlevel)
    type(multifab) ::        rho_tmp(mla%nlevel)
    type(multifab) ::         p_baro(mla%nlevel)
    type(multifab) :: chem_rate_temp(mla%nlevel)
    type(multifab) ::          n_old(mla%nlevel)
    type(multifab) ::          n_new(mla%nlevel)

    type(multifab) ::         gmres_rhs_v(mla%nlevel,mla%dim)
    type(multifab) ::               dumac(mla%nlevel,mla%dim)
    type(multifab) ::           rhotot_fc(mla%nlevel,mla%dim)
    type(multifab) ::              gradpi(mla%nlevel,mla%dim)
    type(multifab) ::              rho_fc(mla%nlevel,mla%dim)
    type(multifab) ::      diff_mass_flux(mla%nlevel,mla%dim)
    type(multifab) ::     total_mass_flux(mla%nlevel,mla%dim)

    type(multifab) :: stoch_mass_fluxdiv_old(mla%nlevel)
    type(multifab) ::    stoch_mass_flux_old(mla%nlevel,mla%dim)


    integer :: i,dm,n,nlevs

    real(kind=dp_t) :: theta_alpha, norm_pre_rhs, gmres_abs_tol_in
    real(kind=dp_t) :: weights(2)

    real(kind=dp_t), parameter :: mattingly_lin_comb_coef(1:2) = (/-1.d0, 2.d0/)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "advance_timestep_overdamped")

    if (use_charged_fluid) then
       call bl_error('advance_timestep_overdamped does not support charges yet')
    end if

    if (nreactions > 0) then
       if (use_Poisson_rng .eq. 2) then
          call bl_error('Error: currently use_Poisson_rng=2 not allowed for algorith_type=2 and nreactions>0')
       end if
    end if

    nlevs = mla%nlevel
    dm = mla%dim

    theta_alpha = 0.d0
    
    do n=1,nlevs
       call multifab_build( rho_update(n),mla%la(n),nspecies,0)
       call multifab_build(  bds_force(n),mla%la(n),nspecies,1)
       call multifab_build(gmres_rhs_p(n),mla%la(n),1       ,0)
       call multifab_build(        dpi(n),mla%la(n),1       ,1)
       call multifab_build(     p_baro(n),mla%la(n),1       ,1)
       do i=1,dm
          call multifab_build_edge(        gmres_rhs_v(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(              dumac(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(             gradpi(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(          rhotot_fc(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(             rho_fc(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge(     diff_mass_flux(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge(    total_mass_flux(n,i),mla%la(n),nspecies,0,i)
       end do
    end do

    if (variance_coef_mass .ne. 0.d0) then
       do n=1,nlevs
          call multifab_build(stoch_mass_fluxdiv_old(n),mla%la(n),nspecies,0)
          do i=1,dm
             call multifab_build_edge(stoch_mass_flux_old(n,i),mla%la(n),nspecies,0,i)
          end do
       end do
    end if

    if (nreactions > 0) then
       do n=1,nlevs
          call multifab_build(chem_rate_temp(n),mla%la(n),nspecies,0)
          call multifab_build(         n_old(n),mla%la(n),nspecies,0)
          call multifab_build(         n_new(n),mla%la(n),nspecies,0)
       end do
    end if

    do n=1,nlevs
       do i=1,dm
          call setval(dumac(n,i),0.d0,all=.true.)
       end do
    end do

    ! average rho and rhotot to faces
    call average_cc_to_face(nlevs,   rho_old,   rho_fc,1,   c_bc_comp,nspecies,the_bc_tower%bc_tower_array)
    call average_cc_to_face(nlevs,rhotot_old,rhotot_fc,1,scal_bc_comp,       1,the_bc_tower%bc_tower_array)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 1 - Predictor Stochastic/Diffusive Fluxes
    ! Step 2 - Predictor Stokes Solve
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    weights(1) = 1.d0
    weights(2) = 0.d0

    ! fill the stochastic multifabs with a new set of random numbers
    ! if this is the first step after initialization or restart then
    ! we already have random numbers from initialization
    if (variance_coef_mass .ne. 0.d0) then
       call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)
       if (use_bl_rng) then
          ! save random state for restart
          call bl_rng_copy_engine(rng_eng_diffusion_chk,rng_eng_diffusion)
       end if
    end if

    ! compute grad pi^{n-1/2}
    call compute_grad(mla,pi,gradpi,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

    if (barodiffusion_type .eq. 2) then
       ! barodiffusion uses lagged grad(pi)
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(gradp_baro(n,i),1,gradpi(n,i),1,1,0)
          end do
       end do
    else if (barodiffusion_type .eq. 3) then
       ! compute p0 from rho0*g
       call compute_HSE_pres(mla,rhotot_old,p_baro,dx,the_bc_tower)
       call compute_grad(mla,p_baro,gradp_baro,dx,1,pres_bc_comp,1,1, &
                         the_bc_tower%bc_tower_array)
    end if

    ! set gmres_rhs_v to -grad pi^{n-1/2}
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_5(gmres_rhs_v(n,i),-1.d0,gradpi(n,i),0.d0,gradpi(n,i))
       end do
    end do

    ! add div(Sigma^(1)) to gmres_rhs_v
    if (variance_coef_mom .ne. 0.d0) then
       ! fill the stochastic multifabs with a new set of random numbers
       call fill_m_stochastic(mla)
       call stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,gmres_rhs_v,.true., &
                                 eta,eta_ed,Temp,Temp_ed,dx,0.5d0*dt,weights)
    end if

    ! add rho^n*g to gmres_rhs_v
    if (any(grav(1:dm) .ne. 0.d0)) then
       call mk_grav_force(mla,gmres_rhs_v,.true.,rhotot_fc,rhotot_fc,the_bc_tower)
    end if

    ! reset inhomogeneous bc condition to deal with reservoirs
    call set_inhomogeneous_vel_bcs(mla,vel_bc_n,vel_bc_t,eta_ed,dx,time, &
                                   the_bc_tower%bc_tower_array)

    ! compute diffusive and stochastic mass fluxes
    ! this computes "-F = rho*W*chi*Gamma*grad(x) - ..."
    call compute_mass_fluxdiv(mla,rho_old,rhotot_old,gradp_baro,Temp, &
                              diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                              diff_mass_flux,stoch_mass_flux, &
                              0.5d0*dt,time,dx,weights,the_bc_tower)

    ! assemble total fluxes to be used in reservoirs
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(total_mass_flux(n,i),1,diff_mass_flux(n,i),1,nspecies,0)
          if (variance_coef_mass .ne. 0.d0) then
             call multifab_plus_plus_c(total_mass_flux(n,i),1,stoch_mass_flux(n,i),1,nspecies,0)
          end if
       end do
    end do

    ! compute chemical rates m_i*R^n_i
    if (nreactions>0) then

       ! convert rho_old (at n) (mass densities rho_i) into n_old (number densities n_i=rho_i/m_i)
       do n=1,nlevs
          call multifab_copy_c(n_old(n),1,rho_old(n),1,nspecies,0)
          do i=1,nspecies
             call multifab_div_div_s_c(n_old(n),i,molmass(i),1,0)
          end do
       end do

       ! compute chemical rates R^n_i (units=[number density]/[time])
       call chemical_rates(mla,n_old,chem_rate,dx,0.5d0*dt)

       ! convert chemical rates R^n_i into m_i*R^n_i (units=[mass density]/[time])
       do n=1,nlevs
          do i=1,nspecies
             call multifab_mult_mult_s_c(chem_rate(n),i,molmass(i),1,0)
          end do
       end do
    end if

    ! set the Dirichlet velocity value on reservoir faces
    call reservoir_bc_fill(mla,total_mass_flux,vel_bc_n,the_bc_tower%bc_tower_array)

    ! compute gmres_rhs_p
    ! put -S = sum_i div(F^n_i)/rhobar_i into gmres_rhs_p (we will later add divu)
    ! if nreactions>0, also add sum_i -(m_i*R^n_i)/rhobar_i
    do n=1,nlevs
       call setval(gmres_rhs_p(n),0.d0,all=.true.)
       do i=1,nspecies
          call saxpy(gmres_rhs_p(n),1,-1.d0/rhobar(i),diff_mass_fluxdiv(n),i,1)
          if (variance_coef_mass .ne. 0.d0) then
             call saxpy(gmres_rhs_p(n),1,-1.d0/rhobar(i),stoch_mass_fluxdiv(n),i,1)
          end if
          if (nreactions > 0) then
             call saxpy(gmres_rhs_p(n),1,-1.d0/rhobar(i),chem_rate(n),i,1)
          end if
       end do
    end do

    ! reset rho_update for all scalars to zero
    ! then, add -F^n_i (plus m_i*R^n_i, if nreactions>0)
    do n=1,nlevs
       ! add fluxes
       call multifab_copy_c(rho_update(n),1,diff_mass_fluxdiv(n),1,nspecies)
       if (variance_coef_mass .ne. 0.d0) then
          call multifab_plus_plus_c(rho_update(n),1,stoch_mass_fluxdiv(n),1,nspecies)
       end if
       if (nreactions > 0) then
          call multifab_plus_plus_c(rho_update(n),1,chem_rate(n),1,nspecies,0)
       end if
    end do

    ! modify umac to respect the boundary conditions we want after the next gmres solve
    ! thus when we add A_0^n v^{n-1/2} to gmres_rhs_v and add div v^{n-1/2} to gmres_rhs_p
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

    ! add A_0^n v^{n-1/2} to gmres_rhs_v
    call diffusive_m_fluxdiv(mla,gmres_rhs_v,.true.,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    ! compute div v^{n-1/2} and add to gmres_rhs_p
    ! now gmres_rhs_p = div v^{n-1/2} - S^n
    ! the sign convention is correct since we solve -div(delta v) = gmres_rhs_p
    call compute_div(mla,umac,gmres_rhs_p,dx,1,1,1,increment_in=.true.)

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
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dpi,rhotot_fc, &
               eta,eta_ed,kappa,theta_alpha,norm_pre_rhs)

    ! for the corrector gmres solve we want the stopping criteria based on the
    ! norm of the preconditioned rhs from the predictor gmres solve.  otherwise
    ! for cases where du in the corrector should be small the gmres stalls
    gmres_abs_tol = max(gmres_abs_tol_in, norm_pre_rhs*gmres_rel_tol)

    ! compute v^* = v^{n-1/2} + delta v
    ! compute pi^* = pi^{n-1/2} + delta pi
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
    ! Step 3 - Scalar Predictor Midpoint Euler Step
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! add A^n for scalars to rho_update
    if (advection_type .ge. 1) then
      do n=1,nlevs
         ! set to zero to make sure ghost cells behind physical boundaries don't have NaNs
         call setval(bds_force(n),0.d0,all=.true.)
         call multifab_copy_c(bds_force(n),1,rho_update(n),1,nspecies,0)
         call multifab_fill_boundary(bds_force(n))
      end do

      if (advection_type .eq. 1 .or. advection_type .eq. 2) then

          ! rho_fc (computed above) and rho_nd_old (computed here) are used to set boundary conditions
          do n=1,nlevs
             call multifab_build_nodal(rho_nd_old(n),mla%la(n),nspecies,1)
          end do
          call average_cc_to_node(nlevs,rho_old,rho_nd_old,1,c_bc_comp,nspecies,the_bc_tower%bc_tower_array)

          ! the input s_tmp needs to have ghost cells filled with multifab_physbc_extrap
          ! instead of multifab_physbc
          do n=1,nlevs
             call multifab_build(rho_tmp(n),mla%la(n),nspecies,rho_old(n)%ng)
             call multifab_copy(rho_tmp(n),rho_old(n),rho_tmp(n)%ng)
             call multifab_physbc_extrap(rho_tmp(n),1,c_bc_comp,nspecies, &
                                         the_bc_tower%bc_tower_array(n))
          end do

          call bds(mla,umac,rho_tmp,rho_update,bds_force,rho_fc,rho_nd_old,dx,0.5d0*dt,1,nspecies, &
                   c_bc_comp,the_bc_tower,proj_type_in=2)

      else if (advection_type .eq. 3 .or. advection_type .eq. 4) then
          call bds_quad(mla,umac,rho_old,rho_update,bds_force,rho_fc,dx,0.5d0*dt,1,nspecies, &
                        c_bc_comp,the_bc_tower,proj_type_in=2)
      end if
    else
       call mk_advective_s_fluxdiv(mla,umac,rho_fc,rho_update,.true.,dx,1,nspecies)
    end if

    ! compute s^{*,n+1/2} = s^n + (dt/2) * (A^n + F^n)
    do n=1,nlevs
       call multifab_saxpy_4(rho_new(n),rho_old(n),0.5d0*dt,rho_update(n))
    end do

    ! compute rhotot from rho in VALID REGION
    call compute_rhotot(mla,rho_new,rhotot_new)

    ! fill rho and rhotot ghost cells
    call fill_rho_rhotot_ghost(mla,rho_new,rhotot_new,dx,the_bc_tower)

    call average_cc_to_face(nlevs,   rho_new,   rho_fc,1,   c_bc_comp,nspecies,the_bc_tower%bc_tower_array)
    call average_cc_to_face(nlevs,rhotot_new,rhotot_fc,1,scal_bc_comp,       1,the_bc_tower%bc_tower_array)

    ! compute (eta,kappa)^{*,n+1/2}
    call compute_eta_kappa(mla,eta,eta_ed,kappa,rho_new,rhotot_new,Temp,dx, &
                           the_bc_tower%bc_tower_array)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 4 - Corrector Stochastic/Diffusive Fluxes
    ! Step 5 - Corrector Stokes Solve
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute grad pi^*
    call compute_grad(mla,pi,gradpi,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

    if (barodiffusion_type .eq. 2) then
       ! barodiffusion uses predicted grad(pi)
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(gradp_baro(n,i),1,gradpi(n,i),1,1,0)
          end do
       end do
    else if (barodiffusion_type .eq. 3) then
       ! compute p0 from rho0*g
       call compute_HSE_pres(mla,rhotot_new,p_baro,dx,the_bc_tower)
       call compute_grad(mla,p_baro,gradp_baro,dx,1,pres_bc_comp,1,1, &
                         the_bc_tower%bc_tower_array)
    end if

    ! set gmres_rhs_v to -grad pi^*
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_5(gmres_rhs_v(n,i),-1.d0,gradpi(n,i),0.d0,gradpi(n,i))
       end do
    end do

    ! add div(Sigma^(2)) to gmres_rhs_v
    if (variance_coef_mom .ne. 0.d0) then
       weights(:) = 1.d0/sqrt(2.d0)
       call stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,gmres_rhs_v,.true., &
                                 eta,eta_ed,Temp,Temp_ed,dx,dt,weights)
    end if

    ! add rho^{*,n+1}*g to gmres_rhs_v
    if (any(grav(1:dm) .ne. 0.d0)) then
       call mk_grav_force(mla,gmres_rhs_v,.true.,rhotot_fc,rhotot_fc,the_bc_tower)
    end if

    ! reset inhomogeneous bc condition to deal with reservoirs
    call set_inhomogeneous_vel_bcs(mla,vel_bc_n,vel_bc_t,eta_ed,dx,time+0.5d0*dt, &
                                   the_bc_tower%bc_tower_array)

    if (midpoint_stoch_mass_flux_type .eq. 1) then
       ! strato

       ! compute diffusive and stochastic mass fluxes
       ! this computes "-F = rho*W*chi*Gamma*grad(x) - ..."
       weights(:) = 1.d0/sqrt(2.d0)
       call compute_mass_fluxdiv(mla,rho_new,rhotot_new,gradp_baro,Temp, &
                                 diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                 diff_mass_flux,stoch_mass_flux, &
                                 dt,time,dx,weights,the_bc_tower)

       ! assemble total fluxes to be used in reservoirs
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(total_mass_flux(n,i),1,diff_mass_flux(n,i),1,nspecies,0)
             if (variance_coef_mass .ne. 0.d0) then
                call multifab_plus_plus_c(total_mass_flux(n,i),1,stoch_mass_flux(n,i),1,nspecies,0)
             end if
          end do
       end do

    else if (midpoint_stoch_mass_flux_type .eq. 2) then
       ! ito

       if (variance_coef_mass .ne. 0.d0) then
          ! for ito interpretation we need to save stoch_mass_fluxdiv_old here
          ! then later add it to stoch_mass_fluxdiv and multiply by 1/2
          do n=1,nlevs
             call multifab_copy_c(stoch_mass_fluxdiv_old(n),1,stoch_mass_fluxdiv(n),1,nspecies,0)
             do i=1,dm
                call multifab_copy_c(stoch_mass_flux_old(n,i),1,stoch_mass_flux(n,i),1,nspecies,0)
             end do
          end do
       end if

       ! compute diffusive and stochastic mass fluxes
       ! this computes "-F = rho*W*chi*Gamma*grad(x) - ..."
       weights(1) = 0.d0
       weights(2) = 1.d0
       call compute_mass_fluxdiv(mla,rho_new,rhotot_new,gradp_baro,Temp, &
                                 diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                 diff_mass_flux,stoch_mass_flux, &
                                 0.5d0*dt,time,dx,weights,the_bc_tower)


       ! assemble total fluxes to be used in reservoirs
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(total_mass_flux(n,i),1,diff_mass_flux(n,i),1,nspecies,0)
             if (variance_coef_mass .ne. 0.d0) then
                call multifab_saxpy_3_cc(total_mass_flux(n,i),1,0.5d0,stoch_mass_flux_old(n,i),1,nspecies)
                call multifab_saxpy_3_cc(total_mass_flux(n,i),1,0.5d0,stoch_mass_flux(n,i),1,nspecies)
             end if
          end do
       end do

       if (variance_coef_mass .ne. 0.d0) then
          ! add stoch_mass_fluxdiv_old to stoch_mass_fluxdiv and multiply by 1/2
          do n=1,nlevs
             call multifab_plus_plus_c(stoch_mass_fluxdiv(n),1,stoch_mass_fluxdiv_old(n),1,nspecies,0)
             call multifab_mult_mult_s_c(stoch_mass_fluxdiv(n),1,0.5d0,nspecies,0)
          end do
       end if

    end if

    ! compute chemical rates m_i*R^{n+1/2}_i
    if (nreactions > 0) then
       ! convert rho_new (at n+1/2) (mass densities rho_i) into n_new (number densities n_i=rho_i/m_i)
       do n=1,nlevs
          call multifab_copy_c(n_new(n),1,rho_new(n),1,nspecies,0)
          do i=1,nspecies
             call multifab_div_div_s_c(n_new(n),i,molmass(i),1,0)
          end do
       end do

       ! compute chemical rates R_i (units=[number density]/[time]) for the second half step
       call chemical_rates(mla,n_old,chem_rate_temp,dx,0.5d0*dt,n_new,mattingly_lin_comb_coef)

       ! convert chemical rates R_i into m_i*R_i (units=[mass density]/[time])
       do n=1,nlevs
          do i=1,nspecies
             call multifab_mult_mult_s_c(chem_rate_temp(n),i,molmass(i),1,0)
          end do
       end do

       ! compute chemical rates m_i*R^{n+1/2}_i for the full step
       do n=1,nlevs
          call multifab_plus_plus_c(chem_rate(n),1,chem_rate_temp(n),1,nspecies,0)
          call multifab_mult_mult_s_c(chem_rate(n),1,0.5d0,nspecies,0)
       end do

       if (use_bl_rng) then
          ! save random state for restart
          call bl_rng_copy_engine(rng_eng_reaction_chk,rng_eng_reaction)
       end if

    end if

    ! set the Dirichlet velocity value on reservoir faces
    call reservoir_bc_fill(mla,total_mass_flux,vel_bc_n,the_bc_tower%bc_tower_array)

    ! compute gmres_rhs_p
    ! put -S = sum_i div(F^{n+1/2}_i)/rhobar_i into gmres_rhs_p (we will later add divu)
    ! if nreactions>0, also add sum_i -(m_i*R^{n+1/2}_i)/rhobar_i
    do n=1,nlevs
       call setval(gmres_rhs_p(n),0.d0,all=.true.)
       do i=1,nspecies
          call saxpy(gmres_rhs_p(n),1,-1.d0/rhobar(i),diff_mass_fluxdiv(n),i,1)
          if (variance_coef_mass .ne. 0.d0) then
             call saxpy(gmres_rhs_p(n),1,-1.d0/rhobar(i),stoch_mass_fluxdiv(n),i,1)
          end if
          if (nreactions > 0) then
             call saxpy(gmres_rhs_p(n),1,-1.d0/rhobar(i),chem_rate(n),i,1)
          end if
       end do
    end do

    ! reset rho_update for all scalars to zero
    ! then, add -F^{n+1/2}_i (plus m_i*R^{n+1/2}_i, if nreactions>0)
    do n=1,nlevs
       call multifab_copy_c(rho_update(n),1,diff_mass_fluxdiv(n),1,nspecies)
       if (variance_coef_mass .ne. 0.d0) then
          call multifab_plus_plus_c(rho_update(n),1,stoch_mass_fluxdiv(n),1,nspecies)
       end if
       if (nreactions > 0) then
          call multifab_plus_plus_c(rho_update(n),1,chem_rate(n),1,nspecies)
       end if
    end do

    ! modify umac to respect the boundary conditions we want after the next gmres solve
    ! thus when we add A_0^* v^* to gmres_rhs_v and add div v^* to gmres_rhs_p
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

    ! add A_0^* v^* to gmres_rhs_v
    call diffusive_m_fluxdiv(mla,gmres_rhs_v,.true.,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    ! compute div v^*and add to gmres_rhs_p
    ! now gmres_rhs_p = div v^* - S^{*,n+1/2}
    ! the sign convention is correct since we solve -div(delta v) = gmres_rhs_p
    call compute_div(mla,umac,gmres_rhs_p,dx,1,1,1,increment_in=.true.)

    ! set the initial guess to zero
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(dumac(n,i),0.d0,all=.true.)
       end do
       call multifab_setval(dpi(n),0.d0,all=.true.)
    end do

    ! call gmres to compute delta v and delta pi
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dpi,rhotot_fc, &
               eta,eta_ed,kappa,theta_alpha)
                              
    gmres_abs_tol = gmres_abs_tol_in ! Restore the desired tolerance         

    ! compute v^{n+1/2} = v^* + dumac
    ! compute pi^{n+1/2} = pi^* + dpi
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
    ! Step 6 - Midpoint Scalar Corrector
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! add A^{n+1/2} for scalars to rho_update
    if (advection_type .ge. 1) then
      do n=1,nlevs
         call multifab_copy_c(bds_force(n),1,rho_update(n),1,2,0)
         call multifab_fill_boundary(bds_force(n))
      end do
      if (advection_type .eq. 1 .or. advection_type .eq. 2) then

          call bds(mla,umac,rho_tmp,rho_update,bds_force,rho_fc,rho_nd_old,dx,dt,1,nspecies, &
                   c_bc_comp,the_bc_tower,proj_type_in=2)

          do n=1,nlevs
             call multifab_destroy(rho_nd_old(n))
             call multifab_destroy(rho_tmp(n))
          end do

      else if (advection_type .eq. 3 .or. advection_type .eq. 4) then
          call bds_quad(mla,umac,rho_old,rho_update,bds_force,rho_fc,dx,dt,1,nspecies, &
                        c_bc_comp,the_bc_tower,proj_type_in=2)
      end if
    else
       call mk_advective_s_fluxdiv(mla,umac,rho_fc,rho_update,.true.,dx,1,nspecies)
    end if

    ! compute s^{n+1} = s^n + dt * (A^{n+1/2} + F^{*,n+1/2})
    do n=1,nlevs
       call multifab_saxpy_4(rho_new(n),rho_old(n),dt,rho_update(n))
    end do

    ! project rho onto eos
    if ( project_eos_int .gt. 0 .and. mod(istep,project_eos_int) .eq. 0) then
       call project_onto_eos(mla,rho_new)
    end if

    ! compute rhotot from rho in VALID REGION
    call compute_rhotot(mla,rho_new,rhotot_new)

    ! fill rho and rhotot ghost cells
    call fill_rho_rhotot_ghost(mla,rho_new,rhotot_new,dx,the_bc_tower)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Compute stuff for plotfile and next time step

    ! compute (eta,kappa)^{n+1}
    call compute_eta_kappa(mla,eta,eta_ed,kappa,rho_new,rhotot_new,Temp,dx, &
                           the_bc_tower%bc_tower_array)

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
          call multifab_destroy(gmres_rhs_v(n,i))
          call multifab_destroy(dumac(n,i))
          call multifab_destroy(gradpi(n,i))
          call multifab_destroy(rho_fc(n,i))
          call multifab_destroy(rhotot_fc(n,i))
          call multifab_destroy(diff_mass_flux(n,i))
          call multifab_destroy(total_mass_flux(n,i))
       end do
    end do

    if (variance_coef_mass .ne. 0.d0) then
       do n=1,nlevs
          call multifab_destroy(stoch_mass_fluxdiv_old(n))
          do i=1,dm
             call multifab_destroy(stoch_mass_flux_old(n,i))
          end do
       end do
    end if

    if (nreactions > 0) then
       do n=1,nlevs
          call multifab_destroy(chem_rate_temp(n))
          call multifab_destroy(n_old(n))
          call multifab_destroy(n_new(n))
       end do
    end if

    call destroy(bpt)

  end subroutine advance_timestep_overdamped

end module advance_timestep_overdamped_module
