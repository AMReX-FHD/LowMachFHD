module advance_timestep_bousq_module

  use ml_layout_module
  use define_bc_module
  use bc_module
  use convert_stag_module
  use mk_advective_s_fluxdiv_module
  use diffusive_m_fluxdiv_module
  use stochastic_m_fluxdiv_module
  use stochastic_mass_fluxdiv_module
  use compute_mass_fluxdiv_module
  use compute_HSE_pres_module
  use convert_m_to_umac_module
  use convert_rhoc_to_c_module
  use mk_advective_m_fluxdiv_module
  use reservoir_bc_fill_module
  use gmres_module
  use div_and_grad_module
  use compute_mixture_properties_module
  use mass_flux_utilities_module
  use multifab_physbc_module
  use multifab_physbc_stag_module
  use zero_edgeval_module
  use fill_rho_ghost_cells_module
  use ml_solve_module
  use bl_rng_module
  use bl_random_module
  use chemical_rates_module
  use fluid_charge_module
  use mk_grav_force_module
  use bds_module
  use reversible_stress_module
  use project_onto_eos_module
  use probin_common_module, only: advection_type, grav, variance_coef_mass, &
                                  variance_coef_mom, barodiffusion_type, project_eos_int, &
                                  molmass, use_bl_rng, nspecies, plot_int, max_step
  use probin_gmres_module, only: gmres_abs_tol, gmres_rel_tol
  use probin_multispecies_module, only: midpoint_stoch_mass_flux_type, is_ideal_mixture, use_multiphase
  use probin_charged_module, only: use_charged_fluid, electroneutral, relxn_param_charge
  use probin_chemistry_module, only: nreactions, use_Poisson_rng, include_discrete_LMA_correction, &
                                     exclude_solvent_comput_rates, use_mole_frac_LMA

  use fabio_module

  implicit none

  private

  public :: advance_timestep_bousq

contains

  subroutine advance_timestep_bousq(mla,umac,rho_old,rho_new,rhotot_old,rhotot_new, &
                                    gradp_baro,pi,eta,eta_ed,kappa,Temp,Temp_ed, &
                                    diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                    stoch_mass_flux,chem_rate, &
                                    dx,dt,time,the_bc_tower,istep, &
                                    grad_Epot_old,grad_Epot_new,charge_old,charge_new, &
                                    Epot,permittivity,total_mass_flux)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: umac(:,:)             ! enters as t^n, leaves as t^{n+1}
    type(multifab) , intent(inout) :: rho_old(:)
    type(multifab) , intent(inout) :: rho_new(:)
    type(multifab) , intent(inout) :: rhotot_old(:)
    type(multifab) , intent(inout) :: rhotot_new(:)
    type(multifab) , intent(inout) :: gradp_baro(:,:)
    type(multifab) , intent(inout) :: pi(:)
    type(multifab) , intent(inout) :: eta(:)                ! enters as t^n, leaves as t^{n+1}
    type(multifab) , intent(inout) :: eta_ed(:,:)           ! enters as t^n, leaves as t^{n+1}
    type(multifab) , intent(inout) :: kappa(:)              ! enters as t^n, leaves as t^{n+1}
    type(multifab) , intent(inout) :: Temp(:)
    type(multifab) , intent(inout) :: Temp_ed(:,:)
    ! We only really pass these non-persistent multifabs in since we have them in main so we don't want to rebuild them
    ! But the values returned are not actually the correct total flux divergencies and these arrays should not be used
    type(multifab) , intent(inout) :: diff_mass_fluxdiv(:)  ! not persistent
    type(multifab) , intent(inout) :: stoch_mass_fluxdiv(:) ! not persistent
    type(multifab) , intent(inout) :: stoch_mass_flux(:,:)  ! not persistent
    type(multifab) , intent(inout) :: chem_rate(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt,time
    type(bc_tower) , intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: istep
    type(multifab) , intent(inout) :: grad_Epot_old(:,:)
    type(multifab) , intent(inout) :: grad_Epot_new(:,:)
    type(multifab) , intent(inout) :: charge_old(:)
    type(multifab) , intent(inout) :: charge_new(:)
    ! These should not be relied upon to have consistent values
    type(multifab) , intent(inout) :: Epot(:)               ! not persistent
    type(multifab) , intent(inout) :: permittivity(:)       ! not persistent
    ! If desired we can return the total fluxes for all species for purposes of computing total currents etc.
    ! Total fluxes consistent with update from time n to n+1
    type(multifab) , intent(inout), optional :: total_mass_flux(:,:)

    ! local
    type(multifab) :: adv_mass_fluxdiv(mla%nlevel)
    type(multifab) ::      gmres_rhs_p(mla%nlevel)
    type(multifab) ::              dpi(mla%nlevel)

    type(multifab) ::             umac_old(mla%nlevel,mla%dim)
    type(multifab) ::                mtemp(mla%nlevel,mla%dim)
    type(multifab) ::  adv_mom_fluxdiv_old(mla%nlevel,mla%dim)
    type(multifab) ::  adv_mom_fluxdiv_new(mla%nlevel,mla%dim)
    type(multifab) :: diff_mom_fluxdiv_old(mla%nlevel,mla%dim)
    type(multifab) :: diff_mom_fluxdiv_new(mla%nlevel,mla%dim)
    type(multifab) ::    stoch_mom_fluxdiv(mla%nlevel,mla%dim)
    type(multifab) ::          gmres_rhs_v(mla%nlevel,mla%dim)
    type(multifab) ::                dumac(mla%nlevel,mla%dim)
    type(multifab) ::               gradpi(mla%nlevel,mla%dim)
    type(multifab) ::               rho_fc(mla%nlevel,mla%dim)
    type(multifab) ::        rhotot_fc_old(mla%nlevel,mla%dim)
    type(multifab) ::        rhotot_fc_new(mla%nlevel,mla%dim)
    type(multifab) ::       diff_mass_flux(mla%nlevel,mla%dim)
    type(multifab) ::       mom_grav_force(mla%nlevel,mla%dim)    

    ! only used when variance_coef_mass>0 and midpoint_stoch_mass_flux_type=2
    type(multifab) :: stoch_mass_fluxdiv_old(mla%nlevel)
    type(multifab) ::    stoch_mass_flux_old(mla%nlevel,mla%dim)

    ! only used when nreactions>0
    type(multifab) :: chem_rate_temp(mla%nlevel)
    type(multifab) ::          n_old(mla%nlevel)
    type(multifab) ::          n_new(mla%nlevel)

    ! only used when use_charged_fluid=T
    type(multifab) :: Lorentz_force(mla%nlevel,mla%dim)

    type(multifab) :: div_reversible_stress(mla%nlevel,mla%dim)

    ! only needed if total_mass_flux is passed in
    ! this a temporary used to hold the advective mass fluxes
    type(multifab) :: adv_mass_flux(mla%nlevel,mla%dim)

    ! only used for bds advection
    type(multifab) ::    rho_tmp(mla%nlevel)
    type(multifab) :: rho_update(mla%nlevel)
    type(multifab) ::  bds_force(mla%nlevel)
    type(multifab) ::     rho_nd(mla%nlevel)
    type(multifab) ::   umac_tmp(mla%nlevel,mla%dim)

    integer :: i,dm,n,nlevs,proj_type,comp
    
    real(kind=dp_t) :: theta_alpha, norm_pre_rhs, gmres_abs_tol_in, relxn_param_charge_in
    real(kind=dp_t) :: weights(2)

    real(kind=dp_t), parameter :: mattingly_lin_comb_coef(1:2) = (/-1.d0, 2.d0/)

    if (barodiffusion_type .ne. 0) then
       call bl_error("advance_timestep_bousq: barodiffusion not supported yet")
    end if

    if (nreactions > 0) then
       if (use_Poisson_rng .eq. 2) then
          call bl_error('Error: currently use_Poisson_rng=2 not allowed for algorith_type=6 and nreactions>0')
       end if

       if (use_mole_frac_LMA) then
          if (.not. is_ideal_mixture) then
             call bl_error('Error: currently use_mole_frac_LMA can be used only with is_ideal_mixture = T')
          end if
          if (exclude_solvent_comput_rates .ne. 0) then
             call bl_error('Error: currently use_mole_frac_LMA can be used only with exclude_solvent_comput_rates=0')
          end if
       end if
    end if
    
    ! For BDS we need to project edge states onto constraints
    if(use_charged_fluid .and. electroneutral) then
      proj_type=4
    else
      proj_type=3
    end if    

    nlevs = mla%nlevel
    dm = mla%dim

    theta_alpha = 1.d0/dt
    
    do n=1,nlevs
       call multifab_build(adv_mass_fluxdiv(n),mla%la(n),nspecies,0)
       call multifab_build(     gmres_rhs_p(n),mla%la(n),1       ,0)
       call multifab_build(             dpi(n),mla%la(n),1       ,1)
       do i=1,dm
          call multifab_build_edge(             umac_old(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(                mtemp(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(  adv_mom_fluxdiv_old(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(  adv_mom_fluxdiv_new(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge( diff_mom_fluxdiv_old(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge( diff_mom_fluxdiv_new(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(    stoch_mom_fluxdiv(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(          gmres_rhs_v(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(                dumac(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(               gradpi(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(               rho_fc(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge(        rhotot_fc_old(n,i),mla%la(n),nspecies,1,i)
          call multifab_build_edge(        rhotot_fc_new(n,i),mla%la(n),nspecies,1,i)
          call multifab_build_edge(       diff_mass_flux(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge(       mom_grav_force(n,i),mla%la(n),1       ,0,i)
       end do
    end do

    ! for ito interpretation we need to save stoch_mass_fluxdiv_old 
    if (variance_coef_mass .ne. 0.d0 .and. midpoint_stoch_mass_flux_type .eq. 2) then
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

    if (use_charged_fluid) then
       do n=1,nlevs
          do i=1,dm
             call multifab_build_edge(Lorentz_force(n,i),mla%la(n),1,0,i)
          end do
       end do
    end if

    if (use_multiphase) then
       do n=1,nlevs
          do i=1,dm
             call multifab_build_edge(div_reversible_stress(n,i),mla%la(n),1,0,i)
          end do
       end do
    end if

    if (present(total_mass_flux)) then
       do n=1,nlevs
          do i=1,dm
             call multifab_build_edge(adv_mass_flux(n,i),mla%la(n),nspecies,0,i)
          enddo
       enddo
    end if
    
    ! make a copy of umac at t^n
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(umac_old(n,i),1,umac(n,i),1,1,1)
       end do
    end do    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 1: solve for v^{n+1,*} and pi^{n+1/2,*} using GMRES
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! average rho_i^n and rho^n to faces
    call average_cc_to_face(nlevs,   rho_old,   rho_fc    ,1,   c_bc_comp,nspecies,the_bc_tower%bc_tower_array)
    call average_cc_to_face(nlevs,rhotot_old,rhotot_fc_old,1,scal_bc_comp,       1,the_bc_tower%bc_tower_array)

    ! compute mtemp = (rho*v)^n
    call convert_m_to_umac(mla,rhotot_fc_old,mtemp,umac,.false.)

    ! set gmres_rhs_v = mtemp / dt
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(gmres_rhs_v(n,i),1,mtemp(n,i),1,1,0)
          call multifab_mult_mult_s_c(gmres_rhs_v(n,i),1,1.d0/dt,1,0)
       end do
    end do

    ! compute adv_mom_fluxdiv_old = -rho*v^n*v^n
    ! save this for use in the corrector GMRES solve
    call mk_advective_m_fluxdiv(mla,umac,mtemp,adv_mom_fluxdiv_old,.false., &
                                dx,the_bc_tower%bc_tower_array)

    ! add -rho*v^n,v^n to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,1.d0,adv_mom_fluxdiv_old(n,i),1,1)
       end do
    end do

    ! compute diff_mom_fluxdiv_old = L_0^n v^n
    ! save this for use in the corrector GMRES solve
    call diffusive_m_fluxdiv(mla,diff_mom_fluxdiv_old,.false.,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    ! add (1/2)*diff_mom_fluxdiv_old to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,diff_mom_fluxdiv_old(n,i),1,1)
       end do
    end do

    if (variance_coef_mom .ne. 0.d0) then

       ! fill the stochastic momentum multifabs with new sets of random numbers
       call fill_m_stochastic(mla)

       ! compute stoch_mom_fluxdiv = div (sqrt(eta^n...) Wbar^n)
       ! save this for use in the corrector GMRES solve
       weights(1) = 1.d0
       call stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,stoch_mom_fluxdiv,.false., &
                                 eta,eta_ed,Temp,Temp_ed,dx,dt,weights)

       ! add stochastic momentum fluxes to gmres_rhs_v
       do n=1,nlevs
          do i=1,dm
             call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,1.d0,stoch_mom_fluxdiv(n,i),1,1)
          end do
       end do

    end if

    ! gravity
    if (any(grav(1:dm) .ne. 0.d0)) then
       call mk_grav_force_bousq(mla,gmres_rhs_v,.true.,rho_fc,the_bc_tower)
    end if

    if (variance_coef_mass .ne. 0.d0) then
       call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)
       if (use_bl_rng) then
          ! save random state for restart
          call bl_rng_copy_engine(rng_eng_diffusion_chk,rng_eng_diffusion)
       end if
    end if

    ! compute diffusive, stochastic, potential mass fluxes
    ! with barodiffusion and thermodiffusion
    ! this computes "-F = rho W chi [Gamma grad x... ]"
    weights(1) = 1.d0
    weights(2) = 0.d0
    ! For electroneutral, we only have to do charge relaxation in the corrector
    relxn_param_charge_in=relxn_param_charge
    relxn_param_charge=0.0 ! Don't correct in predictor
    call compute_mass_fluxdiv(mla,rho_old,rhotot_old,gradp_baro,Temp, &
                              diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                              diff_mass_flux,stoch_mass_flux, &
                              0.5d0*dt,time,dx,weights,the_bc_tower, &
                              charge_old,grad_Epot_old,Epot, &
                              permittivity)

    ! here is a reasonable place to call something to compute in reversible stress term
    ! in this case want to get divergence so it looks like a add to rhs for stokes solver
    if(use_multiphase) then

     !compute reversible stress tensor ---added term
     call compute_div_reversible_stress(mla,div_reversible_stress,rhotot_old,rho_old,dx,the_bc_tower)

      ! add divergence of reversible stress to gmres_rhs_v
      do n=1,nlevs
         do i=1,dm
            call multifab_saxpy_3(gmres_rhs_v(n,i),1.d0,div_reversible_stress(n,i))
         end do
      end do

    end if

    if (use_charged_fluid) then

       ! compute old Lorentz force
       call compute_Lorentz_force(mla,Lorentz_force,grad_Epot_old,permittivity, &
                                  charge_old,dx,the_bc_tower)

       ! add Lorentz force to gmres_rhs_v
       do n=1,nlevs
          do i=1,dm
             call multifab_saxpy_3(gmres_rhs_v(n,i),1.d0,Lorentz_force(n,i))
          end do
       end do

    end if

    ! compute grad pi^{n-1/2}
    call compute_grad(mla,pi,gradpi,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

    ! subtract grad pi^{n-1/2} from gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradpi(n,i),1,1,0)
       end do
    end do

    ! set inhomogeneous velocity bc's to values supplied in inhomogeneous_bc_val
    call set_inhomogeneous_vel_bcs(mla,vel_bc_n,vel_bc_t,eta_ed,dx,time+dt, &
                                   the_bc_tower%bc_tower_array)

    ! modify umac to respect the boundary conditions we want after the next gmres solve
    ! thus when we add L_0^n vbar^n to gmres_rhs_v and add div vbar^n to gmres_rhs_p we
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

    ! compute mtemp = (rho*vbar)^n
    call convert_m_to_umac(mla,rhotot_fc_old,mtemp,umac,.false.)

    ! add -(mtemp/dt) to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,-1.d0/dt,mtemp(n,i),1,1)
       end do
    end do

    ! compute diff_mom_fluxdiv_new = L_0^n vbar^n
    call diffusive_m_fluxdiv(mla,diff_mom_fluxdiv_new,.false.,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    ! add (1/2)*diff_mom_fluxdiv_new to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,diff_mom_fluxdiv_new(n,i),1,1)
       end do
    end do

    ! compute div(vbar^n) and store in gmres_rhs_p
    ! the sign convention is correct since we solve -div(delta v) = div(vbar^n)
    call compute_div(mla,umac,gmres_rhs_p,dx,1,1,1)

    ! multiply eta and kappa by 1/2 to put in proper form for gmres solve
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,0.5d0,1,eta(n)%ng)
       call multifab_mult_mult_s_c(kappa(n),1,0.5d0,1,kappa(n)%ng)
       do i=1,size(eta_ed,dim=2)
          call multifab_mult_mult_s_c(eta_ed(n,i),1,0.5d0,1,eta_ed(n,i)%ng)
       end do
    end do

    ! set the initial guess to zero
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(dumac(n,i),0.d0,all=.true.)
       end do
       call multifab_setval(dpi(n),0.d0,all=.true.)
    end do

    ! zero gmres_rhs_v on physical boundaries
    do n=1,nlevs
       call zero_edgeval_physical(gmres_rhs_v(n,:),1,1,the_bc_tower%bc_tower_array(n))
    end do

    gmres_abs_tol_in = gmres_abs_tol ! Save this 

    ! This relies entirely on relative tolerance and can fail if the rhs is roundoff error only:
    ! gmres_abs_tol = 0.d0 ! It is better to set gmres_abs_tol in namelist to a sensible value

    ! call gmres to compute delta v and delta pi
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dpi,rhotot_fc_old, &
               eta,eta_ed,kappa,theta_alpha,norm_pre_rhs)    

    ! for the corrector gmres solve we want the stopping criteria based on the
    ! norm of the preconditioned rhs from the predictor gmres solve.  otherwise
    ! for cases where du in the corrector should be small the gmres stalls
    gmres_abs_tol = max(gmres_abs_tol_in, norm_pre_rhs*gmres_rel_tol)
       
    ! compute v^{n+1,*} = vbar^n + dumac
    ! compute pi^{n+1/2,*} = pi^{n-1/2} + dpi
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

    ! restore eta and kappa
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,2.d0,1,eta(n)%ng)
       call multifab_mult_mult_s_c(kappa(n),1,2.d0,1,kappa(n)%ng)
       do i=1,size(eta_ed,dim=2)
          call multifab_mult_mult_s_c(eta_ed(n,i),1,2.d0,1,eta_ed(n,i)%ng)
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 2: compute reactions and mass fluxes at t^n
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute chemical rates m_i*R^n_i
    if (nreactions > 0) then
       
       ! convert rho_old (at t^n) (mass densities rho_i) into n_old (number densities n_i=rho_i/m_i)
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

    if (advection_type .ge. 1) then

       ! bds advection
       do n=1,nlevs
          call multifab_build(     rho_tmp(n),mla%la(n),nspecies,rho_old(n)%ng)
          call multifab_build(  rho_update(n),mla%la(n),nspecies,0)
          call multifab_build(   bds_force(n),mla%la(n),nspecies,1)
          call multifab_build_nodal(rho_nd(n),mla%la(n),nspecies,1)
          do i=1,dm
             call multifab_build_edge(umac_tmp(n,i),mla%la(n),1,1,i)
          end do
       end do

       do n=1,nlevs
          do i=1,dm
             ! create average of umac^n and umac^{n+1,*}
             call multifab_copy_c(umac_tmp(n,i),1,umac_old(n,i),1,1,1)
             call multifab_saxpy_3(umac_tmp(n,i),1.d0,umac(n,i),all=.true.)
             call multifab_mult_mult_s(umac_tmp(n,i),0.5d0,1)
          end do
       end do

       ! add the diff/stoch/react terms to rho_update
       do n=1,nlevs
          call multifab_copy_c(rho_update(n),1,diff_mass_fluxdiv(n),1,nspecies)
          if (variance_coef_mass .ne. 0.d0) then
             call multifab_plus_plus_c(rho_update(n),1,stoch_mass_fluxdiv(n),1,nspecies)
          end if
          if (nreactions > 0) then
             call multifab_plus_plus_c(rho_update(n),1,chem_rate(n),1,nspecies,0)
          end if
       end do

       do n=1,nlevs
          ! set to zero to make sure ghost cells behind physical boundaries don't have NaNs
          call setval(bds_force(n),0.d0,all=.true.)
          call multifab_copy_c(bds_force(n),1,rho_update(n),1,nspecies,0)
          call multifab_fill_boundary(bds_force(n))
       end do

       if (advection_type .eq. 1 .or. advection_type .eq. 2) then

          ! the input rho_tmp needs to have ghost cells filled with multifab_physbc_extrap
          ! instead of multifab_physbc
          do n=1,nlevs
             call multifab_copy(rho_tmp(n),rho_old(n),rho_tmp(n)%ng)
             call multifab_physbc_extrap(rho_tmp(n),1,c_bc_comp,nspecies, &
                                         the_bc_tower%bc_tower_array(n))
          end do

          call average_cc_to_node(nlevs,rho_old,rho_nd,1,c_bc_comp,nspecies,the_bc_tower%bc_tower_array)

          ! bds increments rho_update with the advection term
          call bds(mla,umac_tmp,rho_tmp,rho_update,bds_force,rho_fc,rho_nd,dx,0.5d0*dt,1, &
                   nspecies,c_bc_comp,the_bc_tower,proj_type_in=proj_type)

       else if (advection_type .eq. 3 .or. advection_type .eq. 4) then
          call bds_quad(mla,umac_tmp,rho_old,rho_update,bds_force,rho_fc,dx,0.5d0*dt,1,nspecies, &
                        c_bc_comp,the_bc_tower,proj_type_in=proj_type)
       end if

    else

       ! compute adv_mass_fluxdiv = -rho_i^n * v^n and then
       ! increment adv_mass_fluxdiv by -rho_i^n * v^{n+1,*}
       call mk_advective_s_fluxdiv(mla,umac_old,rho_fc,adv_mass_fluxdiv,.false.,dx,1,nspecies)
       call mk_advective_s_fluxdiv(mla,umac    ,rho_fc,adv_mass_fluxdiv,.true. ,dx,1,nspecies)

    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 3: density prediction to t^{n+1/2}
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (advection_type .ge. 1) then

       do n=1,nlevs
          call multifab_saxpy_4(rho_new(n),rho_old(n),0.5d0*dt,rho_update(n))
       end do

    else

       ! compute rho_i^{n+1/2} (store in rho_new)
       ! multiply adv_mass_fluxdiv by (1/4) since it contains -rho_i^n * (v^n + v^{n+1,*})
       do n=1,nlevs
          call multifab_copy_c(rho_new(n),1,rho_old(n),1,nspecies,0)
          call multifab_saxpy_3_cc(rho_new(n),1,0.25d0*dt, adv_mass_fluxdiv(n),1,nspecies)
          call multifab_saxpy_3_cc(rho_new(n),1, 0.5d0*dt,diff_mass_fluxdiv(n),1,nspecies)
          if (variance_coef_mass .ne. 0.d0) then
             call multifab_saxpy_3_cc(rho_new(n),1,0.5d0*dt,stoch_mass_fluxdiv(n),1,nspecies)
          end if
          if (nreactions > 0) then
             call multifab_saxpy_3_cc(rho_new(n),1,0.5d0*dt,chem_rate(n),1,nspecies)
          end if
       end do

    end if

    ! compute rhotot^{n+1/2} from rho^{n+1/2} in VALID REGION
    call compute_rhotot(mla,rho_new,rhotot_new)

    ! fill rho and rhotot ghost cells at t^{n+1/2}
    call fill_rho_rhotot_ghost(mla,rho_new,rhotot_new,dx,the_bc_tower)

    ! average rho_i^{n+1/2} to faces
    call average_cc_to_face(nlevs,rho_new,rho_fc,1,c_bc_comp,nspecies,the_bc_tower%bc_tower_array)

    if (use_charged_fluid) then
       ! compute total charge at t^{n+1/2}
       call dot_with_z(mla,rho_new,charge_new)
       ! compute permittivity at t^{n+1/2}
       call compute_permittivity(mla,permittivity,rho_new,rhotot_new,the_bc_tower)
    end if

    ! compute (eta,kappa)^{n+1/2}
    call compute_eta_kappa(mla,eta,eta_ed,kappa,rho_new,rhotot_new,Temp,dx, &
                           the_bc_tower%bc_tower_array)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 4: compute mass fluxes and reactions at t^{n+1/2}
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute mass fluxes and reactions at t^{n+1/2}
    ! For electroneutral, enable charge correction in the corrector
    relxn_param_charge=relxn_param_charge_in ! Default value is 1  
    if (midpoint_stoch_mass_flux_type .eq. 1) then
       ! strato

       ! compute diffusive, stochastic, potential mass fluxes
       ! with barodiffusion and thermodiffusion
       ! this computes "-F = rho W chi [Gamma grad x... ]"
       weights(:) = 1.d0/sqrt(2.d0)
       call compute_mass_fluxdiv(mla,rho_new,rhotot_new,gradp_baro,Temp, &
                                 diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                 diff_mass_flux,stoch_mass_flux, &
                                 dt,time,dx,weights,the_bc_tower, &
                                 charge_new,grad_Epot_new,Epot, &
                                 permittivity,zero_initial_Epot_in=.false.)

       ! begin to assemble total fluxes (diffusive plus stochastic);
       ! we still need to add add advective flux later
       if(present(total_mass_flux)) then
          do n=1,nlevs
             do i=1,dm
                call multifab_copy_c(total_mass_flux(n,i),1,diff_mass_flux(n,i),1,nspecies,0)
                if (variance_coef_mass .ne. 0.d0) then
                   call multifab_plus_plus_c(total_mass_flux(n,i),1,stoch_mass_flux(n,i),1,nspecies,0)
                end if
             end do
          end do
       end if

    else if (midpoint_stoch_mass_flux_type .eq. 2) then
       ! ito

       if(present(total_mass_flux)) then
          call bl_error('Code cannot presently return total_mass_flux if midpoint_stoch_mass_flux_type!=1')
       end if

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

       ! compute diffusive, stochastic, potential mass fluxes
       ! with barodiffusion and thermodiffusion
       ! this computes "-F = rho W chi [Gamma grad x... ]"
       weights(1) = 0.d0
       weights(2) = 1.d0
       call compute_mass_fluxdiv(mla,rho_new,rhotot_new,gradp_baro,Temp, &
                                 diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                 diff_mass_flux,stoch_mass_flux, &
                                 0.5d0*dt,time,dx,weights,the_bc_tower, &
                                 charge_new,grad_Epot_new,Epot, &
                                 permittivity,zero_initial_Epot_in=.false.)

       if (variance_coef_mass .ne. 0.d0) then
          ! add stoch_mass_fluxdiv_old to stoch_mass_fluxdiv and multiply by 1/2
          do n=1,nlevs
             call multifab_plus_plus_c(stoch_mass_fluxdiv(n),1,stoch_mass_fluxdiv_old(n),1,nspecies,0)
             call multifab_mult_mult_s_c(stoch_mass_fluxdiv(n),1,0.5d0,nspecies,0)
          end do
       end if

    end if
    
    if(use_multiphase) then

       !compute reversible stress tensor ---added term (will add to gmres_rhs_v later)
       call compute_div_reversible_stress(mla,div_reversible_stress,rhotot_new,rho_new,dx,the_bc_tower)

    end if

    if (use_charged_fluid) then

       ! compute Lorentz force (using midpoint value not trapezoidal, will add to gmres_rhs_v later)
       call compute_Lorentz_force(mla,Lorentz_force,grad_Epot_new,permittivity, &
                                  charge_new,dx,the_bc_tower)

    end if

    ! compute chemical rates m_i*R^{n+1/2}_i
    if (nreactions > 0) then
       ! convert rho_old (at n) and rho_new (at n+1/2) (mass densities rho_i)
       ! into n_old and n_new (number densities n_i=rho_i/m_i)
       do n=1,nlevs
          call multifab_copy_c(n_old(n),1,rho_old(n),1,nspecies,0)
          call multifab_copy_c(n_new(n),1,rho_new(n),1,nspecies,0)
          do i=1,nspecies
             call multifab_div_div_s_c(n_old(n),i,molmass(i),1,0)
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

    if (advection_type .ge. 1) then

       if(present(total_mass_flux)) then
          call bl_error('Code cannot presently return total_mass_flux if using BDS advection')
       end if

       ! add the diff/stoch/react terms to rho_update
       do n=1,nlevs
          call multifab_copy_c(rho_update(n),1,diff_mass_fluxdiv(n),1,nspecies)
          if (variance_coef_mass .ne. 0.d0) then
             call multifab_plus_plus_c(rho_update(n),1,stoch_mass_fluxdiv(n),1,nspecies)
          end if
          if (nreactions > 0) then
             call multifab_plus_plus_c(rho_update(n),1,chem_rate(n),1,nspecies,0)
          end if
       end do

       do n=1,nlevs
          ! set to zero to make sure ghost cells behind physical boundaries don't have NaNs
          call setval(bds_force(n),0.d0,all=.true.)
          call multifab_copy_c(bds_force(n),1,rho_update(n),1,nspecies,0)
          call multifab_fill_boundary(bds_force(n))
       end do

       if (advection_type .eq. 1 .or. advection_type .eq. 2) then
          ! bds increments rho_update with the advection term
          call bds(mla,umac_tmp,rho_tmp,rho_update,bds_force,rho_fc,rho_nd,dx,dt,1, &
                   nspecies,c_bc_comp,the_bc_tower,proj_type_in=proj_type)

       else if (advection_type .eq. 3 .or. advection_type .eq. 4) then
          call bds_quad(mla,umac_tmp,rho_old,rho_update,bds_force,rho_fc,dx,dt,1,nspecies, &
                        c_bc_comp,the_bc_tower,proj_type_in=proj_type)
       end if

       do n=1,nlevs
          call multifab_destroy(rho_tmp(n))
          call multifab_destroy(bds_force(n))
          call multifab_destroy(rho_nd(n))
          do i=1,dm
             call multifab_destroy(umac_tmp(n,i))
          end do
       end do

    else

       ! compute adv_mass_fluxdiv = -rho_i^{n+1/2} * v^n and
       ! increment adv_mass_fluxdiv by -rho_i^{n+1/2} * v^{n+1,*}
       call mk_advective_s_fluxdiv(mla,umac_old,rho_fc,adv_mass_fluxdiv,.false.,dx,1,nspecies)
       call mk_advective_s_fluxdiv(mla,umac    ,rho_fc,adv_mass_fluxdiv,.true. ,dx,1,nspecies)

       ! Add advective fluxes to total fluxes and fix sign for fluxes to be actual fluxes not negative of flux ;-)
       if(present(total_mass_flux)) then

          do n=1,nlevs
             do i=1,dm

                ! compute advective mass flux contribution, rho_i^{n+1/2} * v^n
                do comp=1,nspecies
                   call multifab_copy_c(adv_mass_flux(n,i),comp,umac_old(n,i),1,1,0)
                end do
                call multifab_mult_mult_c(adv_mass_flux(n,i),1,rho_fc(n,i),1,nspecies,0)

                ! now *subtract* off half of these advective fluxes from total_mass_flux
                call multifab_saxpy_3_cc(total_mass_flux(n,i),1,-0.5d0,adv_mass_flux(n,i),1,nspecies)

                ! compute advective mass flux contribution, rho_i^{n+1/2} * v^{n+1,*}
                do comp=1,nspecies
                   call multifab_copy_c(adv_mass_flux(n,i),comp,umac(n,i),1,1,0)
                end do
                call multifab_mult_mult_c(adv_mass_flux(n,i),1,rho_fc(n,i),1,nspecies,0)

                ! now *subtract* off half of these advective fluxes from total_mass_flux
                call multifab_saxpy_3_cc(total_mass_flux(n,i),1,-0.5d0,adv_mass_flux(n,i),1,nspecies)
                
                ! multiply the total mass flux w adv by -1, so that we take the 
                ! standard convention that the mass flux has the same sign as the fluid flow 
                call multifab_mult_mult_s_c(total_mass_flux(n,i),1,-1.d0,nspecies,0)

             enddo
          end do

       do n=1,nlevs
          do i=1,dm
                
          enddo
       enddo        

       end if

    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 5: density integration to t^{n+1}
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (advection_type .ge. 1) then

       do n=1,nlevs
          call multifab_saxpy_4(rho_new(n),rho_old(n),dt,rho_update(n))
       end do

       do n=1,nlevs
          call multifab_destroy(rho_update(n))
       end do

    else

       ! compute rho_i^{n+1}
       ! multiply adv_mass_fluxdiv by (1/2) since it contains -rho_i^{n+1/2} * (v^n + v^{n+1,*})
       do n=1,nlevs
          call multifab_copy_c(rho_new(n),1,rho_old(n),1,nspecies,0)
          call multifab_saxpy_3_cc(rho_new(n),1,0.5d0*dt, adv_mass_fluxdiv(n),1,nspecies)
          call multifab_saxpy_3_cc(rho_new(n),1,      dt,diff_mass_fluxdiv(n),1,nspecies)
          if (variance_coef_mass .ne. 0.d0) then
             call multifab_saxpy_3_cc(rho_new(n),1,dt,stoch_mass_fluxdiv(n),1,nspecies)
          end if
          if (nreactions > 0) then
             call multifab_saxpy_3_cc(rho_new(n),1,dt,chem_rate(n),1,nspecies)
          end if
       end do

    end if

    if ( project_eos_int .gt. 0 .and. mod(istep,project_eos_int) .eq. 0) then
       call project_onto_eos(mla,rho_new)
    end if

    ! compute rhotot^{n+1} from rho^{n+1} in VALID REGION
    call compute_rhotot(mla,rho_new,rhotot_new)

    ! fill rho and rhotot ghost cells at t^{n+1}
    call fill_rho_rhotot_ghost(mla,rho_new,rhotot_new,dx,the_bc_tower)

    ! average rho^{n+1} to faces
    call average_cc_to_face(nlevs,rhotot_new,rhotot_fc_new,1,scal_bc_comp,1,the_bc_tower%bc_tower_array)

    if (use_charged_fluid) then
       ! compute total charge at t^{n+1}
       call dot_with_z(mla,rho_new,charge_new)
       ! compute permittivity at t^{n+1}
       call compute_permittivity(mla,permittivity,rho_new,rhotot_new,the_bc_tower)
    end if

    ! compute (eta,kappa)^{n+1}
    call compute_eta_kappa(mla,eta,eta_ed,kappa,rho_new,rhotot_new,Temp,dx, &
                           the_bc_tower%bc_tower_array)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 6: solve for v^{n+1} and pi^{n+1/2} using GMRES
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute mtemp = (rho*v)^n
    call convert_m_to_umac(mla,rhotot_fc_old,mtemp,umac_old,.false.)

    ! set gmres_rhs_v = mtemp / dt
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(gmres_rhs_v(n,i),1,mtemp(n,i),1,1,0)
          call multifab_mult_mult_s_c(gmres_rhs_v(n,i),1,1.d0/dt,1,0)
       end do
    end do

    ! compute mtemp = rho^{n+1}*v^{n+1,*} for adv_mom_fluxdiv computation
    call convert_m_to_umac(mla,rhotot_fc_new,mtemp,umac,.false.)

    ! compute adv_mom_fluxdiv_new = -rho^{n+1}*v^{n+1,*}*v^{n+1,*}
    call mk_advective_m_fluxdiv(mla,umac,mtemp,adv_mom_fluxdiv_new,.false., &
                                dx,the_bc_tower%bc_tower_array)

    ! add (1/2) (-rho^n*v^n*v^n -rho^{n+1}*v^{n+1,*}*v^{n+1,*}) to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,adv_mom_fluxdiv_old(n,i),1,1)
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,adv_mom_fluxdiv_new(n,i),1,1)
       end do
    end do

    ! add (1/2)*diff_mom_fluxdiv_old to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,diff_mom_fluxdiv_old(n,i),1,1)
       end do
    end do

    if (variance_coef_mom .ne. 0.d0) then

       ! increment stoch_mom_fluxdiv by = div (sqrt(eta^{n+1}...) Wbar^n)
       weights(1) = 1.d0
       call stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,stoch_mom_fluxdiv,.true., &
                                 eta,eta_ed,Temp,Temp_ed,dx,dt,weights)

       ! add stochastic momentum fluxes to gmres_rhs_v
       ! note the factor of (1/2) since stoch_mom_fluxdiv contains
       ! div (sqrt(eta^{n+1}...) Wbar^n) + div (sqrt(eta^n...) Wbar^n)
       do n=1,nlevs
          do i=1,dm
             call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,stoch_mom_fluxdiv(n,i),1,1)
          end do
       end do

    end if

    ! gravity
    if (any(grav(1:dm) .ne. 0.d0)) then
       call mk_grav_force_bousq(mla,gmres_rhs_v,.true.,rho_fc,the_bc_tower)
    end if

    if(use_multiphase) then

      ! add divergence of reversible stress to gmres_rhs_v
      do n=1,nlevs
         do i=1,dm
            call multifab_saxpy_3(gmres_rhs_v(n,i),1.d0,div_reversible_stress(n,i))
         end do
      end do

    end if
   

    if (use_charged_fluid) then

       ! add Lorentz force to gmres_rhs_v
       do n=1,nlevs
          do i=1,dm
             call multifab_saxpy_3(gmres_rhs_v(n,i),1.d0,Lorentz_force(n,i))
          end do
       end do

    end if

    ! compute grad pi^{n+1/2,*}
    call compute_grad(mla,pi,gradpi,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

    ! subtract grad pi^{n+1/2,*} from gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradpi(n,i),1,1,0)
       end do
    end do

    ! modify umac to respect the boundary conditions we want after the next gmres solve
    ! thus when we add L_0^{n+1} vbar^{n+1,*} to gmres_rhs_v and add div vbar^{n+1,*} to gmres_rhs_p we
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

    ! compute mtemp = rho^{n+1}*vbar^{n+1,*}
    call convert_m_to_umac(mla,rhotot_fc_new,mtemp,umac,.false.)

    ! add -(mtemp/dt) to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,-1.d0/dt,mtemp(n,i),1,1)
       end do
    end do

    ! compute diff_mom_fluxdiv_new = L_0^{n+1} vbar^{n+1,*}
    call diffusive_m_fluxdiv(mla,diff_mom_fluxdiv_new,.false.,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    ! add (1/2)*diff_mom_fluxdiv_new to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,diff_mom_fluxdiv_new(n,i),1,1)
       end do
    end do

    ! compute div(vbar^{n+1,*}) and store in gmres_rhs_p
    ! the sign convention is correct since we solve -div(delta v) = div(vbar^{n+1,*})
    call compute_div(mla,umac,gmres_rhs_p,dx,1,1,1)

    ! multiply eta and kappa by 1/2 to put in proper form for gmres solve
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,0.5d0,1,eta(n)%ng)
       call multifab_mult_mult_s_c(kappa(n),1,0.5d0,1,kappa(n)%ng)
       do i=1,size(eta_ed,dim=2)
          call multifab_mult_mult_s_c(eta_ed(n,i),1,0.5d0,1,eta_ed(n,i)%ng)
       end do
    end do

    ! set the initial guess to zero
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(dumac(n,i),0.d0,all=.true.)
       end do
       call multifab_setval(dpi(n),0.d0,all=.true.)
    end do

    ! zero gmres_rhs_v on physical boundaries
    do n=1,nlevs
       call zero_edgeval_physical(gmres_rhs_v(n,:),1,1,the_bc_tower%bc_tower_array(n))
    end do

    ! call gmres to compute delta v and delta pi
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dpi,rhotot_fc_new, &
               eta,eta_ed,kappa,theta_alpha,norm_pre_rhs)

    gmres_abs_tol = gmres_abs_tol_in ! Restore the desired tolerance   
       
    ! compute v^{n+1} = vbar^{n+1,*} + dumac
    ! compute pi^{n+1/2}= pi^{n+1/2,*} + dpi
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

    ! restore eta and kappa
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,2.d0,1,eta(n)%ng)
       call multifab_mult_mult_s_c(kappa(n),1,2.d0,1,kappa(n)%ng)
       do i=1,size(eta_ed,dim=2)
          call multifab_mult_mult_s_c(eta_ed(n,i),1,2.d0,1,eta_ed(n,i)%ng)
       end do
    end do

    ! if writing a plotfile for electro-explicit option, compute Epot using new densities
    ! use existing machinery to compute all mass fluxes - yes this is overkill
    ! but better than writing a separate routine for now
    ! diff/stoch_mass_flux and fluxdiv are not persistent so changing values doesn't
    ! matter.  grad_Epot is persistent so we overwrite the old values,
    ! which are thrown away
    if (use_charged_fluid .and. (.not. electroneutral) ) then
    if (plot_int.gt.0 .and. ( (mod(istep,plot_int).eq.0) .or. (istep.eq.max_step)) ) then
       ! compute the new charge and store it in charge_old (it could get modified to
       ! subtract off the average in compute_mass_fluxdiv)
       call dot_with_z(mla,rho_new,charge_old)
       call compute_mass_fluxdiv(mla,rho_new,rhotot_new,gradp_baro,Temp, &
                                 diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                 diff_mass_flux,stoch_mass_flux, &
                                 dt,time,dx,weights,the_bc_tower, &
                                 charge_old,grad_Epot_old,Epot, &
                                 permittivity)
    end if
    end if

    do n=1,nlevs
       call multifab_destroy(adv_mass_fluxdiv(n))
       call multifab_destroy(gmres_rhs_p(n))
       call multifab_destroy(dpi(n))
       do i=1,dm
          call multifab_destroy(umac_old(n,i))
          call multifab_destroy(mtemp(n,i))
          call multifab_destroy(adv_mom_fluxdiv_old(n,i))
          call multifab_destroy(adv_mom_fluxdiv_new(n,i))
          call multifab_destroy(diff_mom_fluxdiv_old(n,i))
          call multifab_destroy(diff_mom_fluxdiv_new(n,i))
          call multifab_destroy(stoch_mom_fluxdiv(n,i))
          call multifab_destroy(gmres_rhs_v(n,i))
          call multifab_destroy(dumac(n,i))
          call multifab_destroy(gradpi(n,i))
          call multifab_destroy(rho_fc(n,i))
          call multifab_destroy(rhotot_fc_old(n,i))
          call multifab_destroy(rhotot_fc_new(n,i))
          call multifab_destroy(diff_mass_flux(n,i))
          call multifab_destroy(mom_grav_force(n,i))
       end do
    end do

    if (nreactions > 0) then
       do n=1,nlevs
          call multifab_destroy(chem_rate_temp(n))
          call multifab_destroy(n_old(n))
          call multifab_destroy(n_new(n))
       end do
    end if

    if (variance_coef_mass .ne. 0.d0 .and. midpoint_stoch_mass_flux_type .eq. 2) then
       do n=1,nlevs
          call multifab_destroy(stoch_mass_fluxdiv_old(n))
          do i=1,dm
             call multifab_destroy(stoch_mass_flux_old(n,i))
          end do
       end do
    end if

    if (use_charged_fluid) then
       do n=1,nlevs
          do i=1,dm
             call multifab_destroy(Lorentz_force(n,i))
          end do
       end do
    end if

    if (use_multiphase) then
       do n=1,nlevs
          do i=1,dm
             call multifab_destroy(div_reversible_stress(n,i))
          end do
       end do
    end if

    if (present(total_mass_flux)) then
       do n=1,nlevs
          do i=1,dm
             call multifab_destroy(adv_mass_flux(n,i))
          end do
       end do
    end if

  end subroutine advance_timestep_bousq

end module advance_timestep_bousq_module
