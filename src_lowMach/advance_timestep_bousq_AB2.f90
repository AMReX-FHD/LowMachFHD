module advance_timestep_bousq_AB2_module

  use ml_layout_module
  use define_bc_module
  use bc_module
  use convert_stag_module
  use mk_advective_s_fluxdiv_module
  use diffusive_m_fluxdiv_module
  use stochastic_m_fluxdiv_module
  use stochastic_mass_fluxdiv_module
  use electrodiffusive_mass_fluxdiv_module
  use compute_mass_fluxdiv_module
  use compute_HSE_pres_module
  use convert_m_to_umac_module
  use convert_rhoc_to_c_module
  use mk_advective_m_fluxdiv_module
  use reservoir_bc_fill_module
  use gmres_module
  use div_and_grad_module
  use mk_grav_force_module
  use compute_mixture_properties_module
  use mass_flux_utilities_module
  use multifab_physbc_module
  use multifab_physbc_stag_module
  use zero_edgeval_module
  use fill_rho_ghost_cells_module
  use ml_solve_module
  use bndry_reg_module
  use bl_rng_module
  use bl_random_module
  use probin_common_module, only: advection_type, grav, rhobar, variance_coef_mass, &
                                  variance_coef_mom, barodiffusion_type, project_eos_int, &
                                  use_bl_rng, nspecies
  use probin_multispecies_module, only: midpoint_stoch_mass_flux_type
  use probin_gmres_module, only: gmres_abs_tol, gmres_rel_tol, mg_verbose
  use probin_charged_module, only: use_charged_fluid
  use probin_chemistry_module, only: nreactions

  use fabio_module

  implicit none

  private

  public :: advance_timestep_bousq_AB2

contains

  subroutine advance_timestep_bousq_AB2(mla,umac,rho_old,rho_new,rho0,adv_mom_fluxdiv_nm1, &
                                        gradp_baro,pi,eta,eta_ed,kappa,Temp,Temp_ed, &
                                        Epot_mass_fluxdiv, &
                                        diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                        stoch_mass_flux, &
                                        dx,dt,time,the_bc_tower,istep, &
                                        grad_Epot_old,grad_Epot_new,charge_old,charge_new, &
                                        Epot,permittivity)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: umac(:,:)             ! persistent - enters as v^n, leaves as v^{n+1}
    type(multifab) , intent(inout) :: rho_old(:)
    type(multifab) , intent(inout) :: rho_new(:)
    real(kind=dp_t), intent(in   ) :: rho0
    type(multifab) , intent(inout) :: adv_mom_fluxdiv_nm1(:,:) ! persistent - enters as ()^{n-1}, leaves as ()^n
    type(multifab) , intent(inout) :: gradp_baro(:,:)
    type(multifab) , intent(inout) :: pi(:)
    type(multifab) , intent(inout) :: eta(:)                ! not persistent
    ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(inout) :: eta_ed(:,:)           ! not persistent
    type(multifab) , intent(inout) :: kappa(:)              ! not persistent
    type(multifab) , intent(inout) :: Temp(:)               
    ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(inout) :: Temp_ed(:,:)          
    type(multifab) , intent(inout) :: Epot_mass_fluxdiv(:)  ! not persistent
    type(multifab) , intent(inout) :: diff_mass_fluxdiv(:)  ! not persistent
    type(multifab) , intent(inout) :: stoch_mass_fluxdiv(:) ! not persistent
    type(multifab) , intent(inout) :: stoch_mass_flux(:,:)  ! not persistent
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt,time
    type(bc_tower) , intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: istep
    type(multifab) , intent(inout) :: grad_Epot_old(:,:)
    type(multifab) , intent(inout) :: grad_Epot_new(:,:)
    type(multifab) , intent(inout) :: charge_old(:)
    type(multifab) , intent(inout) :: charge_new(:)
    type(multifab) , intent(inout) :: Epot(:)
    type(multifab) , intent(inout) :: permittivity(:)

    ! local
    type(multifab) ::           rhotot(mla%nlevel)
    type(multifab) :: adv_mass_fluxdiv(mla%nlevel)
    type(multifab) ::      gmres_rhs_p(mla%nlevel)
    type(multifab) ::              dpi(mla%nlevel)

    type(multifab) ::             mtemp(mla%nlevel,mla%dim)
    type(multifab) :: adv_mom_fluxdiv_n(mla%nlevel,mla%dim)
    type(multifab) ::  diff_mom_fluxdiv(mla%nlevel,mla%dim)
    type(multifab) :: stoch_mom_fluxdiv(mla%nlevel,mla%dim)
    type(multifab) ::       gmres_rhs_v(mla%nlevel,mla%dim)
    type(multifab) ::             dumac(mla%nlevel,mla%dim)
    type(multifab) ::            gradpi(mla%nlevel,mla%dim)
    type(multifab) ::            rho_fc(mla%nlevel,mla%dim)
    type(multifab) ::         rhotot_fc(mla%nlevel,mla%dim)
    type(multifab) ::           rho0_fc(mla%nlevel,mla%dim)
    type(multifab) ::    diff_mass_flux(mla%nlevel,mla%dim)
    type(multifab) ::    mom_grav_force(mla%nlevel,mla%dim)

    ! only used when variance_coef_mass>0 and midpoint_stoch_mass_flux_type=2
    type(multifab) :: stoch_mass_fluxdiv_old(mla%nlevel)
    type(multifab) :: stoch_mass_flux_old   (mla%nlevel,mla%dim)

    integer :: i,dm,n,nlevs,comp

    real(kind=dp_t) :: theta_alpha, norm_pre_rhs, gmres_abs_tol_in

    real(kind=dp_t) :: weights(2), sum

    weights(1) = 1.d0
    weights(2) = 0.d0

    if (use_charged_fluid) then
       call bl_error("advance_timestep_bousq_AB2 does not support charges yet")
    end if

    if (nreactions .gt. 0) then
       call bl_error("advance_timestep_bousq_AB2 does not support reactions yet")
    end if

    if (barodiffusion_type .ne. 0) then
       call bl_error("advance_timestep_bousq_AB2: barodiffusion not supported yet")
    end if
       
    nlevs = mla%nlevel
    dm = mla%dim

    theta_alpha = 1.d0/dt
    
    do n=1,nlevs
       call multifab_build(          rhotot(n),mla%la(n),1       ,rho_old(n)%ng)
       call multifab_build(adv_mass_fluxdiv(n),mla%la(n),nspecies,0)
       call multifab_build(     gmres_rhs_p(n),mla%la(n),1       ,0)
       call multifab_build(             dpi(n),mla%la(n),1       ,1)
       do i=1,dm
          call multifab_build_edge(            mtemp(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(adv_mom_fluxdiv_n(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge( diff_mom_fluxdiv(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(stoch_mom_fluxdiv(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(      gmres_rhs_v(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(            dumac(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(           gradpi(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(          rho0_fc(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(           rho_fc(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge(        rhotot_fc(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(   diff_mass_flux(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge(   mom_grav_force(n,i),mla%la(n),1       ,0,i)
       end do
    end do

    if (variance_coef_mass .ne. 0.d0 .and. midpoint_stoch_mass_flux_type .eq. 2) then
       do n=1,nlevs
          call multifab_build(stoch_mass_fluxdiv_old(n),mla%la(n),nspecies,0)
          do i=1,dm
             call multifab_build_edge(stoch_mass_flux_old(n,i),mla%la(n),nspecies,0,i)
          end do
       end do
    end if

    ! create a multifab rho0_fc for use in gmres solver
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(rho0_fc(n,i),rho0,all=.true.)
       end do
    end do

    ! average rho_i^n to faces
    call average_cc_to_face(nlevs,rho_old,rho_fc,1,c_bc_comp,nspecies,the_bc_tower%bc_tower_array)

    ! compute advective mass flux at t^n
    call mk_advective_s_fluxdiv(mla,umac,rho_fc,adv_mass_fluxdiv,.false.,dx,1,nspecies)


    ! compute diffusive, stochastic, potential mass fluxes
    ! with barodiffusion and thermodiffusion
    ! this computes "-F = rho W chi [Gamma grad x... ]"
    call compute_mass_fluxdiv(mla,rho_new,rhotot,gradp_baro,Temp, &
                              diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                              diff_mass_flux,stoch_mass_flux, &
                              dt,time,dx,weights,the_bc_tower, &
                              charge_old,grad_Epot_old,Epot, &
                              permittivity)

    ! FIXME compute reactions at t^n
    !
    !

    ! compute rho_i^{n+1/2} (store in rho_new)
    do n=1,nlevs
       call multifab_copy_c(rho_new(n),1,rho_old(n),1,nspecies,0)
       call multifab_saxpy_3_cc(rho_new(n),1,0.5d0*dt, adv_mass_fluxdiv(n),1,nspecies)
       call multifab_saxpy_3_cc(rho_new(n),1,0.5d0*dt,diff_mass_fluxdiv(n),1,nspecies)
       if (variance_coef_mass .ne. 0.d0) then
          call multifab_saxpy_3_cc(rho_new(n),1,0.5d0*dt,stoch_mass_fluxdiv(n),1,nspecies)
       end if
    end do

    ! compute rhotot^{n+1/2} from rho^{n+1/2} in VALID REGION
    call compute_rhotot(mla,rho_new,rhotot)

    ! fill rho and rhotot ghost cells
    call fill_rho_rhotot_ghost(mla,rho_new,rhotot,dx,the_bc_tower)

    ! average rho_i^{n+1/2} to faces
    call average_cc_to_face(nlevs,rho_new,rho_fc,1,c_bc_comp,nspecies,the_bc_tower%bc_tower_array)

    ! compute (eta,kappa)^{n+1/2}
    call compute_eta_kappa(mla,eta,eta_ed,kappa,rho_new,rhotot,Temp,dx, &
                           the_bc_tower%bc_tower_array)

    !!!!!!!!!
    ! set up GMRES solve for v^{n+1} and pi^{n+1/2}

    ! compute mtemp = (rho0*v)^n
    call convert_m_to_umac(mla,rho0_fc,mtemp,umac,.false.)

    ! build gmres_rhs_v
    ! first set gmres_rhs_v = mtemp / dt
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(gmres_rhs_v(n,i),1,mtemp(n,i),1,1,0)
          call multifab_mult_mult_s_c(gmres_rhs_v(n,i),1,1.d0/dt,1,0)
       end do
    end do

    ! compute adv_mom_fluxdiv_n
    call mk_advective_m_fluxdiv(mla,umac,mtemp,adv_mom_fluxdiv_n,.false., &
                                dx,the_bc_tower%bc_tower_array)

    ! add momentum advection terms to gmres_rhs_v
    if (istep .eq. 1) then
       ! first time step, use forward Euler
       call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,-1.d0,adv_mom_fluxdiv_n(n,i),1,1)
    else
       ! use AB2
       do n=1,nlevs
          do i=1,dm
             call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1, 1.5d0,adv_mom_fluxdiv_n  (n,i),1,1)
             call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,-0.5d0,adv_mom_fluxdiv_nm1(n,i),1,1)
          end do
       end do
    end if

    ! compute diff_mom_fluxdiv = L_0^n v^n
    call diffusive_m_fluxdiv(mla,diff_mom_fluxdiv,.false.,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    ! add (1/2)*diff_mom_fluxdiv to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,diff_mom_fluxdiv(n,i),1,1)
       end do
    end do


    if (variance_coef_mom .ne. 0.d0) then

       ! fill the stochastic momentum multifabs with new sets of random numbers
       call fill_m_stochastic(mla)

       ! compute stoch_mom_fluxdiv = div (sqrt() (W + W^T)^{n:n+1})
       ! we save this so we can re-use this term for the predictor/corrector first time step
       call stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,stoch_mom_fluxdiv,.false., &
                                 eta,eta_ed,Temp,Temp_ed,dx,dt,weights)

       ! add stochastic momentum fluxes to gmres_rhs_v
       do n=1,nlevs
          do i=1,dm
             call multifab_plus_plus_c(gmres_rhs_v(n,i),1,stoch_mom_fluxdiv(n,i),1,1,0)
          end do
       end do

    end if

    ! gravity (rho^{n+1/2}-rho0) * g
    if (any(grav(1:dm) .ne. 0.d0)) then

       ! put rho^{n+1/2} on faces (store in mom_grav_force)
       call average_cc_to_face(nlevs,rhotot,mom_grav_force,1,scal_bc_comp,1,the_bc_tower%bc_tower_array)
       do n=1,nlevs
          do i=1,dm
             ! compute (rho^{n+1/2}-rho0)*g
             call multifab_sub_sub_s(mom_grav_force(n,i),rho0)
             call multifab_mult_mult_s(mom_grav_force(n,i),grav(i))
          end do
       end do

       ! add gravity force to gmres_rhs_v
       do n=1,nlevs
          do i=1,dm
             call multifab_plus_plus_c(gmres_rhs_v(n,i),1,mom_grav_force(n,i),1,1,0)
          end do
       end do

    end if

    ! compute grad pi^n
    call compute_grad(mla,pi,gradpi,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

    ! subtract grad pi^n from gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradpi(n,i),1,1,0)
       end do
    end do

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

    ! compute mtemp = (rho0*vbar)^n
    call convert_m_to_umac(mla,rho0_fc,mtemp,umac,.false.)

    ! add -(mtemp/dt) to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,-1.d0/dt,mtemp(n,i),1,1)
       end do
    end do

    ! compute diff_mom_fluxdiv = L_0^n vbar^n
    call diffusive_m_fluxdiv(mla,diff_mom_fluxdiv,.false.,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    ! add (1/2)*diff_mom_fluxdiv to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,diff_mom_fluxdiv(n,i),1,1)
       end do
    end do

    ! compute div(vbar^n) and store in gmres_rhs_p
    ! the sign convention is correct since we solve -div(delta v) = div(vbar^n)
    call compute_div(mla,umac,gmres_rhs_p,dx,1,1,1)

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

    do n=1,nlevs
       call zero_edgeval_physical(gmres_rhs_v(n,:),1,1,the_bc_tower%bc_tower_array(n))
    end do

    ! call gmres to compute delta v and delta pi
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dpi,rho0_fc, &
               eta,eta_ed,kappa,theta_alpha,norm_pre_rhs)

    ! restore eta and kappa
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,2.d0,1,eta(n)%ng)
       call multifab_mult_mult_s_c(kappa(n),1,2.d0,1,kappa(n)%ng)
       do i=1,size(eta_ed,dim=2)
          call multifab_mult_mult_s_c(eta_ed(n,i),1,2.d0,1,eta_ed(n,i)%ng)
       end do
    end do
       
    ! compute v^{n+1} = v^n + dumac
    ! compute pi^{n+1}= pi^n + dpi
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

    ! do corrector velocity solve if this is the first time step
    if (istep .eq. 1) then

       ! subtract (1/2)*adv_mom_fluxdiv from gmres_rhs_v


       ! compute adv_mom_fluxdiv = -rho*v^{n+1,*},v^{n+1,*}


       ! add (1/2)*adv_mom_fluxdiv to gmres_rhs_v




    end if

    ! copy momentum advective fluxes into the "nm1" multifab for use in the next time step
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(adv_mom_fluxdiv_nm1(n,i),1,adv_mom_fluxdiv_n(n,i),1,1,0)
       end do
    end do

    ! compute mass fluxes and reactions at t^{n+1/2}
    if (midpoint_stoch_mass_flux_type .eq. 1) then
       ! strato

       ! compute diffusive, stochastic, potential mass fluxes
       ! with barodiffusion and thermodiffusion
       ! this computes "-F = rho W chi [Gamma grad x... ]"
       weights(:) = 1.d0/sqrt(2.d0)
       call compute_mass_fluxdiv(mla,rho_new,rhotot,gradp_baro,Temp, &
                                 diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                 diff_mass_flux,stoch_mass_flux, &
                                 dt,time,dx,weights,the_bc_tower, &
                                 charge_new,grad_Epot_new,Epot, &
                                 permittivity)

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

       ! compute diffusive, stochastic, potential mass fluxes
       ! with barodiffusion and thermodiffusion
       ! this computes "-F = rho W chi [Gamma grad x... ]"
       weights(1) = 0.d0
       weights(2) = 1.d0
       call compute_mass_fluxdiv(mla,rho_new,rhotot,gradp_baro,Temp, &
                                 diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                 diff_mass_flux,stoch_mass_flux, &
                                 0.5d0*dt,time,dx,weights,the_bc_tower, &
                                 charge_new,grad_Epot_new,Epot, &
                                 permittivity)

       if (variance_coef_mass .ne. 0.d0) then
          ! add stoch_mass_fluxdiv_old to stoch_mass_fluxdiv and multiply by 1/2
          do n=1,nlevs
             call multifab_plus_plus_c(stoch_mass_fluxdiv(n),1,stoch_mass_fluxdiv_old(n),1,nspecies,0)
             call multifab_mult_mult_s_c(stoch_mass_fluxdiv(n),1,0.5d0,nspecies,0)
          end do
       end if

    end if

    ! compute rho_i^{n+1}
    do n=1,nlevs
       call multifab_copy_c(rho_new(n),1,rho_old(n),1,nspecies,0)
       call multifab_saxpy_3_cc(rho_new(n),1,dt, adv_mass_fluxdiv(n),1,nspecies)
       call multifab_saxpy_3_cc(rho_new(n),1,dt,diff_mass_fluxdiv(n),1,nspecies)
       if (variance_coef_mass .ne. 0.d0) then
          call multifab_saxpy_3_cc(rho_new(n),1,dt,stoch_mass_fluxdiv(n),1,nspecies)
       end if
    end do

    ! compute rhotot^{n+1} from rho^{n+1} in VALID REGION
    call compute_rhotot(mla,rho_new,rhotot)

    ! fill rho and rhotot ghost cells
    call fill_rho_rhotot_ghost(mla,rho_new,rhotot,dx,the_bc_tower)

    ! average rho^{n+1} to faces
    call average_cc_to_face(nlevs,rho_new,rho_fc,1,c_bc_comp,nspecies,the_bc_tower%bc_tower_array)

    do n=1,nlevs
       call multifab_destroy(rhotot(n))
       call multifab_destroy(adv_mass_fluxdiv(n))
       call multifab_destroy(gmres_rhs_p(n))
       call multifab_destroy(dpi(n))
       do i=1,dm
          call multifab_destroy(mtemp(n,i))
          call multifab_destroy(adv_mom_fluxdiv_n(n,i))
          call multifab_destroy(diff_mom_fluxdiv(n,i))
          call multifab_destroy(stoch_mom_fluxdiv(n,i))
          call multifab_destroy(gmres_rhs_v(n,i))
          call multifab_destroy(dumac(n,i))
          call multifab_destroy(gradpi(n,i))
          call multifab_destroy(rho_fc(n,i))
          call multifab_destroy(rhotot_fc(n,i))
          call multifab_destroy(rho0_fc(n,i))
          call multifab_destroy(diff_mass_flux(n,i))
          call multifab_destroy(mom_grav_force(n,i))
       end do
    end do

    if (variance_coef_mass .ne. 0.d0 .and. midpoint_stoch_mass_flux_type .eq. 2) then
       do n=1,nlevs
          call multifab_destroy(stoch_mass_fluxdiv_old(n))
          do i=1,dm
             call multifab_destroy(stoch_mass_flux_old(n,i))
          end do
       end do
    end if

  end subroutine advance_timestep_bousq_AB2

end module advance_timestep_bousq_AB2_module
