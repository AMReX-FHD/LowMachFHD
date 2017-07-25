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
  use bl_random_module
  use probin_common_module, only: advection_type, grav, rhobar, variance_coef_mass, &
                                  variance_coef_mom, barodiffusion_type, project_eos_int, &
                                  use_bl_rng, nspecies
  use probin_multispecies_module, only: midpoint_stoch_mass_flux_type
  use probin_charged_module, only: use_charged_fluid
  use probin_chemistry_module, only: nreactions

  use fabio_module

  implicit none

  private

  public :: advance_timestep_bousq

contains

  subroutine advance_timestep_bousq(mla,umac,rho_old,rho_new,rhotot_old,rhotot_new, &
                                    rho0,gradp_baro,pi,eta,eta_ed,kappa,Temp,Temp_ed, &
                                    Epot_mass_fluxdiv, &
                                    diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                    stoch_mass_flux, &
                                    dx,dt,time,the_bc_tower,istep, &
                                    grad_Epot_old,grad_Epot_new,charge_old,charge_new, &
                                    Epot,permittivity)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: umac(:,:)             ! enters as t^n, leaves as t^{n+1}
    type(multifab) , intent(inout) :: rho_old(:)
    type(multifab) , intent(inout) :: rho_new(:)
    type(multifab) , intent(inout) :: rhotot_old(:)
    type(multifab) , intent(inout) :: rhotot_new(:)
    real(kind=dp_t), intent(in   ) :: rho0
    type(multifab) , intent(inout) :: gradp_baro(:,:)
    type(multifab) , intent(inout) :: pi(:)
    type(multifab) , intent(inout) :: eta(:)                ! enters as t^n, leaves as t^{n+1}
    type(multifab) , intent(inout) :: eta_ed(:,:)           ! enters as t^n, leaves as t^{n+1}
    type(multifab) , intent(inout) :: kappa(:)              ! enters as t^n, leaves as t^{n+1}
    type(multifab) , intent(inout) :: Temp(:)
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
    type(multifab) , intent(inout) :: Epot(:)               ! not persistent
    type(multifab) , intent(inout) :: permittivity(:)       ! not persistent

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
    type(multifab) ::       diff_mass_flux(mla%nlevel,mla%dim)
    type(multifab) ::       mom_grav_force(mla%nlevel,mla%dim)

    ! only used when variance_coef_mass>0 and midpoint_stoch_mass_flux_type=2
    type(multifab) :: stoch_mass_fluxdiv_old(mla%nlevel)
    type(multifab) ::    stoch_mass_flux_old(mla%nlevel,mla%dim)

    ! need a face-centered multifab of rho0 for gmres solver
    type(multifab) :: rho0_fc(mla%nlevel,mla%dim)

    integer :: i,dm,n,nlevs

    real(kind=dp_t) :: theta_alpha, norm_pre_rhs
    real(kind=dp_t) :: weights(2)

    weights(1) = 1.d0
    weights(2) = 0.d0

    if (use_charged_fluid) then
       call bl_error("advance_timestep_bousq does not support charges yet")
    end if

    if (nreactions .gt. 0) then
       call bl_error("advance_timestep_bousq does not support reactions yet")
    end if

    if (barodiffusion_type .ne. 0) then
       call bl_error("advance_timestep_bousq: barodiffusion not supported yet")
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

    ! create a multifab rho0_fc for use in gmres solver
    do n=1,nlevs
       do i=1,dm
          call multifab_build_edge(rho0_fc(n,i),mla%la(n),1,1,i)
          call multifab_setval(rho0_fc(n,i),rho0,all=.true.)
       end do
    end do

    ! make a copy of umac at t^n
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(umac_old(n,i),1,umac(n,i),1,1,1)
       end do
    end do    

    ! Step 1: solve for v^{n+1,*} and pi^{n+1/2,*} using GMRES

    ! compute mtemp = (rho0*v)^n
    call convert_m_to_umac(mla,rho0_fc,mtemp,umac,.false.)

    ! set gmres_rhs_v = mtemp / dt
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(gmres_rhs_v(n,i),1,mtemp(n,i),1,1,0)
          call multifab_mult_mult_s_c(gmres_rhs_v(n,i),1,1.d0/dt,1,0)
       end do
    end do

    ! compute adv_mom_fluxdiv_old = -rho0*v^n*v^n
    ! save this for use in the corrector GMRES solve
    call mk_advective_m_fluxdiv(mla,umac,mtemp,adv_mom_fluxdiv_old,.false., &
                                dx,the_bc_tower%bc_tower_array)

    ! add -rho0*v^n,v^n to gmres_rhs_v
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
       call stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,stoch_mom_fluxdiv,.false., &
                                 eta,eta_ed,Temp,Temp_ed,dx,dt,weights)

       ! add stochastic momentum fluxes to gmres_rhs_v
       do n=1,nlevs
          do i=1,dm
             call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,1.d0,stoch_mom_fluxdiv(n,i),1,1)
          end do
       end do

    end if

    ! gravity (rho^n-rho0) * g
    if (any(grav(1:dm) .ne. 0.d0)) then

       ! put rho^n on faces (store in mom_grav_force)
       call average_cc_to_face(nlevs,rhotot_old,mom_grav_force,1,scal_bc_comp,1,the_bc_tower%bc_tower_array)
       do n=1,nlevs
          do i=1,dm
             ! compute (rho^n-rho0)*g
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

    ! compute grad pi^{n-1/2}
    call compute_grad(mla,pi,gradpi,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

    ! subtract grad pi^{n-1/2} from gmres_rhs_v
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

    do n=1,nlevs
       call zero_edgeval_physical(gmres_rhs_v(n,:),1,1,the_bc_tower%bc_tower_array(n))
    end do

    ! call gmres to compute delta v and delta pi
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dpi,rho0_fc, &
               eta,eta_ed,kappa,theta_alpha,norm_pre_rhs)
       
    ! compute v^{n+1} = v^n + dumac
    ! no need to update pi yet
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(umac(n,i),1,dumac(n,i),1,1,0)
       end do
    end do
       
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

    ! restore eta and kappa
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,2.d0,1,eta(n)%ng)
       call multifab_mult_mult_s_c(kappa(n),1,2.d0,1,kappa(n)%ng)
       do i=1,size(eta_ed,dim=2)
          call multifab_mult_mult_s_c(eta_ed(n,i),1,2.d0,1,eta_ed(n,i)%ng)
       end do
    end do

    ! Step 2: compute mass fluxes and reactions at t^n

    ! average rho_i^n to faces
    call average_cc_to_face(nlevs,rho_old,rho_fc,1,c_bc_comp,nspecies,the_bc_tower%bc_tower_array)

    ! compute adv_mass_fluxdiv = -rho_i^n * v^n
    call mk_advective_s_fluxdiv(mla,umac_old,rho_fc,adv_mass_fluxdiv,.false.,dx,1,nspecies)

    ! increment adv_mass_fluxdiv by -rho_i^n * v^{n+1,*}
    call mk_advective_s_fluxdiv(mla,umac,rho_fc,adv_mass_fluxdiv,.true.,dx,1,nspecies)

    ! compute diffusive, stochastic, potential mass fluxes
    ! with barodiffusion and thermodiffusion
    ! this computes "-F = rho W chi [Gamma grad x... ]"
    call compute_mass_fluxdiv(mla,rho_old,rhotot_old,gradp_baro,Temp, &
                              diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                              diff_mass_flux,stoch_mass_flux, &
                              dt,time,dx,weights,the_bc_tower, &
                              charge_old,grad_Epot_old,Epot, &
                              permittivity)

    ! FIXME compute reactions at t^n
    !
    !

    ! Step 3: density prediction to t^{n+1/2}

    ! compute rho_i^{n+1/2} (store in rho_new)
    ! multiply adv_mass_fluxdiv by (1/4) since it contains -rho_i^n * (v^n + v^{n+1,*})
    do n=1,nlevs
       call multifab_copy_c(rho_new(n),1,rho_old(n),1,nspecies,0)
       call multifab_saxpy_3_cc(rho_new(n),1,0.25d0*dt, adv_mass_fluxdiv(n),1,nspecies)
       call multifab_saxpy_3_cc(rho_new(n),1, 0.5d0*dt,diff_mass_fluxdiv(n),1,nspecies)
       if (variance_coef_mass .ne. 0.d0) then
          call multifab_saxpy_3_cc(rho_new(n),1,0.5d0*dt,stoch_mass_fluxdiv(n),1,nspecies)
       end if
    end do

    ! compute rhotot^{n+1/2} from rho^{n+1/2} in VALID REGION
    call compute_rhotot(mla,rho_new,rhotot_new)

    ! fill rho and rhotot ghost cells
    call fill_rho_rhotot_ghost(mla,rho_new,rhotot_new,dx,the_bc_tower)

    ! compute (eta,kappa)^{n+1/2}
    call compute_eta_kappa(mla,eta,eta_ed,kappa,rho_new,rhotot_new,Temp,dx, &
                           the_bc_tower%bc_tower_array)

    ! Step 4: compute mass fluxes and reactions at t^{n+1/2}

    ! average rho_i^{n+1/2} to faces
    call average_cc_to_face(nlevs,rho_new,rho_fc,1,c_bc_comp,nspecies,the_bc_tower%bc_tower_array)

    ! compute adv_mass_fluxdiv = -rho_i^{n+1/2} * v^n
    call mk_advective_s_fluxdiv(mla,umac_old,rho_fc,adv_mass_fluxdiv,.false.,dx,1,nspecies)

    ! increment adv_mass_fluxdiv by -rho_i^{n+1/2} * v^{n+1,*}
    call mk_advective_s_fluxdiv(mla,umac,rho_fc,adv_mass_fluxdiv,.true.,dx,1,nspecies)


    ! compute mass fluxes and reactions at t^{n+1/2}
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
       call compute_mass_fluxdiv(mla,rho_new,rhotot_new,gradp_baro,Temp, &
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

    ! FIXME compute reactions at t^{n+1}
    !
    !

    ! Step 5: density integration to t^{n+1}

    ! compute rho_i^{n+1}
    ! multiply adv_mass_fluxdiv by (1/2) since it contains -rho_i^{n+1/2} * (v^n + v^{n+1,*})
    do n=1,nlevs
       call multifab_copy_c(rho_new(n),1,rho_old(n),1,nspecies,0)
       call multifab_saxpy_3_cc(rho_new(n),1,0.5d0*dt, adv_mass_fluxdiv(n),1,nspecies)
       call multifab_saxpy_3_cc(rho_new(n),1,      dt,diff_mass_fluxdiv(n),1,nspecies)
       if (variance_coef_mass .ne. 0.d0) then
          call multifab_saxpy_3_cc(rho_new(n),1,dt,stoch_mass_fluxdiv(n),1,nspecies)
       end if
    end do

    ! compute rhotot^{n+1} from rho^{n+1} in VALID REGION
    call compute_rhotot(mla,rho_new,rhotot_new)

    ! fill rho and rhotot ghost cells
    call fill_rho_rhotot_ghost(mla,rho_new,rhotot_new,dx,the_bc_tower)

    ! compute (eta,kappa)^{n+1}
    call compute_eta_kappa(mla,eta,eta_ed,kappa,rho_new,rhotot_new,Temp,dx, &
                           the_bc_tower%bc_tower_array)

    ! Step 6: solve for v^{n+1} and pi^{n+1/2} using GMRES

    ! set gmres_rhs_v = (rho0*v)^n / dt
    do n=1,nlevs
       do i=1,dm
          call convert_m_to_umac(mla,rho0_fc,gmres_rhs_v,umac_old,.false.)
          call multifab_mult_mult_s_c(gmres_rhs_v(n,i),1,1.d0/dt,1,0)
       end do
    end do

    ! compute mtemp = (rho0*v)^{n+1,*}
    call convert_m_to_umac(mla,rho0_fc,mtemp,umac,.false.)

    ! compute adv_mom_fluxdiv_new = -rho0*v^{n+1,*}*v^{n+1,*}
    call mk_advective_m_fluxdiv(mla,umac,mtemp,adv_mom_fluxdiv_new,.false., &
                                dx,the_bc_tower%bc_tower_array)

    ! add (1/2) (-rho0*v^n*v^n -rho0*v^{n+1,*}*v^{n+1,*}) to gmres_rhs_v
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
       call stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,stoch_mom_fluxdiv,.false., &
                                 eta,eta_ed,Temp,Temp_ed,dx,dt,weights)

       ! add stochastic momentum fluxes to gmres_rhs_v
       do n=1,nlevs
          do i=1,dm
             call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,stoch_mom_fluxdiv(n,i),1,1)
          end do
       end do

    end if

    ! gravity [ (rho^n+rho^{n+1})/2 - rho0 ] * g
    if (any(grav(1:dm) .ne. 0.d0)) then

       ! rhotot_old is no needed anymore, so we put (1/2)*(rho_old + rho_new) in it
       do n=1,nlevs
          call multifab_plus_plus_c(rhotot_old(n),1,rhotot_new(n),1,1,rhotot_old(n)%ng)
          call multifab_mult_mult_s_c(rhotot_old(n),1,0.5d0,1,rhotot_old(n)%ng)
       end do

       ! put (rho^n+rho^{n+1})/2 on faces (store in mom_grav_force)
       call average_cc_to_face(nlevs,rhotot_old,mom_grav_force,1,scal_bc_comp,1,the_bc_tower%bc_tower_array)
       do n=1,nlevs
          do i=1,dm
             ! compute [(rho^n+rho^{n+1})/2-rho0]*g
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

    ! subtract grad pi^{n-1/2} from gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradpi(n,i),1,1,0)
       end do
    end do

    ! copy umac_old back into umac
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(umac(n,i),1,umac_old(n,i),1,1,1)
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

    ! compute diff_mom_fluxdiv_new = L_0^{n+1} vbar^n
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

    do n=1,nlevs
       call zero_edgeval_physical(gmres_rhs_v(n,:),1,1,the_bc_tower%bc_tower_array(n))
    end do

    ! call gmres to compute delta v and delta pi
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dpi,rho0_fc, &
               eta,eta_ed,kappa,theta_alpha,norm_pre_rhs)
       
    ! compute v^{n+1} = vbar^n + dumac
    ! compute pi^{n+1/2}= pi^{n-1/2} + dpi
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

    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(rho0_fc(n,i))
       end do
    end do

  end subroutine advance_timestep_bousq

end module advance_timestep_bousq_module
