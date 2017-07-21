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
  use probin_gmres_module, only: gmres_abs_tol, gmres_rel_tol, mg_verbose
  use probin_charged_module, only: use_charged_fluid
  use probin_chemistry_module, only: nreactions

  use fabio_module

  implicit none

  private

  public :: advance_timestep_bousq_AB2

contains

  subroutine advance_timestep_bousq_AB2(mla,umac,rho_old,rho_new,rhotot_old,rhotot_new, &
                                        gradp_baro,pi,eta,eta_ed,kappa,Temp,Temp_ed, &
                                        Epot_mass_fluxdiv, &
                                        diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                        stoch_mass_flux, &
                                        dx,dt,time,the_bc_tower,istep, &
                                        grad_Epot_old,grad_Epot_new,charge_old,charge_new, &
                                        Epot,permittivity)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: umac(:,:)             ! enters are v^n, leaves as v^{n+1}
    type(multifab) , intent(inout) :: rho_old(:)
    type(multifab) , intent(inout) :: rho_new(:)
    type(multifab) , intent(inout) :: rhotot_old(:)         ! store rho0 in here?
    type(multifab) , intent(inout) :: rhotot_new(:)         ! shore rho0 in here?
    type(multifab) , intent(inout) :: gradp_baro(:,:)
    type(multifab) , intent(inout) :: pi(:)
    type(multifab) , intent(inout) :: eta(:)                ! not persistent - compute t^n value in this routine
    ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(inout) :: eta_ed(:,:)           ! not persistent - compute t^n value in this routine
    type(multifab) , intent(inout) :: kappa(:)              ! not persistent - compute t^n value in this routine
    type(multifab) , intent(inout) :: Temp(:)               ! constant
    ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(inout) :: Temp_ed(:,:)          ! constant
    type(multifab) , intent(inout) :: Epot_mass_fluxdiv(:)  ! not persistent - compute t^n value in this routine
    type(multifab) , intent(inout) :: diff_mass_fluxdiv(:)  ! not persistent - compute t^n value in this routine
    type(multifab) , intent(inout) :: stoch_mass_fluxdiv(:) ! not persistent - compute t^n value in this routine
    type(multifab) , intent(inout) :: stoch_mass_flux(:,:)  
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
    type(multifab) ::  diff_mass_fluxdiv_old(mla%nlevel)
    type(multifab) :: stoch_mass_fluxdiv_old(mla%nlevel)
    type(multifab) ::   adv_mass_fluxdiv_old(mla%nlevel)
    type(multifab) ::   adv_mass_fluxdiv    (mla%nlevel)

    type(multifab) ::          gmres_rhs_p(mla%nlevel)
    type(multifab) ::                  dpi(mla%nlevel)

    type(multifab) ::                  mold(mla%nlevel,mla%dim)
    type(multifab) ::                 mtemp(mla%nlevel,mla%dim)
    type(multifab) ::   adv_mom_fluxdiv_new(mla%nlevel,mla%dim)
    type(multifab) ::   adv_mom_fluxdiv_old(mla%nlevel,mla%dim)
    type(multifab) ::  diff_mom_fluxdiv_new(mla%nlevel,mla%dim)
    type(multifab) ::  diff_mom_fluxdiv_old(mla%nlevel,mla%dim)
    type(multifab) :: stoch_mom_fluxdiv_new(mla%nlevel,mla%dim)
    type(multifab) :: stoch_mom_fluxdiv_old(mla%nlevel,mla%dim)
    type(multifab) ::           gmres_rhs_v(mla%nlevel,mla%dim)
    type(multifab) ::                 dumac(mla%nlevel,mla%dim)
    type(multifab) ::                gradpi(mla%nlevel,mla%dim)
    type(multifab) ::                rho_fc(mla%nlevel,mla%dim)
    type(multifab) ::             rhotot_fc(mla%nlevel,mla%dim)
    type(multifab) ::        diff_mass_flux(mla%nlevel,mla%dim)

    type(multifab) :: m_grav_force_old(mla%nlevel,mla%dim)
    type(multifab) :: m_grav_force_new(mla%nlevel,mla%dim)

    integer :: i,dm,n,nlevs,comp

    real(kind=dp_t) :: theta_alpha, norm_pre_rhs, gmres_abs_tol_in

    real(kind=dp_t) :: weights(1), sum

    weights(1) = 1.d0

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
       call multifab_build( diff_mass_fluxdiv_old(n),mla%la(n),nspecies,0) 
       call multifab_build(stoch_mass_fluxdiv_old(n),mla%la(n),nspecies,0) 
       call multifab_build(adv_mass_fluxdiv_old(n),mla%la(n),nspecies,0)
       call multifab_build(adv_mass_fluxdiv    (n),mla%la(n),nspecies,0)

       call multifab_build(gmres_rhs_p(n),mla%la(n),1       ,0)
       call multifab_build(        dpi(n),mla%la(n),1       ,1)
       do i=1,dm
          call multifab_build_edge(                 mold(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(                mtemp(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(  adv_mom_fluxdiv_new(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(  adv_mom_fluxdiv_old(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge( diff_mom_fluxdiv_new(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge( diff_mom_fluxdiv_old(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(stoch_mom_fluxdiv_new(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(stoch_mom_fluxdiv_old(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(          gmres_rhs_v(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(                dumac(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(               gradpi(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(            rhotot_fc(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(               rho_fc(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge(       diff_mass_flux(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge(     m_grav_force_old(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(     m_grav_force_new(n,i),mla%la(n),1       ,0,i)
       end do
    end do

    ! average rho_i^n and rho^n to faces
    call average_cc_to_face(nlevs,   rho_old,   rho_fc,1,   c_bc_comp,nspecies,the_bc_tower%bc_tower_array)
    call average_cc_to_face(nlevs,rhotot_old,rhotot_fc,1,scal_bc_comp,       1,the_bc_tower%bc_tower_array)

    ! make copies of 
    ! diff_mass_fluxdiv_old = div (rho W chi Gamma grad x + ...)^n
    ! stoch_mass_fluxdiv_old = div sqrt(...) B^n Z^{n:n+1}
    do n=1,nlevs
       call multifab_copy_c( diff_mass_fluxdiv_old(n),1, diff_mass_fluxdiv(n),1,nspecies,0)
       call multifab_copy_c(stoch_mass_fluxdiv_old(n),1,stoch_mass_fluxdiv(n),1,nspecies,0)
    end do

    ! compute adv_mass_fluxdiv_old = -div(rho v w)^n
    call mk_advective_s_fluxdiv(mla,umac,rho_fc,adv_mass_fluxdiv_old,.false.,dx,1,nspecies)

    ! compute diff_mom_fluxdiv_old = L_0^n v^n
    call diffusive_m_fluxdiv(mla,diff_mom_fluxdiv_old,.false.,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    if (variance_coef_mom .ne. 0.d0) then   
       ! fill the stochastic momentum multifabs with new sets of random numbers
       call fill_m_stochastic(mla)

       ! compute stoch_mom_fluxdiv_old = div (sqrt() (W + W^T)^{n:n+1})
       call stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,stoch_mom_fluxdiv_old,.false., &
                                 eta,eta_ed,Temp,Temp_ed,dx,dt,weights)
    end if

    ! compute "old" gravity force
    if (any(grav(1:dm) .ne. 0.d0)) then
       call mk_grav_force(mla,m_grav_force_old,.false.,rhotot_fc,rhotot_fc,the_bc_tower)
    end if

    ! compute momentum, mold = rho^n v^n
    call convert_m_to_umac(mla,rhotot_fc,mold,umac,.false.)

    ! compute advective flux divergence, adv_mom_fluxdiv_old = div(-rho v v)^n
    call mk_advective_m_fluxdiv(mla,umac,mold,adv_mom_fluxdiv_old,.false., &
                                dx,the_bc_tower%bc_tower_array)

    ! start building rho_new - first piece is
    ! (rho w)^n - (dt/2) div (F_a^n + F_d^n)
    do n=1,nlevs
       call multifab_copy_c(rho_new(n),1,rho_old(n),1,nspecies,0)
       call multifab_saxpy_3_cc(rho_new(n),1,0.5d0*dt, adv_mass_fluxdiv_old(n),1,nspecies)
       call multifab_saxpy_3_cc(rho_new(n),1,0.5d0*dt,diff_mass_fluxdiv_old(n),1,nspecies)
       if (variance_coef_mass .ne. 0.d0) then
          call multifab_saxpy_3_cc(rho_new(n),1,0.5d0*dt,stoch_mass_fluxdiv_old(n),1,nspecies)
       end if
    end do

    ! compute rhotot^{n+1,*} from rho^{n+1,*} in VALID REGION
    call compute_rhotot(mla,rho_new,rhotot_new)

    ! fill rho and rhotot ghost cells
    call fill_rho_rhotot_ghost(mla,rho_new,rhotot_new,dx,the_bc_tower)

    ! average (rho,rhotot)^{n+1,*} to faces
    call average_cc_to_face(nlevs,   rho_new,   rho_fc,1,   c_bc_comp,nspecies,the_bc_tower%bc_tower_array)
    call average_cc_to_face(nlevs,rhotot_new,rhotot_fc,1,scal_bc_comp,       1,the_bc_tower%bc_tower_array)

    ! compute rho^{n+1,*}*g
    if (any(grav(1:dm) .ne. 0.d0)) then
       call mk_grav_force(mla,m_grav_force_new,.false.,rhotot_fc,rhotot_fc,the_bc_tower)
    end if

    ! compute (eta,kappa)^{n+1,*}
    call compute_eta_kappa(mla,eta,eta_ed,kappa,rho_new,rhotot_new,Temp,dx, &
                           the_bc_tower%bc_tower_array)

    ! diff_mass_fluxdiv = -div(F) = div (rho W chi Gamma grad x + ...)^{n+1,*}
    ! stoch_mass_fluxdiv = -div(F) = div sqrt(...) B^{n+1,*} Z^{n:n+1}
    call compute_mass_fluxdiv(mla,rho_new,rhotot_new,gradp_baro,Temp, &
                              diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                              diff_mass_flux,stoch_mass_flux, &
                              dt,time,dx,weights,the_bc_tower)

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

    ! compute div(vbar^n) and store in gmres_rhs_p
    ! the sign convention is correct since we solve -div(delta v) = div(vbar^n)
    call compute_div(mla,umac,gmres_rhs_p,dx,1,1,1)

    ! compute mtemp = rho^{n+1,*} vbar^n
    call convert_m_to_umac(mla,rhotot_fc,mtemp,umac,.false.)

    ! build gmres_rhs_v
    ! first set gmres_rhs_v = (mold - mtemp) / dt
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(gmres_rhs_v(n,i),1,mold(n,i),1,1,0)
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,mtemp(n,i),1,1,0)
          call multifab_mult_mult_s_c(gmres_rhs_v(n,i),1,1.d0/dt,1,0)
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

    ! add advective momentum flux divergence, adv_mom_fluxdiv_old = div(-rho v v)^n
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,1.d0,adv_mom_fluxdiv_old(n,i),1,1)
       end do
    end do

    ! compute diff_mom_fluxdiv_new = L_0^{n+1,*} vbar, where vbar = v^n with t^{n+1} bc's
    call diffusive_m_fluxdiv(mla,diff_mom_fluxdiv_new,.false.,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    ! add (1/2) (L_0^n v^n + L_0^{n+1,*} vbar^n)
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,diff_mom_fluxdiv_old(n,i),1,1)
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,diff_mom_fluxdiv_new(n,i),1,1)
       end do
    end do

    ! compute stoch_mom_fluxdiv_new = div (sqrt() (W + W^T)^{n:n+1})
    ! (these should only differ from the t^n stochastic fluxdiv because of eta^{n+1,*})
    if (variance_coef_mom .ne. 0.d0) then
       call stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,stoch_mom_fluxdiv_new,.false., &
                                 eta,eta_ed,Temp,Temp_ed,dx,dt,weights)
    end if

    ! add (1/2) (t^n + t^{n+1,*}) stochastic fluxes
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,stoch_mom_fluxdiv_old(n,i),1,1)
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,stoch_mom_fluxdiv_new(n,i),1,1)
       end do
    end do

    ! add (1/2) (t^n + t^{n+1,*}) gravitational force
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_grav_force_old(n,i),1,1)
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_grav_force_new(n,i),1,1)
       end do
    end do

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

    gmres_abs_tol_in = gmres_abs_tol ! Save this; reset after corrector

    ! call gmres to compute delta v and delta pi
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dpi,rhotot_fc, &
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
       
    ! compute v^{n+1,*} = v^n + dumac
    ! compute pi^{n+1,*}= pi^n + dpi
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
    ! end of predictor
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute mtemp = rho^{n+1,*} v^{n+1,*}
    call convert_m_to_umac(mla,rhotot_fc,mtemp,umac,.false.)

    ! compute -div(rho v w)^{n+1,*}
    call mk_advective_s_fluxdiv(mla,umac,rho_fc,adv_mass_fluxdiv,.false.,dx,1,nspecies)

    ! compute advective flux divergence, adv_mom_fluxdiv_new = div(-rho v v)^{n+1,*}
    call mk_advective_m_fluxdiv(mla,umac,mtemp,adv_mom_fluxdiv_new,.false., &
                                dx,the_bc_tower%bc_tower_array)

    ! continue rho_new
    ! (rho w)^n - (dt/2) div (F_a^n + F_d^n + F_e^{n+1,*})
    !           - (dt/2) div (F_a^{n+1,*} + F_d^{n+1,*})
    do n=1,nlevs
       call multifab_saxpy_3_cc(rho_new(n),1,0.5d0*dt, adv_mass_fluxdiv(n),1,nspecies)
       call multifab_saxpy_3_cc(rho_new(n),1,0.5d0*dt,diff_mass_fluxdiv(n),1,nspecies)
       if (variance_coef_mass .ne. 0.d0) then
          call multifab_saxpy_3_cc(rho_new(n),1,0.5d0*dt,stoch_mass_fluxdiv(n),1,nspecies)
       end if
    end do

    ! compute rhotot^{n+1} from rho^{n+1} in VALID REGION
    call compute_rhotot(mla,rho_new,rhotot_new)

    ! fill rho and rhotot ghost cells
    call fill_rho_rhotot_ghost(mla,rho_new,rhotot_new,dx,the_bc_tower)

    ! average (rho,rhotot)^{n+1} to faces
    call average_cc_to_face(nlevs,   rho_new,   rho_fc,1,   c_bc_comp,nspecies,the_bc_tower%bc_tower_array)
    call average_cc_to_face(nlevs,rhotot_new,rhotot_fc,1,scal_bc_comp,       1,the_bc_tower%bc_tower_array)

    ! compute rho^{n+1}*g
    if (any(grav(1:dm) .ne. 0.d0)) then
       call mk_grav_force(mla,m_grav_force_new,.false.,rhotot_fc,rhotot_fc,the_bc_tower)
    end if

    ! compute (eta,kappa)^{n+1}
    call compute_eta_kappa(mla,eta,eta_ed,kappa,rho_new,rhotot_new,Temp,dx, &
                           the_bc_tower%bc_tower_array)

    ! fill the stochastic mass multifabs with new sets of random numbers
    if (variance_coef_mass .ne. 0.d0) then
       if (use_bl_rng) then
          call bl_rng_copy_engine(rng_eng_diffusion_chk,rng_eng_diffusion)
       end if
       call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)
    end if

    ! mass fluxes for the next time step
    ! diff_mass_fluxdiv = -div(F) = div (rho W chi Gamma grad x + ...)^{n+1}
    ! stoch_mass_fluxdiv = -div(F) = div sqrt(...) B^{n+1} Z^{n+1:n+2}
    call compute_mass_fluxdiv(mla,rho_new,rhotot_new,gradp_baro,Temp, &
                              diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                              diff_mass_flux,stoch_mass_flux, &
                              dt,time,dx,weights,the_bc_tower)

    ! modify umac to respect the boundary conditions we want after the next gmres solve
    ! thus when we add L_0^n vbar^{n+1,*} to gmres_rhs_v and add div vbar^n to gmres_rhs_p we
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

    ! compute div(vbar^{n+1,*}) and store in gmres_rhs_p
    ! the sign convention is correct since we solve -div(delta v) = div(vbar^n)
    call compute_div(mla,umac,gmres_rhs_p,dx,1,1,1)

    ! compute mtemp = rho^{n+1} vbar^{n+1,*}
    call convert_m_to_umac(mla,rhotot_fc,mtemp,umac,.false.)

    ! build gmres_rhs_v
    ! first set gmres_rhs_v = (mold - mtemp) / dt
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(gmres_rhs_v(n,i),1,mold(n,i),1,1,0)
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,mtemp(n,i),1,1,0)
          call multifab_mult_mult_s_c(gmres_rhs_v(n,i),1,1.d0/dt,1,0)
       end do
    end do

    ! compute grad pi^{n+1,*}
    call compute_grad(mla,pi,gradpi,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

    ! subtract grad pi^{n+1,*} from gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradpi(n,i),1,1,0)
       end do
    end do

    ! add advective momentum flux divergence, (1/2) (div(-rho v v)^n + div(-rho v v)^{n+1,*})
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,adv_mom_fluxdiv_old(n,i),1,1)
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,adv_mom_fluxdiv_new(n,i),1,1)
       end do
    end do

    ! compute diff_mom_fluxdiv_new = L_0^{n+1} vbar, where vbar = v^{n+1,*} with t^{n+1} bc's
    call diffusive_m_fluxdiv(mla,diff_mom_fluxdiv_new,.false.,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    ! add (1/2) (L_0^n v^n + L_0^{n+1} vbar^{n+1,*})
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,diff_mom_fluxdiv_old(n,i),1,1)
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,diff_mom_fluxdiv_new(n,i),1,1)
       end do
    end do

    if (variance_coef_mom .ne. 0.d0) then
       ! compute stoch_mom_fluxdiv_new = div (sqrt() (W + W^T)^{n:n+1})
       ! (these should only differ from the t^n stochastic fluxdiv because of eta^{n+1})
       call stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,stoch_mom_fluxdiv_new,.false., &
                                 eta,eta_ed,Temp,Temp_ed,dx,dt,weights)
    end if

    ! add (1/2) (t^n + t^{n+1}) stochastic fluxes
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,stoch_mom_fluxdiv_old(n,i),1,1)
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,stoch_mom_fluxdiv_new(n,i),1,1)
       end do
    end do

    ! add (1/2) (t^n + t^{n+1}) gravitational force
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_grav_force_old(n,i),1,1)
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_grav_force_new(n,i),1,1)
       end do
    end do

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
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dpi,rhotot_fc, &
               eta,eta_ed,kappa,theta_alpha,norm_pre_rhs)

    ! restore eta and kappa
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,2.d0,1,eta(n)%ng)
       call multifab_mult_mult_s_c(kappa(n),1,2.d0,1,kappa(n)%ng)
       do i=1,size(eta_ed,dim=2)
          call multifab_mult_mult_s_c(eta_ed(n,i),1,2.d0,1,eta_ed(n,i)%ng)
       end do
    end do
       
    ! compute v^{n+1} = v^{n+1,*} + dumac
    ! compute pi^{n+1}= pi^{n+1,*} + dpi
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


    gmres_abs_tol = gmres_abs_tol_in ! Restore the desired tolerance

    do n=1,nlevs
       call multifab_destroy(diff_mass_fluxdiv_old(n))
       call multifab_destroy(stoch_mass_fluxdiv_old(n))
       call multifab_destroy(adv_mass_fluxdiv_old(n))
       call multifab_destroy(adv_mass_fluxdiv(n))
       call multifab_destroy(gmres_rhs_p(n))
       call multifab_destroy(dpi(n))
       do i=1,dm
          call multifab_destroy(mold(n,i))
          call multifab_destroy(mtemp(n,i))
          call multifab_destroy(adv_mom_fluxdiv_new(n,i))
          call multifab_destroy(adv_mom_fluxdiv_old(n,i))
          call multifab_destroy(diff_mom_fluxdiv_new(n,i))
          call multifab_destroy(diff_mom_fluxdiv_old(n,i))
          call multifab_destroy(stoch_mom_fluxdiv_new(n,i))
          call multifab_destroy(stoch_mom_fluxdiv_old(n,i))
          call multifab_destroy(gmres_rhs_v(n,i))
          call multifab_destroy(dumac(n,i))
          call multifab_destroy(rhotot_fc(n,i))
          call multifab_destroy(gradpi(n,i))
          call multifab_destroy(rho_fc(n,i))
          call multifab_destroy(diff_mass_flux(n,i))
          call multifab_destroy(m_grav_force_old(n,i))
          call multifab_destroy(m_grav_force_new(n,i))
       end do
    end do

  end subroutine advance_timestep_bousq_AB2

end module advance_timestep_bousq_AB2_module
