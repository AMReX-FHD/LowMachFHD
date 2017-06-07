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

  subroutine advance_timestep_inertial(mla,umac_old,umac_new,rho_old,rho_new,rhotot_old,rhotot_new, &
                                       rhoh_old,rhoh_new,Temp_old,Temp_new,p0_old,p0_new, &
                                       gradp_baro,pi,eta,eta_ed,kappa, &
                                       diff_mass_fluxdiv, &
                                       dx,dt,time,the_bc_tower,istep)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: umac_old(:,:)
    type(multifab) , intent(inout) :: umac_new(:,:)
    type(multifab) , intent(in   ) :: rho_old(:)
    type(multifab) , intent(inout) :: rho_new(:)
    type(multifab) , intent(in   ) :: rhotot_old(:)
    type(multifab) , intent(inout) :: rhotot_new(:)
    type(multifab) , intent(in   ) :: rhoh_old(:)
    type(multifab) , intent(inout) :: rhoh_new(:)
    type(multifab) , intent(in   ) :: Temp_new(:)
    type(multifab) , intent(inout) :: Temp_new(:)
    real(kind=dp_t), intent(in   ) :: p0_old
    real(kind=dp_t), intent(inout) :: p0_new
    type(multifab) , intent(in   ) :: gradp_baro(:,:) ! not used, but required argument
    type(multifab) , intent(inout) :: pi(:)
    ! eta and kappa need to enter consistent with old and leave consistent with new
    type(multifab) , intent(inout) :: eta(:)
    type(multifab) , intent(inout) :: eta_ed(:,:) ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(inout) :: kappa(:)
    type(multifab) , intent(inout) :: diff_mass_fluxdiv(:) 
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt,time
    type(bc_tower) , intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: istep

    ! local
    type(multifab) :: gmres_rhs_p(mla%nlevel)
    type(multifab) :: gmres_rhs_v(mla%nlevel,mla%dim)

    type(multifab) :: dpi  (mla%nlevel)
    type(multifab) :: dumac(mla%nlevel,mla%dim)

    type(multifab) :: deltaT             (mla%nlevel)
    type(multifab) :: deltaT_rhs         (mla%nlevel)
    type(multifab) :: deltaT_solver_alpha(mla%nlevel)
    type(multifab) :: deltaT_solver_beta (mla%nlevel,mla%dim)

    type(multifab) :: mom_old(mla%nlevel,mla%dim)
    type(multifab) :: mom_new(mla%nlevel,mla%dim)
    type(multifab) ::   mtemp(mla%nlevel,mla%dim)

    type(multifab) :: conc_old(mla%nlevel,mla%dim)
    type(multifab) :: conc_new(mla%nlevel,mla%dim)

    type(multifab) ::        rho_fc_old(mla%nlevel,mla%dim)
    type(multifab) ::        rho_fc_new(mla%nlevel,mla%dim)
    type(multifab) ::     rhotot_fc_old(mla%nlevel,mla%dim)
    type(multifab) ::     rhotot_fc_new(mla%nlevel,mla%dim)
    type(multifab) ::       rhoh_fc_old(mla%nlevel,mla%dim)
    type(multifab) ::       rhoh_fc_new(mla%nlevel,mla%dim)

    type(multifab) :: gradpi(mla%nlevel,mla%dim)

    type(multifab) :: Peos(mla%nlevel)

    type(multifab) :: diff_mass_flux(mla%nlevel,mla%dim)

    type(multifab) :: mass_update_old(mla%nlevel,mla%dim)
    type(multifab) :: mass_update_new(mla%nlevel,mla%dim)

    type(multifab) :: rhoh_update_old(mla%nlevel,mla%dim)
    type(multifab) :: rhoh_update_new(mla%nlevel,mla%dim)

    type(multifab) :: mom_adv_update_old(mla%nlevel,mla%dim)
    type(multifab) :: mom_adv_update_new(mla%nlevel,mla%dim)

    type(multifab) :: mom_diff_update_old(mla%nlevel,mla%dim)

    type(multifab) ::       Scorr(mla%nlevel)
    type(multifab) :: delta_Scorr(mla%nlevel)

    real(kind=dp_t) :: Sbar_old, Sbar_new
    real(kind=dp_t) :: thetabar_old, thetabar_new
    real(kind=dp_t) :: Scorrbar
    real(kind=dp_t) :: peosbar, driftbar


    integer :: i,dm,n,nlevs

    real(kind=dp_t) :: theta_alpha, norm_pre_rhs, gmres_abs_tol_in

    real(kind=dp_t) :: weights(1)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "advance_timestep_inertial")

    weights(1) = 1.d0

    nlevs = mla%nlevel
    dm = mla%dim

    theta_alpha = 1.d0/dt
    
    do n=1,nlevs
       call multifab_build(gmres_rhs_p(n),mla%la(n),1       ,0)
       call multifab_build(        dpi(n),mla%la(n),1       ,1)
       do i=1,dm
          call multifab_build_edge(         mom_old(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(           mtemp(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge( adv_mom_fluxdiv(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(diff_mom_fluxdiv(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(     gmres_rhs_v(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(           dumac(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(          gradpi(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(      rho_fc_old(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge(      rho_fc_new(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge(   rhotot_fc_old(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(   rhotot_fc_new(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(     rhoh_fc_old(n,i),mla%la(n),       1,0,i)
          call multifab_build_edge(     rhoh_fc_new(n,i),mla%la(n),       1,0,i)
          call multifab_build_edge(  diff_mass_flux(n,i),mla%la(n),nspecies,0,i)
       end do
    end do

    ! compute grad(pi^n)
    call compute_grad(mla,pi,gradpi,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

    ! compute conc_old
    call convert_rhoc_to_c(mla,rho_old,rhotot_old,conc_old,.true.) ! valid region only
    call fill_c_ghost_cells(mla,conc_old,dx,the_bc_tower)

    ! compute rhotot, rho, and rhoh on faces at t^n
    call average_cc_to_face(nlevs,rhotot_old,rhotot_fc_old,1,scal_bc_comp,       1,the_bc_tower%bc_tower_array)
    call average_cc_to_face(nlevs,   rho_old,   rho_fc_old,1,   c_bc_comp,nspecies,the_bc_tower%bc_tower_array)
    call average_cc_to_face(nlevs,  rhoh_old,  rhoh_fc_old,1,   h_bc_comp,       1,the_bc_tower%bc_tower_array)

    ! compute diffusive mass and rhoh flux divergences at t^n
    ! (store result in mass_update_old and rhoh_update_old)
    ! also compute eta_old and kappa_old
    call diffusive_fluxdiv(mla,rho,rhotot,gradp_baro,Temp,p0,eta,eta_ed,kappa,
                           mass_update_old,diff_mass_flux,rhoh_update_old, &
                           dx,the_bc_tower)

    ! compute S^n and theta^n (store them in delta_S_old and delta_theta_old)
    call compute_S_theta(mla,delta_S_old,delta_theta_old,mass_udpate_old, &
                         rhoh_update_old,conc_old,Temp_old,rhotot_old)

    ! compute Scorr^n = (rho_P / rho) (Peos - p0) / dt
    ! first compute Scorr = Peos
    call compute_p(mla,rhotot_old,Temp_old,conc_old,Scorr)
    do n=1,nlevs
       ! compute Scorr = Peos - p0
       call multifab_sub_sub_s_c(Scorr(n),1,p0_old,1,0)
       ! compute Scorr = (rho_P / rho) (Peos - p0) / dt
       call scale_deltaP(mla,Scorr,rhotot_old,Temp,conc_old,dt,1.d0)
    end do

    ! split S^n, theta^n, Scorr into average and perturbational pieces
    do n=1,nlevs
       Sbar_old     = multifab_sum_c(delta_S_old(n)    ,1,1) / total_volume
       thetabar_old = multifab_sum_c(delta_theta_old(n),1,1) / total_volume
       Scorrbar     = multifab_sum_c(Scorr(n),1,1) / total_volume
       call multifab_sub_sub_s_c(delta_S_old(n)    ,1,Sbar_old    ,1,0)
       call multifab_sub_sub_s_c(delta_theta_old(n),1,thetabar_old,1,0)
       call multifab_copy_c(delta_Scorr(n),1,Scorr(n),1,1,0)
       call multifab_sub_sub_s_c(delta_Scorr(n),1,Scorrbar,1,0)
    end do

    ! compute p0_update_old
    p0_update_old = (Sbar_old + Scorrbar)/thetabar_old

    ! add advective terms to mass_update_old
    ! mass_update_old = [-div(rho_i*v) + div(F)]^n
    call mk_advective_s_fluxdiv(mla,umac_old,rho_fc_old,mass_update_old,.true.,dx,1,nspecies)

    ! add advective terms to rhoh_update_old
    call mk_advective_s_fluxdiv(mla,umac_old,rhoh_fc_old,rhoh_update_old,.true.,dx,1,1)

    ! add dp0/dt (p0_update_old) to rhoh_update_old
    ! rhoh_update_old = [-div(rhoh*v) + p0_update + div(Q) + sum(div(hk*Fk)) + rho*Hext]^n
    do n=1,nlevs
       call multifab_plus_plus_s_c(rhoh_update_old(n),1,p0_update_old,1,0)
    end do

    ! compute old momentum
    call convert_m_to_umac(mla,rhotot_fc_old,mom_old,umac_old,.false.)

    ! compute -div(rho*v*v)
    call mk_advective_m_fluxdiv(mla,umac_old,mom_old,mom_adv_update_old,.false., &
                                dx,the_bc_tower%bc_tower_array)

    ! compute A_0^n * v^n
    call diffusive_m_fluxdiv(mla,mom_diff_update_old,.false.,umac_old, &
                             eta,eta_ed,kappa,dx,the_bc_tower%bc_tower_array)

    ! loop over dpdt_iters
    do k=1,dpdt_iters

       if (k .eq. 1) then
          ! if this is the first iteration, copy from the old
          
          do n=1,nlevs
             call multifab_copy_c(rho_new(n)   ,1,rho_old(n)   ,1,nspecies,rho_new(n)%ng)
             call multifab_copy_c(rhotot_new(n),1,rhotot_old(n),1,1       ,rhotot_new(n)%ng)
             call multifab_copy_c(rhoh_new(n)  ,1,rhoh_old(n)  ,1,1       ,rhoh_new(n)%ng)
             do i=1,dm
                call multifab_copy_c(umac_new(n,i)     ,1,umac_old(n,i)     ,1,1       ,umac_new(n,i)%ng)
                call multifab_copy_c(rhotot_fc_new(n,i),1,rhotot_fc_old(n,i),1,1       ,rhotot_fc_new(n,i)%ng)
                call multifab_copy_c(rho_fc_new(n,i)   ,1,rho_fc_old(n,i)   ,1,nspecies,rho_fc_new(n,i)%ng)
                call multifab_copy_c(rhoh_fc_new(n,i)  ,1,rhoh_fc_old(n,i)  ,1,1       ,rhoh_fc_new(n,i)%ng)
             end do
             call multifab_copy_c(eta_new(n)         ,1,eta_old(n)         ,1,1          ,1)
             call multifab_copy_c(kappa_new(n)       ,1,kappa_old(n)       ,1,1          ,1)
             call multifab_copy_c(mass_fluxdiv_new(n),1,mass_fluxdiv_old(n),1,nspecies   ,0)
             call multifab_copy_c(rhoh_fluxdiv_new(n),1,rhoh_fluxdiv_old(n),1,1          ,0)
             call multifab_copy_c(delta_S_new(n)     ,1,delta_S_old(n)     ,1,1          ,0)
             call multifab_copy_c(delta_theta_new(n) ,1,delta_theta_old(n) ,1,1          ,0)

             ! rhoh_update_new
          end do

          p0_new = p0_old
          p0_update_new = p0_update_old
          Sbar_new = Sbar_old
          thetabar_new = thetabar_old

       else
          ! otherwise, generate new stuff not done at the end of the dpdt_iters loop


          ! compute new-time densities on faces

          ! update densities to new time


          ! update total density


          ! update pressure


       end if




       ! loop over deltaT_iters for temperature solve
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
          call compute_cp(mla,cc_solver_alpha,conc_new,Temp)
          do n=1,nlevs
             call multifab_mult_mult_c(cc_solver_alpha(n),1,rhotot_new(n),1,1,0)
             call multifab_mult_mult_s_c(cc_solver_alpha(n),1,1.d0/dt,1,0)
          end do


       end do



    end do











!!!!!!!!!!!!!!!!!!!!!!
    
    ! add D^n and St^n to rho_update
    do n=1,nlevs
       call multifab_copy_c(rho_update(n),1,diff_mass_fluxdiv(n),1,nspecies,0)
    end do

    ! add A^n to rho_update
    ! compute A^n for scalars using centered advection
    call mk_advective_s_fluxdiv(mla,umac,rho_fc_old,rho_update,.true.,dx,1,nspecies)

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
    call average_cc_to_face(nlevs,   rho_new,   rho_fc_new,1,   c_bc_comp,nspecies,the_bc_tower%bc_tower_array)
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
                              diff_mass_fluxdiv, &
                              diff_mass_flux, &
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
    call mk_advective_s_fluxdiv(mla,umac,rho_fc_new,rho_update,.true.,dx,1,nspecies)

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
    call average_cc_to_face(nlevs,   rho_new,   rho_fc_new,1,   c_bc_comp,nspecies,the_bc_tower%bc_tower_array)
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
                              diff_mass_fluxdiv, &
                              diff_mass_flux, &
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
       call multifab_destroy(gmres_rhs_p(n))
       call multifab_destroy(dpi(n))
       do i=1,dm
          call multifab_destroy(mold(n,i))
          call multifab_destroy(mtemp(n,i))
          call multifab_destroy(adv_mom_fluxdiv(n,i))
          call multifab_destroy(diff_mom_fluxdiv(n,i))
          call multifab_destroy(gmres_rhs_v(n,i))
          call multifab_destroy(dumac(n,i))
          call multifab_destroy(umac_tmp(n,i))
          call multifab_destroy(rho_fc_old(n,i))
          call multifab_destroy(rho_fc_new(n,i))
          call multifab_destroy(rhotot_fc_old(n,i))
          call multifab_destroy(rhotot_fc_new(n,i))
          call multifab_destroy(rhoh_fc_old(n,i))
          call multifab_destroy(rhoh_fc_new(n,i))
          call multifab_destroy(gradpi(n,i))
          call multifab_destroy(diff_mass_flux(n,i))
       end do
    end do

    call destroy(bpt)

  end subroutine advance_timestep_inertial

end module advance_timestep_inertial_module
