module advance_timestep_inertial_module

  use ml_layout_module
  use define_bc_module
  use bc_module
  use convert_stag_module
  use mk_advective_s_fluxdiv_module
  use diffusive_m_fluxdiv_module
  use stochastic_m_fluxdiv_module
  use stochastic_mass_fluxdiv_module
  use compute_mass_fluxdiv_module
  use reservoir_bc_fill_module
  use bds_module
  use gmres_module
  use div_and_grad_module
  use eos_check_module
  use mk_grav_force_module
  use compute_mixture_properties_module
  use mass_flux_utilities_module
  use multifab_physbc_module
  use multifab_physbc_stag_module
  use probin_common_module, only: advection_type, grav, rhobar, variance_coef_mass, &
                                  variance_coef_mom, restart, algorithm_type
  use probin_gmres_module, only: gmres_abs_tol, gmres_rel_tol
  use probin_multispecies_module, only: nspecies, rho_part_bc_comp
  use analysis_module

  implicit none

  private

  public :: advance_timestep_inertial

  ! special inhomogeneous boundary condition multifab
  ! vel_bc_n(nlevs,dm) are the normal velocities
  ! in 2D, vel_bc_t(nlevs,2) respresents
  !   1. y-velocity bc on x-faces (nodal)
  !   2. x-velocity bc on y-faces (nodal)
  ! in 3D, vel_bc_t(nlevs,6) represents
  !   1. y-velocity bc on x-faces (nodal in y and x)
  !   2. z-velocity bc on x-faces (nodal in z and x)
  !   3. x-velocity bc on y-faces (nodal in x and y)
  !   4. z-velocity bc on y-faces (nodal in z and y)
  !   5. x-velocity bc on z-faces (nodal in x and z)
  !   6. y-velocity bc on z-faces (nodal in y and z)
  type(multifab), allocatable, save :: vel_bc_n(:,:)
  type(multifab), allocatable, save :: vel_bc_t(:,:)

contains

  subroutine advance_timestep_inertial(mla,umac,rho_old,rho_new,rhotot_old,rhotot_new, &
                                         pres,eta,eta_ed,kappa,Temp,Temp_ed, &
                                         diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                         dx,dt,time,the_bc_tower,istep)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(inout) :: rho_old(:)
    type(multifab) , intent(inout) :: rho_new(:)
    type(multifab) , intent(inout) :: rhotot_old(:)
    type(multifab) , intent(inout) :: rhotot_new(:)
    type(multifab) , intent(inout) :: pres(:)
    ! eta and kappa need to enter consistent with old and leave consistent with new
    type(multifab) , intent(inout) :: eta(:)
    type(multifab) , intent(inout) :: eta_ed(:,:) ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(inout) :: kappa(:)
    type(multifab) , intent(inout) :: Temp(:)
    type(multifab) , intent(inout) :: Temp_ed(:,:) ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(inout) :: diff_mass_fluxdiv(:)
    type(multifab) , intent(inout) :: stoch_mass_fluxdiv(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt,time
    type(bc_tower) , intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: istep

    ! local
    type(multifab) ::  rho_update(mla%nlevel)
    type(multifab) ::   bds_force(mla%nlevel)
    type(multifab) :: gmres_rhs_p(mla%nlevel)
    type(multifab) ::          dp(mla%nlevel)
    type(multifab) ::        divu(mla%nlevel)

    type(multifab) :: gmres_rhs_v(mla%nlevel,mla%dim)
    type(multifab) ::       dumac(mla%nlevel,mla%dim)
    type(multifab) ::       gradp(mla%nlevel,mla%dim)
    type(multifab) ::      rho_fc(mla%nlevel,mla%dim)
    type(multifab) ::   rhotot_fc(mla%nlevel,mla%dim)
    type(multifab) :: flux_total(mla%nlevel,mla%dim)

    integer :: i,dm,n,nlevs

    real(kind=dp_t) :: theta_alpha, norm_pre_rhs, gmres_abs_tol_in

    real(kind=dp_t) :: weights(algorithm_type)

    weights(:) = 0.d0
    weights(1) = 1.d0

    nlevs = mla%nlevel
    dm = mla%dim

    theta_alpha = 0.d0
    
    call build_bc_multifabs(mla)
    
    do n=1,nlevs
       call multifab_build(   rho_update(n),mla%la(n),nspecies,0)
       call multifab_build(  bds_force(n),mla%la(n),nspecies,1)
       call multifab_build(gmres_rhs_p(n),mla%la(n),1,0)
       call multifab_build(         dp(n),mla%la(n),1,1)
       call multifab_build(       divu(n),mla%la(n),1,0)
       do i=1,dm
          call multifab_build_edge(gmres_rhs_v(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(      dumac(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(      gradp(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(     rho_fc(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge(  rhotot_fc(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge( flux_total(n,i),mla%la(n),nspecies,0,i)
       end do
    end do

    do n=1,nlevs
       call setval(rho_update(n),0.d0,all=.true.)
       call setval(bds_force(n),0.d0,all=.true.)
       do i=1,dm
          call setval(dumac(n,i),0.d0,all=.true.)
       end do
    end do

    ! average rho and rhotot to faces
    call average_cc_to_face(nlevs,   rho_old,   rho_fc,1,rho_part_bc_comp,nspecies,the_bc_tower%bc_tower_array)
    call average_cc_to_face(nlevs,rhotot_old,rhotot_fc,1,    scal_bc_comp,       1,the_bc_tower%bc_tower_array)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 1 - Predictor Stochastic/Diffusive Fluxes
    ! Step 2 - Predictor Stokes Solve
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! fill the stochastic multifabs with a new set of random numbers
    ! if this is the first step after initialization or restart then
    ! we already have random numbers from initialization
    if (istep .ne. 1 .and. istep .ne. restart+1) then
       call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)
       call fill_m_stochastic(mla)
    end if

    ! build up rhs_v for gmres solve
    do n=1,nlevs
       do i=1,dm
          call setval(gmres_rhs_v(n,i),0.d0,all=.true.)
       end do
    end do

    ! compute grad p^{n-1/2}
    call compute_grad(mla,pres,gradp,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

    ! subtract grad p^{n-1/2} from gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradp(n,i),1,1,0)
       end do
    end do

    ! add div(Sigma^(1)) to gmres_rhs_v
    if (variance_coef_mom .ne. 0.d0) then
       if (algorithm_type .eq. 1) then
          call stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,gmres_rhs_v, &
                                    eta,eta_ed,Temp,Temp_ed,dx,dt,weights)
       else if (algorithm_type .eq. 2) then
          call stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,gmres_rhs_v, &
                                    eta,eta_ed,Temp,Temp_ed,dx,0.5d0*dt,weights)
       end if
    end if

    ! add rho^n*g to gmres_rhs_v
    if (any(grav(1:dm) .ne. 0.d0)) then
       call mk_grav_force(mla,gmres_rhs_v,rhotot_fc,rhotot_fc,the_bc_tower)
    end if

    ! initialize rhs_p for gmres solve to zero
    do n=1,nlevs
       call setval(gmres_rhs_p(n),0.d0,all=.true.)
    end do

    ! reset inhomogeneous bc condition to deal with reservoirs
    call set_inhomogeneous_vel_bcs(mla,vel_bc_n,vel_bc_t,eta_ed,dx, &
                                   the_bc_tower%bc_tower_array)

    do n=1,nlevs
       do i=1,dm
          call setval(flux_total(n,i),0.d0,all=.true.)
       end do
    end do

    ! compute diffusive and stochastic mass fluxes
    ! this computes "-F" so we later multiply by -1
    if (algorithm_type .eq. 1) then
       call compute_mass_fluxdiv_wrapper(mla,rho_old,rhotot_old, &
                                         diff_mass_fluxdiv,stoch_mass_fluxdiv,Temp, &
                                         flux_total,dt,time,dx,weights, &
                                         the_bc_tower%bc_tower_array)
    else if (algorithm_type .eq. 2) then
       call compute_mass_fluxdiv_wrapper(mla,rho_old,rhotot_old, &
                                         diff_mass_fluxdiv,stoch_mass_fluxdiv,Temp, &
                                         flux_total,0.5d0*dt,time,dx,weights, &
                                         the_bc_tower%bc_tower_array)
    end if

    do n=1,nlevs
       call multifab_mult_mult_s_c(diff_mass_fluxdiv(n),1,-1.d0,nspecies,0)
       if (variance_coef_mass .ne. 0) then
          call multifab_mult_mult_s_c(stoch_mass_fluxdiv(n),1,-1.d0,nspecies,0)
       end if
       do i=1,dm
          call multifab_mult_mult_s_c(flux_total(n,i),1,-1.d0,nspecies,0)
       end do
    end do

    ! set the Dirichlet velocity value on reservoir faces
    call reservoir_bc_fill(mla,flux_total,vel_bc_n,the_bc_tower%bc_tower_array)

    do n=1,nlevs
       do i=1,nspecies
          call saxpy(gmres_rhs_p(n),1,-1.d0/rhobar(i), diff_mass_fluxdiv(n),i,1)
          if (variance_coef_mass .ne. 0.d0) then
             call saxpy(gmres_rhs_p(n),1,-1.d0/rhobar(i),stoch_mass_fluxdiv(n),i,1)
          end if
       end do
    end do

    ! reset rho_update for all scalars to zero
    ! then, set rho_update for rho1 to F^n = div(rho*chi grad c)^n + div(Psi^(1))
    do n=1,nlevs
       call multifab_setval_c(rho_update(n),0.d0,1,nspecies,all=.true.)
       ! add fluxes
       call multifab_plus_plus_c(rho_update(n),1, diff_mass_fluxdiv(n),1,nspecies)
       if (variance_coef_mass .ne. 0.d0) then
          call multifab_plus_plus_c(rho_update(n),1,stoch_mass_fluxdiv(n),1,nspecies)
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
    call diffusive_m_fluxdiv(mla,gmres_rhs_v,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    ! compute div v^{n-1/2}
    call compute_div(mla,umac,divu,dx,1,1,1)

    ! add div v^{n-1/2} to gmres_rhs_p
    ! now gmres_rhs_p = div v^{n-1/2} - S^n
    ! the sign convention is correct since we solve -div(delta v) = gmres_rhs_p
    do n=1,nlevs
       call multifab_plus_plus_c(gmres_rhs_p(n),1,divu(n),1,1,0)
    end do

    ! set the initial guess to zero
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(dumac(n,i),0.d0,all=.true.)
       end do
       call multifab_setval(dp(n),0.d0,all=.true.)
    end do

    gmres_abs_tol_in = gmres_abs_tol ! Save this  
    ! This relies entirely on relative tolerance and can fail if the rhs is roundoff error only:
    !gmres_abs_tol = 0.d0 ! It is better to set gmres_abs_tol in namelist to a sensible value

    ! call gmres to compute delta v and delta p
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dp,rhotot_fc, &
               eta,eta_ed,kappa,theta_alpha,norm_pre_rhs)

    ! for the corrector gmres solve we want the stopping criteria based on the
    ! norm of the preconditioned rhs from the predictor gmres solve.  otherwise
    ! for cases where du in the corrector should be small the gmres stalls
    gmres_abs_tol = max(gmres_abs_tol_in, norm_pre_rhs*gmres_rel_tol)

    ! compute v^* = v^{n-1/2} + delta v
    ! compute p^* = p^{n-1/2} + delta p
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(umac(n,i),1,dumac(n,i),1,1,0)
       end do
       call multifab_plus_plus_c(pres(n),1,dp(n),1,1,0)
    end do

    do n=1,nlevs
       ! presure ghost cells
       call multifab_fill_boundary(pres(n))
       call multifab_physbc(pres(n),1,pres_bc_comp,1,the_bc_tower%bc_tower_array(n),dx(n,:))
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
         call multifab_copy_c(bds_force(n),1,rho_update(n),1,nspecies,0)
         call multifab_fill_boundary(bds_force(n))
      end do

      if (advection_type .eq. 1 .or. advection_type .eq. 2) then
          call bds(mla,umac,rho_old,rho_update,bds_force,rho_fc,dx,dt,1,nspecies,the_bc_tower)
      else
          call bds_quad(mla,umac,rho_old,rho_update,bds_force,rho_fc,dx,dt,1,nspecies,the_bc_tower)
      end if
    else
       call mk_advective_s_fluxdiv(mla,umac,rho_fc,rho_update,dx,1,nspecies)
    end if

    ! compute s^{*,n+1/2} = s^n + (dt/2) * (A^n + F^n)
    ! store result in snew
    do n=1,nlevs
       call multifab_mult_mult_s_c(rho_update(n),1,0.5d0*dt,nspecies,0)
       call multifab_copy_c(rho_new(n),1,rho_old(n),1,nspecies,0)
       call multifab_plus_plus_c(rho_new(n),1,rho_update(n),1,nspecies,0)
       ! fill ghost cells for two adjacent grids including periodic boundary ghost cells
       call multifab_fill_boundary(rho_new(n))
       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(rho_new(n),1,rho_part_bc_comp,nspecies,the_bc_tower%bc_tower_array(n),dx(n,:))
    end do

    !call eos_check(mla,rho_new)
    call compute_rhotot(mla,rho_new,rhotot_new)

    call average_cc_to_face(nlevs,   rho_new,   rho_fc,1,rho_part_bc_comp,nspecies,the_bc_tower%bc_tower_array)
    call average_cc_to_face(nlevs,rhotot_new,rhotot_fc,1,    scal_bc_comp,       1,the_bc_tower%bc_tower_array)

    ! update eta and kappa here
    call compute_eta(mla,eta,eta_ed,rho_new,rhotot_new,Temp,pres,dx,the_bc_tower%bc_tower_array)
    call compute_kappa(mla,kappa)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 4 - Corrector Stochastic/Diffusive Fluxes
    ! Step 5 - Corrector Stokes Solve
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! build up rhs_v for gmres solve
    do n=1,nlevs
       do i=1,dm
          call setval(gmres_rhs_v(n,i),0.d0,all=.true.)
       end do
    end do

    ! compute grad p^*
    call compute_grad(mla,pres,gradp,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

    ! subtract grad p^* from gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradp(n,i),1,1,0)
       end do
    end do

    if (algorithm_type .eq. 2) then
       weights(:) = 1.d0/sqrt(2.d0)
    end if

    ! add div(Sigma^(2)) to gmres_rhs_v
    if (variance_coef_mom .ne. 0.d0) then
       call stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,gmres_rhs_v, &
                                 eta,eta_ed,Temp,Temp_ed,dx,dt,weights)
    end if

    ! add rho^{*,n+1}*g to gmres_rhs_v
    if (any(grav(1:dm) .ne. 0.d0)) then
       call mk_grav_force(mla,gmres_rhs_v,rhotot_fc,rhotot_fc,the_bc_tower)
    end if

    ! initialize rhs_p for gmres solve to zero
    do n=1,nlevs
       call setval(gmres_rhs_p(n),0.d0,all=.true.)
    end do

    ! reset inhomogeneous bc condition to deal with reservoirs
    call set_inhomogeneous_vel_bcs(mla,vel_bc_n,vel_bc_t,eta_ed,dx, &
                                   the_bc_tower%bc_tower_array)

    do n=1,nlevs
       do i=1,dm
          call setval(flux_total(n,i),0.d0,all=.true.)
       end do
    end do

    ! compute diffusive and stochastic mass fluxes
    ! this computes "-F" so we later multiply by -1
    call compute_mass_fluxdiv_wrapper(mla,rho_new,rhotot_new, &
                                      diff_mass_fluxdiv,stoch_mass_fluxdiv,Temp, &
                                      flux_total,dt,time,dx,weights, &
                                      the_bc_tower%bc_tower_array)

    do n=1,nlevs
       call multifab_mult_mult_s_c(diff_mass_fluxdiv(n),1,-1.d0,nspecies,0)
       if (variance_coef_mass .ne. 0.d0) then
          call multifab_mult_mult_s_c(stoch_mass_fluxdiv(n),1,-1.d0,nspecies,0)
       end if
       do i=1,dm
          call multifab_mult_mult_s_c(flux_total(n,i),1,-1.d0,nspecies,0)
       end do
    end do

    ! set the Dirichlet velocity value on reservoir faces
    call reservoir_bc_fill(mla,flux_total,vel_bc_n,the_bc_tower%bc_tower_array)

    do n=1,nlevs
       do i=1,nspecies
          call saxpy(gmres_rhs_p(n),1,-1.d0/rhobar(i), diff_mass_fluxdiv(n),i,1)
          if (variance_coef_mass .ne. 0.d0) then
             call saxpy(gmres_rhs_p(n),1,-1.d0/rhobar(i),stoch_mass_fluxdiv(n),i,1)
          end if
       end do
    end do

    ! reset rho_update for all scalars to zero
    ! then, set rho_update for rho1 to F^{*,n+1/2} = div(rho*chi grad c)^{*,n+1/2} + div(Psi^(2))
    do n=1,nlevs
       call multifab_setval(rho_update(n),0.d0,all=.true.)
       ! add fluxes
       call multifab_plus_plus_c(rho_update(n),1, diff_mass_fluxdiv(n),1,nspecies)
       if (variance_coef_mass .ne. 0.d0) then
          call multifab_plus_plus_c(rho_update(n),1,stoch_mass_fluxdiv(n),1,nspecies)
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
    call diffusive_m_fluxdiv(mla,gmres_rhs_v,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    ! compute div v^*
    call compute_div(mla,umac,divu,dx,1,1,1)

    ! add div v^* to gmres_rhs_p
    ! now gmres_rhs_p = div v^* - S^{*,n+1/2}
    ! the sign convention is correct since we solve -div(delta v) = gmres_rhs_p
    do n=1,nlevs
       call multifab_plus_plus_c(gmres_rhs_p(n),1,divu(n),1,1,0)
    end do

    ! set the initial guess to zero
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(dumac(n,i),0.d0,all=.true.)
       end do
       call multifab_setval(dp(n),0.d0,all=.true.)
    end do

    ! call gmres to compute delta v and delta p
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dp,rhotot_fc, &
               eta,eta_ed,kappa,theta_alpha)
                              
    gmres_abs_tol = gmres_abs_tol_in ! Restore the desired tolerance         

    ! compute v^{n+1/2} = v^* + delta v
    ! compute p^{n+1/2} = p^* + delta p
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(umac(n,i),1,dumac(n,i),1,1,0)
       end do
       call multifab_plus_plus_c(pres(n),1,dp(n),1,1,0)
    end do

    do n=1,nlevs
       ! presure ghost cells
       call multifab_fill_boundary(pres(n))
       call multifab_physbc(pres(n),1,pres_bc_comp,1,the_bc_tower%bc_tower_array(n),dx(n,:))
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
          call bds(mla,umac,rho_old,rho_update,bds_force,rho_fc,dx,dt,1,nspecies,the_bc_tower)
      else if (advection_type .eq. 3) then
          call bds_quad(mla,umac,rho_old,rho_update,bds_force,rho_fc,dx,dt,1,nspecies,the_bc_tower)
      end if
    else
       call mk_advective_s_fluxdiv(mla,umac,rho_fc,rho_update,dx,1,nspecies)
    end if

    ! compute s^{n+1} = s^n + dt * (A^{n+1/2} + F^{*,n+1/2})
    do n=1,nlevs
       call multifab_mult_mult_s_c(rho_update(n),1,dt,nspecies,0)
       call multifab_copy_c(rho_new(n),1,rho_old(n),1,nspecies,0)
       call multifab_plus_plus_c(rho_new(n),1,rho_update(n),1,nspecies,0)
       ! fill ghost cells for two adjacent grids including periodic boundary ghost cells
       call multifab_fill_boundary(rho_new(n))
       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(rho_new(n),1,rho_part_bc_comp,nspecies,the_bc_tower%bc_tower_array(n),dx(n,:))
    end do

    !call eos_check(mla,rho_new)
    call compute_rhotot(mla,rho_new,rhotot_new)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Compute stuff for plotfile and next time step

    ! update eta and kappa here
    call compute_eta(mla,eta,eta_ed,rho_new,rhotot_new,Temp,pres,dx,the_bc_tower%bc_tower_array)
    call compute_kappa(mla,kappa)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End Time-Advancement
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call destroy_bc_multifabs(mla)

    do n=1,nlevs
       call multifab_destroy(rho_update(n))
       call multifab_destroy(bds_force(n))
       call multifab_destroy(gmres_rhs_p(n))
       call multifab_destroy(dp(n))
       call multifab_destroy(divu(n))
       do i=1,dm
          call multifab_destroy(gmres_rhs_v(n,i))
          call multifab_destroy(dumac(n,i))
          call multifab_destroy(gradp(n,i))
          call multifab_destroy(rho_fc(n,i))
          call multifab_destroy(rhotot_fc(n,i))
          call multifab_destroy(flux_total(n,i))
       end do
    end do

  end subroutine advance_timestep_inertial

  subroutine build_bc_multifabs(mla)

    type(ml_layout), intent(in   ) :: mla

    integer :: dm,i,n,nlevs
    logical :: nodal_temp(3)

    dm = mla%dim
    nlevs = mla%nlevel

    allocate(vel_bc_n(nlevs,dm))
    if (dm .eq. 2) then
       allocate(vel_bc_t(nlevs,2))
    else if (dm .eq. 3) then
       allocate(vel_bc_t(nlevs,6))
    end if

    do n=1,nlevs
       ! boundary conditions
       do i=1,dm
          call multifab_build_edge(vel_bc_n(n,i),mla%la(n),1,0,i)
       end do
       if (dm .eq. 2) then
          ! y-velocity bc on x-faces (nodal)
          call multifab_build_nodal(vel_bc_t(n,1),mla%la(n),1,0)
          ! x-velocity bc on y-faces (nodal)
          call multifab_build_nodal(vel_bc_t(n,2),mla%la(n),1,0)
       else
          ! y-velocity bc on x-faces (nodal in y and x)
          nodal_temp(1) = .true.
          nodal_temp(2) = .true.
          nodal_temp(3) = .false.
          call multifab_build(vel_bc_t(n,1),mla%la(n),1,0,nodal_temp)
          ! z-velocity bc on x-faces (nodal in z and x)
          nodal_temp(1) = .true.
          nodal_temp(2) = .false.
          nodal_temp(3) = .true.
          call multifab_build(vel_bc_t(n,2),mla%la(n),1,0,nodal_temp)
          ! x-velocity bc on y-faces (nodal in x and y)
          nodal_temp(1) = .true.
          nodal_temp(2) = .true.
          nodal_temp(3) = .false.
          call multifab_build(vel_bc_t(n,3),mla%la(n),1,0,nodal_temp)
          ! z-velocity bc on y-faces (nodal in z and y)
          nodal_temp(1) = .false.
          nodal_temp(2) = .true.
          nodal_temp(3) = .true.
          call multifab_build(vel_bc_t(n,4),mla%la(n),1,0,nodal_temp)
          ! x-velocity bc on z-faces (nodal in x and z)
          nodal_temp(1) = .true.
          nodal_temp(2) = .false.
          nodal_temp(3) = .true.
          call multifab_build(vel_bc_t(n,5),mla%la(n),1,0,nodal_temp)
          ! y-velocity bc on z-faces (nodal in y and z)
          nodal_temp(1) = .false.
          nodal_temp(2) = .true.
          nodal_temp(3) = .true.
          call multifab_build(vel_bc_t(n,6),mla%la(n),1,0,nodal_temp)
       end if

       do i=1,dm
          call multifab_setval(vel_bc_n(n,i),0.d0,all=.true.)
       end do
       do i=1,size(vel_bc_t,dim=2)
          call multifab_setval(vel_bc_t(n,i),0.d0,all=.true.)
       end do

    end do

  end subroutine build_bc_multifabs

  subroutine destroy_bc_multifabs(mla)

    type(ml_layout), intent(in   ) :: mla

    integer :: dm,i,n,nlevs

    dm = mla%dim
    nlevs = mla%nlevel

    do n=1,nlevs
       do i=1,dm          
          call multifab_destroy(vel_bc_n(n,i))
       end do
       do i=1,size(vel_bc_t,dim=2)
          call multifab_destroy(vel_bc_t(n,i))
       end do
    end do

    deallocate(vel_bc_n,vel_bc_t)

  end subroutine destroy_bc_multifabs

end module advance_timestep_inertial_module
