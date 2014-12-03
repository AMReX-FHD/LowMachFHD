module advance_timestep_overdamped_module

  use ml_layout_module
  use multifab_module
  use define_bc_module
  use convert_stag_module
  use convert_variables_module
  use convert_m_to_umac_module
  use mk_advective_s_fluxdiv_module
  use mk_advective_m_fluxdiv_module
  use diffusive_rhoc_fluxdiv_module
  use diffusive_m_fluxdiv_module
  use mk_external_force_module
  use mk_grav_force_module
  use stochastic_m_fluxdiv_module
  use stochastic_rhoc_fluxdiv_module
  use mk_baro_fluxdiv_module
  use bds_module
  use gmres_module
  use init_module
  use div_and_grad_module
  use bc_module
  use multifab_physbc_module
  use multifab_physbc_extrap_module
  use multifab_physbc_stag_module
  use fill_rho_ghost_cells_module
  use probin_common_module, only: advection_type, grav, rhobar, algorithm_type, &
                                  barodiffusion_type
  use probin_gmres_module, only: gmres_abs_tol, gmres_rel_tol

  implicit none

  private

  public :: advance_timestep_overdamped

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

  subroutine advance_timestep_overdamped(mla,mnew,umac,sold,snew,s_fc,prim,pres, &
                                         chi,chi_fc,eta,eta_ed,kappa,gradp_baro, &
                                         dx,dt,time,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: mnew(:,:) ! only a diagnostic for plotfile purposes
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    ! s_fc and prim need to enter consistent with sold and leave consistent with snew
    type(multifab) , intent(inout) :: s_fc(:,:)
    type(multifab) , intent(inout) :: prim(:)
    type(multifab) , intent(inout) :: pres(:)
    ! chi, eta, and kappa need to enter consistent with sold and leave consistent with snew
    type(multifab) , intent(inout) :: chi(:)
    type(multifab) , intent(inout) :: chi_fc(:,:)
    type(multifab) , intent(inout) :: eta(:)
    type(multifab) , intent(inout) :: eta_ed(:,:) ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(inout) :: kappa(:)
    type(multifab) , intent(inout) :: gradp_baro(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt,time
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    type(multifab) ::    s_update(mla%nlevel)
    type(multifab) ::   bds_force(mla%nlevel)
    type(multifab) :: gmres_rhs_p(mla%nlevel)
    type(multifab) ::          dp(mla%nlevel)
    type(multifab) ::        divu(mla%nlevel)

    type(multifab) ::  gmres_rhs_v(mla%nlevel,mla%dim)
    type(multifab) ::        dumac(mla%nlevel,mla%dim)
    type(multifab) ::        gradp(mla%nlevel,mla%dim)
    type(multifab) ::     s_nd_old(mla%nlevel)
    type(multifab) ::        s_tmp(mla%nlevel)

    integer :: i,dm,n,nlevs

    real(kind=dp_t) :: S_fac, theta_alpha, norm_pre_rhs

    real(kind=dp_t) :: weights(algorithm_type)

    weights(:) = 0.d0
    weights(1) = 1.d0

    nlevs = mla%nlevel
    dm = mla%dim

    theta_alpha = 0.d0

    S_fac = (1.d0/rhobar(1) - 1.d0/rhobar(2))
    
    call build_bc_multifabs(mla)
    
    do n=1,nlevs
       call multifab_build(   s_update(n),mla%la(n),2    ,0)
       call multifab_build(  bds_force(n),mla%la(n),2    ,1)
       call multifab_build(gmres_rhs_p(n),mla%la(n),1    ,0)
       call multifab_build(         dp(n),mla%la(n),1    ,1)
       call multifab_build(       divu(n),mla%la(n),1    ,0)
       do i=1,dm
          call multifab_build_edge(    gmres_rhs_v(n,i),mla%la(n),1    ,0,i)
          call multifab_build_edge(          dumac(n,i),mla%la(n),1    ,1,i)
          call multifab_build_edge(          gradp(n,i),mla%la(n),1    ,0,i)
       end do
    end do

    do n=1,nlevs
       call setval(s_update(n),0.d0,all=.true.)
       call setval(bds_force(n),0.d0,all=.true.)
       do i=1,dm
          call setval(dumac(n,i),0.d0,all=.true.)
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 1 - Predictor Stochastic/Diffusive Fluxes
    ! Step 2 - Predictor Stokes Solve
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! fill the stochastic multifabs with a new set of random numbers
    ! if this is the first step we already have random numbers from the initial projection
    if (time .ne. 0.d0) then
       call fill_m_stochastic(mla)
       call fill_rhoc_stochastic(mla)
    end if

    ! build up rhs_v for gmres solve
    do n=1,nlevs
       do i=1,dm
          call setval(gmres_rhs_v(n,i),0.d0,all=.true.)
       end do
    end do

    ! compute grad p^{n-1/2}
    call compute_grad(mla,pres,gradp,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

    ! barodiffusion uses lagged pressure
    if (barodiffusion_type .eq. 2) then
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(gradp_baro(n,i),1,gradp(n,i),1,1,0)
          end do
       end do
    end if

    ! subtract grad p^{n-1/2} from gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradp(n,i),1,1,0)
       end do
    end do

    ! add div(Sigma^(1)) to gmres_rhs_v
    if (algorithm_type .eq. 1) then
       call stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,gmres_rhs_v, &
                                 eta,eta_ed,dx,dt,weights)
    else if (algorithm_type .eq. 2) then
       call stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,gmres_rhs_v, &
                                 eta,eta_ed,dx,0.5d0*dt,weights)
    end if

    ! add rho^n*g to gmres_rhs_v
    if (any(grav(1:dm) .ne. 0.d0)) then
       call mk_grav_force(mla,gmres_rhs_v,s_fc,s_fc,the_bc_tower)
    end if

    ! add external forcing to gmres_rhs_v
    call mk_external_m_force(mla,gmres_rhs_v,dx,time+0.5d0*dt,the_bc_tower)

    ! initialize rhs_p for gmres solve to zero
    do n=1,nlevs
       call setval(gmres_rhs_p(n),0.d0,all=.true.)
    end do

    ! reset inhomogeneous bc condition to deal with reservoirs
    call set_inhomogeneous_vel_bcs(mla,vel_bc_n,vel_bc_t,eta_ed,dx,time, &
                                   the_bc_tower%bc_tower_array)

    ! add div(rho*chi grad c)^n to rhs_p
    call diffusive_rhoc_fluxdiv(mla,gmres_rhs_p,1,prim,s_fc,chi_fc,dx, &
                                the_bc_tower%bc_tower_array,vel_bc_n)

    if (barodiffusion_type .gt. 0) then
       ! add baro-diffusion flux diveregnce to rhs_p
       call mk_baro_fluxdiv(mla,gmres_rhs_p,1,s_fc,chi_fc,gradp_baro,dx, &
                            the_bc_tower%bc_tower_array,vel_bc_n)
    end if

    ! add div(Psi^(1)) to rhs_p
    if (algorithm_type .eq. 1) then
       call stochastic_rhoc_fluxdiv(mla,the_bc_tower%bc_tower_array,gmres_rhs_p,s_fc, &
                                    chi_fc,dx,dt,vel_bc_n,weights)
    else if (algorithm_type .eq. 2) then
       call stochastic_rhoc_fluxdiv(mla,the_bc_tower%bc_tower_array,gmres_rhs_p,s_fc, &
                                    chi_fc,dx,0.5d0*dt,vel_bc_n,weights)  
    end if

    ! add external forcing for rho*c
    call mk_external_s_force(mla,gmres_rhs_p,dx,time,1)

    ! s_update is zero
    ! set s_update for rho1 to F^n = div(rho*chi grad c)^n + div(Psi^(1))
    do n=1,nlevs
       call multifab_copy_c(s_update(n),2,gmres_rhs_p(n),1,1,0)
    end do

    ! multiply gmres_rhs_p by -S_fac, so gmres_rhs_p = -S^n
    do n=1,nlevs
       call multifab_mult_mult_s_c(gmres_rhs_p(n),1,-S_fac,1,0)
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
          call multifab_setval(   dp(n)  ,0.d0,all=.true.)
       end do
    end do

    gmres_abs_tol = 0.d0

    ! call gmres to compute delta v and delta p
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dp,s_fc, &
               eta,eta_ed,kappa,theta_alpha,norm_pre_rhs)

    ! for the corrector gmres solve we want the stopping criteria based on the
    ! norm of the preconditioned rhs from the predictor gmres solve.  otherwise
    ! for cases where du in the corrector should be small the gmres stalls
    gmres_abs_tol = norm_pre_rhs*gmres_rel_tol

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
       call multifab_physbc(pres(n),1,pres_bc_comp,1,the_bc_tower%bc_tower_array(n), &
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

    ! add A^n for scalars to s_update
    if (advection_type .ge. 1) then
       do n=1,nlevs
          call multifab_copy_c(bds_force(n),1,s_update(n),1,2,0)
          call multifab_fill_boundary(bds_force(n))
       end do

       if (advection_type .eq. 1 .or. advection_type .eq. 2) then

          ! s_fc (computed above) and s_nd_old (computed here) are used to set boundary conditions
          do n=1,nlevs
             call multifab_build_nodal(s_nd_old(n),mla%la(n),2,1)
          end do
          call average_cc_to_node(nlevs,sold,s_nd_old,1,scal_bc_comp,2,the_bc_tower%bc_tower_array)

          ! the input s_tmp needs to have ghost cells filled with multifab_physbc_extrap
          ! instead of multifab_physbc
          do n=1,nlevs
             call multifab_build(s_tmp(n),mla%la(n),2,sold(n)%ng)
             call multifab_copy(s_tmp(n),sold(n),s_tmp(n)%ng)
             call multifab_physbc_extrap(s_tmp(n),1,scal_bc_comp,2, &
                                         the_bc_tower%bc_tower_array(n))
          end do

          call bds(mla,umac,s_tmp,s_update,bds_force,s_fc,s_nd_old,dx,0.5d0*dt,1,2,scal_bc_comp, &
                   the_bc_tower,proj_type_in=1)

       else if (advection_type .eq. 3 .or. advection_type .eq. 4) then

          call bds_quad(mla,umac,sold,s_update,bds_force,s_fc,dx,0.5d0*dt,1,2,scal_bc_comp, &
                        the_bc_tower,proj_type_in=1)

       end if
    else
       call mk_advective_s_fluxdiv(mla,umac,s_fc,s_update,dx,1,2)
    end if

    ! compute s^{*,n+1/2} = s^n + (dt/2) * (A^n + F^n)
    ! store result in snew
    do n=1,nlevs
       call multifab_mult_mult_s_c(s_update(n),1,0.5d0*dt,2,0)
       call multifab_copy_c(snew(n),1,sold(n),1,2,0)
       call multifab_plus_plus_c(snew(n),1,s_update(n),1,2,0)
    end do

    ! convert s^{*,n+1/2} to prim
    call convert_cons_to_prim(mla,snew,prim,.true.)

    ! fill ghost cells for prim
    do n=1,nlevs
       call multifab_fill_boundary(prim(n))
       call multifab_physbc(prim(n),2,scal_bc_comp+1,1,the_bc_tower%bc_tower_array(n), &
                            dx_in=dx(n,:))
       call fill_rho_ghost_cells(prim(n),the_bc_tower%bc_tower_array(n))
    end do

    ! convert prim to s^{*,n+1/2} in valid and ghost region
    ! now s^{*,n+1/2} properly filled ghost cells
    call convert_cons_to_prim(mla,snew,prim,.false.)

    ! average s^{*,n+1/2} to faces
    call average_cc_to_face(nlevs,snew,s_fc,1,scal_bc_comp,2,the_bc_tower%bc_tower_array)

    ! compute (chi,eta,kappa)^{*,n+1/2}
    call compute_chi(mla,chi,chi_fc,prim,dx,the_bc_tower%bc_tower_array)
    call compute_eta(mla,eta,eta_ed,prim,dx,the_bc_tower%bc_tower_array)
    call compute_kappa(mla,kappa,prim,dx)

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

    ! barodiffusion
    if (barodiffusion_type .eq. 2) then
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(gradp_baro(n,i),1,gradp(n,i),1,1,0)
          end do
       end do
    end if

    ! subtract grad p^* from gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradp(n,i),1,1,0)
       end do
    end do

    if (algorithm_type .eq. 2) then
       weights = 1.d0/sqrt(2.d0)
    end if

    ! add div(Sigma^(2)) to gmres_rhs_v
    call stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,gmres_rhs_v, &
                              eta,eta_ed,dx,dt,weights)

    ! add rho^{*,n+1/2}*g gmres_rhs_v
    if (any(grav(1:dm) .ne. 0.d0)) then
       call mk_grav_force(mla,gmres_rhs_v,s_fc,s_fc,the_bc_tower)
    end if

    ! add external forcing to gmres_rhs_v
    call mk_external_m_force(mla,gmres_rhs_v,dx,time+0.5d0*dt,the_bc_tower)

    ! reset inhomogeneous bc condition to deal with reservoirs
    call set_inhomogeneous_vel_bcs(mla,vel_bc_n,vel_bc_t,eta_ed,dx,time+0.5d0*dt, &
                                   the_bc_tower%bc_tower_array)

    ! initialize rhs_p for gmres solve to zero
    do n=1,nlevs
       call setval(gmres_rhs_p(n),0.d0,all=.true.)
    end do

    ! add div(rho*chi grad c)^{*,n+1/2} to rhs_p
    call diffusive_rhoc_fluxdiv(mla,gmres_rhs_p,1,prim,s_fc,chi_fc,dx, &
                                the_bc_tower%bc_tower_array,vel_bc_n)

    if (barodiffusion_type .gt. 0) then
       ! compute baro-diffusion flux divergence
       call mk_baro_fluxdiv(mla,gmres_rhs_p,1,s_fc,chi_fc,gradp_baro,dx, &
                            the_bc_tower%bc_tower_array,vel_bc_n)
    end if

    ! add div(Psi^(2)) to rhs_p
    call stochastic_rhoc_fluxdiv(mla,the_bc_tower%bc_tower_array,gmres_rhs_p,s_fc, &
                                 chi_fc,dx,dt,vel_bc_n,weights)

    ! add external forcing for rho*c
    call mk_external_s_force(mla,gmres_rhs_p,dx,time+0.5d0*dt,1)

    ! set s_update for rho to zero
    ! set s_update for rho1 to F^{*,n+1/2} = div(rho*chi grad c)^{*,n+1/2} + div(Psi^(2))
    do n=1,nlevs
       call multifab_setval_c(s_update(n),0.d0,1,1,all=.true.)
       call multifab_copy_c(s_update(n),2,gmres_rhs_p(n),1,1,0)
    end do

    ! multiply gmres_rhs_p by -S_fac, so gmres_rhs_p = -S^{*,n+1/2}
    do n=1,nlevs
       call multifab_mult_mult_s_c(gmres_rhs_p(n),1,-S_fac,1,0)
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
          call multifab_setval(   dp(n)  ,0.d0,all=.true.)
       end do
    end do

    ! call gmres to compute delta v and delta p
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dp,s_fc, &
               eta,eta_ed,kappa,theta_alpha)

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
       call multifab_physbc(pres(n),1,pres_bc_comp,1,the_bc_tower%bc_tower_array(n), &
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

    ! add A^{n+1/2} for scalars to s_update
    if (advection_type .ge. 1) then
       do n=1,nlevs
          call multifab_copy_c(bds_force(n),1,s_update(n),1,2,0)
          call multifab_fill_boundary(bds_force(n))
       end do
       if (advection_type .eq. 1 .or. advection_type .eq. 2) then

          call bds(mla,umac,s_tmp,s_update,bds_force,s_fc,s_nd_old,dx,dt,1,2,scal_bc_comp, &
                   the_bc_tower,proj_type_in=1)

          do n=1,nlevs
             call multifab_destroy(s_nd_old(n))
             call multifab_destroy(s_tmp(n))
          end do

       else if (advection_type .eq. 3 .or. advection_type .eq. 4) then
          call bds_quad(mla,umac,sold,s_update,bds_force,s_fc,dx,dt,1,2,scal_bc_comp, &
                        the_bc_tower,proj_type_in=1)
       end if
    else
       call mk_advective_s_fluxdiv(mla,umac,s_fc,s_update,dx,1,2)
    end if

    ! compute s^{n+1} = s^n + dt * (A^{n+1/2} + F^{*,n+1/2})
    do n=1,nlevs
       call multifab_mult_mult_s_c(s_update(n),1,dt,2,0)
       call multifab_copy_c(snew(n),1,sold(n),1,2,0)
       call multifab_plus_plus_c(snew(n),1,s_update(n),1,2,0)
    end do

    ! convert s^{n+1} to prim
    call convert_cons_to_prim(mla,snew,prim,.true.)

    ! fill ghost cells for prim
    do n=1,nlevs
       call multifab_fill_boundary(prim(n))
       call multifab_physbc(prim(n),2,scal_bc_comp+1,1,the_bc_tower%bc_tower_array(n), &
                            dx_in=dx(n,:))
       call fill_rho_ghost_cells(prim(n),the_bc_tower%bc_tower_array(n))
    end do

    ! convert prim to s^{n+1} in valid and ghost region
    ! now s^{n+1} properly filled ghost cells
    call convert_cons_to_prim(mla,snew,prim,.false.)

    ! compute s^{n+1} to faces
    call average_cc_to_face(nlevs,snew,s_fc,1,scal_bc_comp,2,the_bc_tower%bc_tower_array)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Compute stuff for plotfile and next time step
    
    ! mnew actually holds rho^{n+1} v^{*,n+1/2}, not (1/2)(rho^n + rho^{n+1}) v^{*,n+1/2}
    ! mnew is just a diagnostic, it does not enter the algorithm so for now this is fine
    call convert_m_to_umac(mla,s_fc,mnew,umac,.false.)

    ! compute (chi,eta,kappa)^{n+1}
    call compute_chi(mla,chi,chi_fc,prim,dx,the_bc_tower%bc_tower_array)
    call compute_eta(mla,eta,eta_ed,prim,dx,the_bc_tower%bc_tower_array)
    call compute_kappa(mla,kappa,prim,dx)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End Time-Advancement
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call destroy_bc_multifabs(mla)

    do n=1,nlevs
       call multifab_destroy(s_update(n))
       call multifab_destroy(bds_force(n))
       call multifab_destroy(gmres_rhs_p(n))
       call multifab_destroy(dp(n))
       call multifab_destroy(divu(n))
       do i=1,dm
          call multifab_destroy(gmres_rhs_v(n,i))
          call multifab_destroy(dumac(n,i))
          call multifab_destroy(gradp(n,i))
       end do
    end do

  end subroutine advance_timestep_overdamped

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

end module advance_timestep_overdamped_module
