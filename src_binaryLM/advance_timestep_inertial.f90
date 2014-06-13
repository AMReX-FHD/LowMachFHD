module advance_timestep_inertial_module

  use ml_layout_module
  use multifab_module
  use define_bc_module
  use convert_stag_module
  use convert_variables_module
  use convert_m_to_umac_module
  use convert_to_homogeneous_module
  use mk_advective_s_fluxdiv_module
  use mk_advective_m_fluxdiv_module
  use mk_baro_fluxdiv_module
  use diffusive_m_fluxdiv_module
  use diffusive_rhoc_fluxdiv_module
  use mk_grav_force_module
  use mk_external_force_module
  use stochastic_m_fluxdiv_module
  use stochastic_rhoc_fluxdiv_module
  use gmres_module
  use init_module
  use div_and_grad_module
  use bds_module
  use bc_module
  use multifab_physbc_module
  use multifab_physbc_stag_module
  use zero_edgeval_module
  use fill_rho_ghost_cells_module
  use probin_binarylm_module, only: barodiffusion_type
  use probin_common_module, only: advection_type, grav, rhobar
  use probin_gmres_module, only: gmres_abs_tol, gmres_rel_tol

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

  subroutine advance_timestep_inertial(mla,mold,mnew,umac,sold,snew,s_fc,prim,pres,chi,chi_fc, &
                              eta,eta_ed,kappa,rhoc_fluxdiv, &
                              gp_fc,dx,dt,time,the_bc_tower)

    ! rhoc_fluxdiv enters with t^n diffusive, baro-diffusion, and stochastic rho*c
    ! mass fluxes and exits with t^n+1 fluxes

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: mold(:,:)
    type(multifab) , intent(inout) :: mnew(:,:)
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(inout) :: s_fc(:,:)
    type(multifab) , intent(inout) :: prim(:)
    type(multifab) , intent(inout) :: pres(:)
    type(multifab) , intent(inout) :: chi(:)
    type(multifab) , intent(inout) :: chi_fc(:,:)
    type(multifab) , intent(inout) :: eta(:)
    type(multifab) , intent(inout) :: eta_ed(:,:) ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(inout) :: kappa(:)
    type(multifab) , intent(inout) :: rhoc_fluxdiv(:)
    type(multifab) , intent(inout) :: gp_fc(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt,time
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    type(multifab) ::    s_update(mla%nlevel)
    type(multifab) ::   bds_force(mla%nlevel)
    type(multifab) :: gmres_rhs_p(mla%nlevel)
    type(multifab) ::          dp(mla%nlevel)
    type(multifab) ::        divu(mla%nlevel)

    type(multifab) ::        mtemp(mla%nlevel,mla%dim)
    type(multifab) ::  gmres_rhs_v(mla%nlevel,mla%dim)
    type(multifab) ::  m_a_fluxdiv(mla%nlevel,mla%dim)
    type(multifab) ::  m_d_fluxdiv(mla%nlevel,mla%dim)
    type(multifab) ::  m_s_fluxdiv(mla%nlevel,mla%dim)
    type(multifab) ::        dumac(mla%nlevel,mla%dim)
    type(multifab) ::     umac_tmp(mla%nlevel,mla%dim)
    type(multifab) ::        gradp(mla%nlevel,mla%dim)
    type(multifab) ::     s_fc_old(mla%nlevel,mla%dim)
    type(multifab) ::   chi_fc_old(mla%nlevel,mla%dim)

    integer :: i,dm,n,nlevs

    real(kind=dp_t) :: S_fac, theta_alpha, norm_pre_rhs

    real(kind=dp_t), allocatable :: weights(:)

    allocate(weights(1))
    weights(1) = 1.d0

    nlevs = mla%nlevel
    dm = mla%dim

    theta_alpha = 1.d0/dt

    S_fac = (1.d0/rhobar(1) - 1.d0/rhobar(2))
    
    call build_bc_multifabs(mla)

    do n=1,nlevs
       call multifab_build(   s_update(n),mla%la(n),2,0)
       call multifab_build(  bds_force(n),mla%la(n),2,1)
       call multifab_build(gmres_rhs_p(n),mla%la(n),1,0)
       call multifab_build(         dp(n),mla%la(n),1,1)
       call multifab_build(       divu(n),mla%la(n),1,0)
       do i=1,dm
          call multifab_build_edge(         mtemp(n,i),mla%la(n),1,0,i)
          call multifab_build_edge(   gmres_rhs_v(n,i),mla%la(n),1,0,i)
          call multifab_build_edge(   m_a_fluxdiv(n,i),mla%la(n),1,0,i)
          call multifab_build_edge(   m_d_fluxdiv(n,i),mla%la(n),1,0,i)
          call multifab_build_edge(   m_s_fluxdiv(n,i),mla%la(n),1,0,i)
          call multifab_build_edge(         dumac(n,i),mla%la(n),1,1,i)
          call multifab_build_edge(      umac_tmp(n,i),mla%la(n),1,1,i)
          call multifab_build_edge(         gradp(n,i),mla%la(n),1,0,i)
          call multifab_build_edge(      s_fc_old(n,i),mla%la(n),2,1,i)
          call multifab_build_edge(    chi_fc_old(n,i),mla%la(n),1,0,i)
       end do

    end do

    do n=1,nlevs
       call setval(s_update(n),0.d0,all=.true.)
       call setval(bds_force(n),0.d0,all=.true.)
       do i=1,dm
          call setval(m_a_fluxdiv(n,i),0.d0,all=.true.)
          call setval(m_d_fluxdiv(n,i),0.d0,all=.true.)
          call setval(m_s_fluxdiv(n,i),0.d0,all=.true.)
          call setval(      dumac(n,i),0.d0,all=.true.)
       end do
    end do

    ! make copies of old quantities
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(  umac_tmp(n,i),1,  umac(n,i),1,1    ,1)
          call multifab_copy_c(  s_fc_old(n,i),1,  s_fc(n,i),1,2    ,1)
          call multifab_copy_c(chi_fc_old(n,i),1,chi_fc(n,i),1,1    ,0)
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 1 - Calculate Predictor Diffusive and Stochastic Fluxes
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! rhoc_fluxdiv contains the diffusive, stochastic, and baro-diffusion fluxes
    ! this was already done in Step 0 (initialization) or Step 6 from the previous time step

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 2 - Predictor Euler Step
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! add D^n  for rho1 to s_update
    ! add St^n for rho1 to s_update
    ! add baro-diffusion^n to s_update
    do n=1,nlevs
       call multifab_plus_plus_c(s_update(n),2,rhoc_fluxdiv(n),1,1,0)
    end do

    ! add external forcing for rho*c
    call mk_external_s_force(mla,s_update,dx,time,2)

    ! add A^n for scalars to s_update
    if (advection_type .ge. 1) then

       do n=1,nlevs
          call multifab_copy_c(bds_force(n),1,s_update(n),1,2,0)
          call multifab_fill_boundary(bds_force(n))
       end do

       if (advection_type .eq. 1 .or. advection_type .eq. 2) then
          call bds(mla,umac,sold,s_update,bds_force,s_fc,dx,dt,1,2,the_bc_tower)
       else if (advection_type .eq. 3) then
          call bds_quad(mla,umac,sold,s_update,bds_force,s_fc,dx,dt,1,2,the_bc_tower)
       end if

    else

       call mk_advective_s_fluxdiv(mla,umac,s_fc,s_update,dx,1,2)

    end if

    ! set snew = s^{*,n+1} = s^n + dt * (A^n + D^n + St^n)
    do n=1,nlevs
       call multifab_mult_mult_s_c(s_update(n),1,dt,2,0)
       call multifab_copy_c(snew(n),1,sold(n),1,2,0)
       call multifab_plus_plus_c(snew(n),1,s_update(n),1,2,0)
    end do

    ! compute prim^{*,n+1} from s^{*,n+1} in valid region
    call convert_cons_to_prim(mla,snew,prim,.true.)

    ! fill ghost cells for prim^{*,n+1}
    do n=1,nlevs
       call multifab_fill_boundary(prim(n))
       call multifab_physbc(prim(n),2,scal_bc_comp+1,1,the_bc_tower%bc_tower_array(n),dx(n,:))
       call fill_rho_ghost_cells(prim(n),the_bc_tower%bc_tower_array(n))
    end do

    ! convert prim^{*,n+1} to s^{*,n+1} in valid and ghost region
    ! now snew = s^{*,n+1} has properly filled ghost cells
    call convert_cons_to_prim(mla,snew,prim,.false.)

    ! compute s^{*,n+1} to faces
    call average_cc_to_face(nlevs,snew,s_fc,1,scal_bc_comp,2,the_bc_tower%bc_tower_array)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 3 - Calculate Corrector Diffusive and Stochastic Fluxes
    ! Step 4 - Predictor Crank-Nicolson Step
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! build up rhs_v for gmres solve: first set gmres_rhs_v to mold/dt
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(gmres_rhs_v(n,i),1,mold(n,i),1,1,0)
          call multifab_mult_mult_s_c(gmres_rhs_v(n,i),1,1.d0/dt,1,0)
       end do
    end do

    ! compute grad p^n
    call compute_grad(mla,pres,gradp,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

    ! subtract grad p^n from gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradp(n,i),1,1,0)
       end do
    end do

    ! compute m_a_fluxdiv = A^n for momentum
    call mk_advective_m_fluxdiv(mla,umac,mold,m_a_fluxdiv,dx, &
                                the_bc_tower%bc_tower_array)

    ! add A^n for momentum to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_a_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! compute m_d_fluxdiv = A_0^n v^n
    call diffusive_m_fluxdiv(mla,m_d_fluxdiv,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    ! add (1/2) A_0^n v^n to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_s_c(m_d_fluxdiv(n,i),1,0.5d0,1,0)
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_d_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! compute m_s_fluxdiv = div(Sigma^n)
    call stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,m_s_fluxdiv,eta,eta_ed,dx,dt,weights)

    ! add div(Sigma^n) to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_s_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! add gravity term
    if (any(grav(1:dm) .ne. 0.d0)) then
       call mk_grav_force(mla,gmres_rhs_v,s_fc_old,s_fc_old,the_bc_tower)
    end if

    ! initialize rhs_p for gmres solve to zero
    do n=1,nlevs
       call setval(gmres_rhs_p(n),0.d0,all=.true.)
    end do

    ! compute (chi,eta,kappa)^{*,n+1}
    call compute_chi(mla,chi,chi_fc,prim,dx,the_bc_tower%bc_tower_array)
    call compute_eta(mla,eta,eta_ed,prim,dx,the_bc_tower%bc_tower_array)
    call compute_kappa(mla,kappa,prim,dx)

    ! reset inhomogeneous bc condition to deal with reservoirs
    call set_inhomogeneous_vel_bcs(mla,vel_bc_n,vel_bc_t,eta_ed,dx, &
                                   the_bc_tower%bc_tower_array)

    ! add div(rho*chi grad c)^{*,n+1} to rhs_p
    call diffusive_rhoc_fluxdiv(mla,gmres_rhs_p,1,prim,s_fc,chi_fc,dx, &
                                the_bc_tower%bc_tower_array,vel_bc_n)

    if (barodiffusion_type .gt. 0) then
       ! add baro-diffusion flux diveregnce to rhs_p
       call mk_baro_fluxdiv(mla,gmres_rhs_p,1,s_fc,chi_fc,gp_fc,dx, &
                            the_bc_tower%bc_tower_array,vel_bc_n)
    end if

    ! add div(Psi^n') to rhs_p
    call stochastic_rhoc_fluxdiv(mla,the_bc_tower%bc_tower_array,gmres_rhs_p,s_fc, &
                                 chi_fc,dx,dt,vel_bc_n,weights)

    ! add external forcing for rho*c
    call mk_external_s_force(mla,gmres_rhs_p,dx,time+dt,1)

    ! modify umac to respect the boundary conditions we want after the next gmres solve
    ! thus when we add (1/2) A_0^n vbar^n and -rho^{*,n+1}*vbar^n/dt to gmres_rhs_v 
    ! and div vbar^n to gmres_rhs_p we
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
   call convert_m_to_umac(mla,s_fc,mtemp,umac,.false.)

   do n=1,nlevs
      do i=1,dm

         ! multiply mtemp by 1/dt
         call multifab_mult_mult_s_c(mtemp(n,i),1,1.d0/dt,1,0)

         ! subtract rho^{*,n+1} * vbar^n / dt from gmres_rhs_v
         call multifab_sub_sub_c(gmres_rhs_v(n,i),1,mtemp(n,i),1,1,0)

      end do
   end do

   do n=1,nlevs
      do i=1,dm
         ! reset mtemp
         call multifab_setval(mtemp(n,i),0.d0,all=.true.)
      end do
   end do

    ! compute mtemp = A_0^n vbar^n
    call diffusive_m_fluxdiv(mla,mtemp,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    ! add (1/2) A_0^n vbar^n to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_s_c(mtemp(n,i),1,0.5d0,1,0)
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,mtemp(n,i),1,1,0)
       end do
    end do

    ! reset s_update for all scalars to zero
    ! then, set s_update for rho1 to F^{*,n+1} = div(rho*chi grad c)^{*,n+1} + div(Psi^n)
    ! it is used in Step 5 below
    do n=1,nlevs
       call multifab_setval_c(s_update(n),0.d0,1,1,all=.true.)
       call multifab_copy_c(s_update(n),2,gmres_rhs_p(n),1,1,0)
    end do

    ! multiply gmres_rhs_p by -S_fac
    ! now gmres_rhs_p = -S^{*,n+1}
    do n=1,nlevs
       call multifab_mult_mult_s_c(gmres_rhs_p(n),1,-S_fac,1,0)
    end do

    ! compute div(vbar^n)
    call compute_div(mla,umac,divu,dx,1,1,1)

    ! add div(vbar^n) to gmres_rhs_p
    ! now gmres_rhs_p = div(vbar^n) - S^{*,n+1}
    ! the sign convention is correct since we solve -div(delta v) = gmres_rhs_p
    do n=1,nlevs
       call multifab_plus_plus_c(gmres_rhs_p(n),1,divu(n),1,1,0)
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
       call multifab_setval(dp(n),0.d0,all=.true.)
    end do

    do n=1,nlevs
       call zero_edgeval_physical(gmres_rhs_v(n,:),1,1,the_bc_tower%bc_tower_array(n))
    end do

    gmres_abs_tol_in = gmres_abs_tol ! Save this 

    ! call gmres to compute delta v and delta p
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dp,s_fc, &
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

    ! compute v^{*,n+1} = v^n + delta v
    ! compute p^{*,n+1}= p^n + delta p
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

    if (barodiffusion_type .eq. 2) then
       ! compute grad p^{n+1,*}
       call compute_grad(mla,pres,gp_fc,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)
    end if

    ! convert v^{*,n+1} to rho^{*,n+1}v^{*,n+1} in valid and ghost region
    ! now mnew has properly filled ghost cells
    call convert_m_to_umac(mla,s_fc,mnew,umac,.false.)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 5 - Trapezoidal Scalar Corrector
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! s_update already contains D^{*,n+1} + St^{*,n+1} for rho1 from above
    ! add A^{*,n+1} for s to s_update
    if (advection_type .ge. 1) then

       do n=1,nlevs
          ! bds force currently contains D^n + St^n
          ! add D^{*,n+1} + St^{*,n+1} and then multiply by 1/2
          call multifab_plus_plus_c(bds_force(n),1,s_update(n),1,2,0)
          call multifab_mult_mult_s_c(bds_force(n),1,0.5d0,2,0)
          call multifab_fill_boundary(bds_force(n))
          do i=1,dm
             call multifab_plus_plus_c(umac_tmp(n,i),1,umac(n,i),1,1,1)
             call multifab_mult_mult_s_c(umac_tmp(n,i),1,0.5d0,1,1)
          end do
          call setval(s_update(n),0.d0,all=.true.)
       end do

       if (advection_type .eq. 1 .or. advection_type .eq. 2) then
          call bds(mla,umac_tmp,sold,s_update,bds_force,s_fc,dx,dt,1,2,the_bc_tower)
       else if (advection_type .eq. 3) then
          call bds_quad(mla,umac_tmp,sold,s_update,bds_force,s_fc,dx,dt,1,2,the_bc_tower)
       end if    

       ! snew = s^n + dt * A^{n+1/2} + (dt/2) * (D^n + D^{n+1,*} + S^n + S^{n+1,*})
       do n=1,nlevs
          call multifab_copy_c(snew(n),1,sold(n),1,2,0)
          call multifab_mult_mult_s_c(s_update(n),1,dt,2,0)
          call multifab_mult_mult_s_c(bds_force(n),1,dt,2,0)
          call multifab_plus_plus_c(snew(n),1,s_update(n),1,2,0)
          call multifab_plus_plus_c(snew(n),1,bds_force(n),1,2,0)
       end do

    else

       ! compute A^{*,n+1} for scalars
       call mk_advective_s_fluxdiv(mla,umac,s_fc,s_update,dx,1,2)

       ! snew = s^{n+1} 
       !      = (1/2)*s^n + (1/2)*s^{*,n+1} + (dt/2)*(A^{*,n+1} + D^{*,n+1} + St^{*,n+1})
       do n=1,nlevs
          call multifab_plus_plus_c(snew(n),1,sold(n),1,2,0)
          call multifab_mult_mult_s_c(snew(n),1,0.5d0,2,0)
          call multifab_mult_mult_s_c(s_update(n),1,dt/2.d0,2,0)
          call multifab_plus_plus_c(snew(n),1,s_update(n),1,2,0)
       end do

    end if

    ! compute prim^{n+1} from s^{n+1} in valid region
    call convert_cons_to_prim(mla,snew,prim,.true.)

    ! fill ghost cells for prim^{n+1}
    do n=1,nlevs
       call multifab_fill_boundary(prim(n))
       call multifab_physbc(prim(n),2,scal_bc_comp+1,1,the_bc_tower%bc_tower_array(n),dx(n,:))
       call fill_rho_ghost_cells(prim(n),the_bc_tower%bc_tower_array(n))
    end do

    ! convert prim^{n+1} to s^{n+1} in valid and ghost region
    ! now snew = s^{n+1} has properly filled ghost cells
    call convert_cons_to_prim(mla,snew,prim,.false.)

    ! compute s^{n+1} to faces
    call average_cc_to_face(nlevs,snew,s_fc,1,scal_bc_comp,2,the_bc_tower%bc_tower_array)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 6 - Calculate Diffusive and Stochastic Fluxes
    ! Step 7 - Corrector Crank-Nicolson Step
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! build up rhs_v for gmres solve: first set gmres_rhs_v to mold/dt
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(gmres_rhs_v(n,i),1,mold(n,i),1,1,0)
          call multifab_mult_mult_s_c(gmres_rhs_v(n,i),1,1.d0/dt,1,0)
       end do
    end do

    ! compute grad p^{*,n+1}
    call compute_grad(mla,pres,gradp,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

    ! subtract grad p^{*,n+1} from gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradp(n,i),1,1,0)
       end do
    end do

    ! m_a_fluxdiv already contains A^n for momentum
    ! add A^{*,n+1} = -rho^{*,n+1} v^{*,n+1} v^{*,n+1} for momentum to m_a_fluxdiv
    call mk_advective_m_fluxdiv(mla,umac,mnew,m_a_fluxdiv,dx, &
                                the_bc_tower%bc_tower_array)

    ! add (1/2) m_a_fluxdiv to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_s_c(m_a_fluxdiv(n,i),1,0.5d0,1,0)
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_a_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! m_d_fluxdiv already contains (1/2) A_0^n v^n
    ! add (1/2) A_0^n v^n to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_d_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! compute div(Sigma^n') by incrementing existing stochastic flux and dividing by 2
    call stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,m_s_fluxdiv,eta,eta_ed,dx,dt,weights)
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_s_c(m_s_fluxdiv(n,i),1,0.5d0,1,0)
       end do
    end do

    ! add div(Sigma^n') to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_s_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! add gravity term
    if (any(grav(1:dm) .ne. 0.d0)) then
       call mk_grav_force(mla,gmres_rhs_v,s_fc_old,s_fc,the_bc_tower)
    end if

    ! initialize rhs_p for gmres solve to zero
    do n=1,nlevs
       call setval(gmres_rhs_p(n),0.d0,all=.true.)
    end do

    ! compute (chi,eta,kappa)^{n+1}
    call compute_chi(mla,chi,chi_fc,prim,dx,the_bc_tower%bc_tower_array)
    call compute_eta(mla,eta,eta_ed,prim,dx,the_bc_tower%bc_tower_array)
    call compute_kappa(mla,kappa,prim,dx)

    ! reset inhomogeneous bc condition to deal with reservoirs
    call set_inhomogeneous_vel_bcs(mla,vel_bc_n,vel_bc_t,eta_ed,dx, &
                                   the_bc_tower%bc_tower_array)

    ! reset to zero since we only add to them
    ! we store these for use at beginning of next time step
    do n=1,nlevs
       call setval(rhoc_fluxdiv(n),0.d0,all=.true.)
    end do

    ! set rhoc_d_fluxdiv to div(rho*chi grad c)^{n+1}
    call diffusive_rhoc_fluxdiv(mla,rhoc_fluxdiv,1,prim,s_fc,chi_fc,dx, &
                                the_bc_tower%bc_tower_array,vel_bc_n)

    if (barodiffusion_type .gt. 0) then
       ! compute baro-diffusion flux divergence
       call mk_baro_fluxdiv(mla,rhoc_fluxdiv,1,s_fc,chi_fc,gp_fc,dx, &
                            the_bc_tower%bc_tower_array,vel_bc_n)
    end if

    ! fill the stochastic multifabs with a new set of random numbers
    call fill_m_stochastic(mla)
    call fill_rhoc_stochastic(mla)

    ! create div(Psi^{n+1})
    call stochastic_rhoc_fluxdiv(mla,the_bc_tower%bc_tower_array,rhoc_fluxdiv, &
                                 s_fc,chi_fc,dx,dt,vel_bc_n,weights)

    ! add diffusive, baro-diffusion, and stochastic rho*c fluxes to rhs_p
    do n=1,nlevs
       call multifab_plus_plus_c(gmres_rhs_p(n),1,rhoc_fluxdiv(n),1,1,0)
    end do

    ! add external forcing for rho*c
    call mk_external_s_force(mla,gmres_rhs_p,dx,time+dt,1)

    ! modify umac to respect the boundary conditions we want after the next gmres solve
    ! thus when we add (1/2) A_0^{n+1} vbar^{*,n+1} and -rho^{n+1}*vbar^{*,n+1}/dt to gmres_rhs_v and 
    ! div vbar^{*,n+1} to gmres_rhs_p
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

   ! compute mtemp = rho^{*,n+1} * vbar^{*,n+1}
   call convert_m_to_umac(mla,s_fc,mtemp,umac,.false.)

   do n=1,nlevs
      do i=1,dm

         ! multiply mtemp by 1/dt
         call multifab_mult_mult_s_c(mtemp(n,i),1,1.d0/dt,1,0)

         ! subtract rho^{*,n+1} * vbar^{*,n+1} / dt from gmres_rhs_v
         call multifab_sub_sub_c(gmres_rhs_v(n,i),1,mtemp(n,i),1,1,0)


      end do
   end do

    ! set m_d_fluxdiv = A_0^{n+1} vbar^{n+1,*}
    do n=1,nlevs
       do i=1,dm
          call setval(m_d_fluxdiv(n,i),0.d0,all=.true.)
       end do
    end do
    call diffusive_m_fluxdiv(mla,m_d_fluxdiv,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    ! add (1/2) A_0^{n+1} vbar^{n+1,*} to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_s_c(m_d_fluxdiv(n,i),1,0.5d0,1,0)
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_d_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! compute div(v^{n+1,*})
    call compute_div(mla,umac,divu,dx,1,1,1)


    ! multiply gmres_rhs_p -S_fac
    do n=1,nlevs
       call multifab_mult_mult_s_c(gmres_rhs_p(n),1,-S_fac,1,0)
    end do

    ! add div(v^{n+1,*}) to gmres_rhs_p
    ! now gmres_rhs_p = div(v^{n+1,*}) - S^{n+1}
    ! the sign convention is correct since we solve -div(delta v) = gmres_rhs_p
    do n=1,nlevs
       call multifab_plus_plus_c(gmres_rhs_p(n),1,divu(n),1,1,0)
    end do

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
          call multifab_setval(dp(n),0.d0,all=.true.)
    end do

    do n=1,nlevs
       call zero_edgeval_physical(gmres_rhs_v(n,:),1,1,the_bc_tower%bc_tower_array(n))
    end do

    ! call gmres to compute delta v and delta p
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dp,s_fc, &
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
    ! compute p^{n+1} = p^{n+1,*} + dp
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

    ! convert v^{n+1} to m^{n+1} in valid and ghost region
    ! now mnew has properly filled ghost cells
    call convert_m_to_umac(mla,s_fc,mnew,umac,.false.)

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
          call multifab_destroy(mtemp(n,i))
          call multifab_destroy(gmres_rhs_v(n,i))
          call multifab_destroy(m_a_fluxdiv(n,i))
          call multifab_destroy(m_d_fluxdiv(n,i))
          call multifab_destroy(m_s_fluxdiv(n,i))
          call multifab_destroy(dumac(n,i))
          call multifab_destroy(umac_tmp(n,i))
          call multifab_destroy(gradp(n,i))
          call multifab_destroy(s_fc_old(n,i))
          call multifab_destroy(chi_fc_old(n,i))
       end do
    end do

    deallocate(weights)

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
