module advance_timestep_module

  use ml_layout_module
  use multifab_module
  use define_bc_module
  use convert_stag_module
  use convert_variables_module
  use convert_to_homogeneous_module
  use mk_advective_s_fluxdiv_module
  use mk_advective_m_fluxdiv_module
  use mk_baro_fluxdiv_module
  use mk_diffusive_fluxdiv_module
  use mk_grav_force_module
  use mk_stochastic_fluxdiv_module
  use gmres_module
  use init_module
  use div_and_grad_module
  use bc_module
  use multifab_physbc_module
  use multifab_physbc_stag_module
  use probin_lowmach_module, only: nscal, rhobar, diff_coef, visc_coef, grav
  use probin_common_module, only: fixed_dt
  use probin_module, only: use_barodiffusion

  use analysis_module

  implicit none

  private

  public :: advance_timestep

contains

  subroutine advance_timestep(mla,mold,mnew,umac,sold,snew,s_fc,prim,pold,pnew,chi,chi_fc, &
                              eta,eta_ed,kappa,rhoc_d_fluxdiv,rhoc_s_fluxdiv,rhoc_b_fluxdiv, &
                              gp_fc,dx,the_bc_tower,vel_bc_n,vel_bc_t)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: mold(:,:)
    type(multifab) , intent(inout) :: mnew(:,:)
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(inout) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(inout) :: s_fc(:,:)
    type(multifab) , intent(inout) :: prim(:)
    type(multifab) , intent(inout) :: pold(:)
    type(multifab) , intent(inout) :: pnew(:)
    type(multifab) , intent(inout) :: chi(:)
    type(multifab) , intent(inout) :: chi_fc(:,:)
    type(multifab) , intent(inout) :: eta(:)
    type(multifab) , intent(inout) :: eta_ed(:,:) ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(inout) :: kappa(:)
    type(multifab) , intent(inout) :: rhoc_d_fluxdiv(:)
    type(multifab) , intent(inout) :: rhoc_s_fluxdiv(:)
    type(multifab) , intent(inout) :: rhoc_b_fluxdiv(:)
    type(multifab) , intent(inout) :: gp_fc(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(inout) :: vel_bc_n(:,:)
    type(multifab) , intent(inout) :: vel_bc_t(:,:)

    ! local
    type(multifab) ::    s_update(mla%nlevel)
    type(multifab) :: gmres_rhs_p(mla%nlevel)
    type(multifab) ::          dp(mla%nlevel)
    type(multifab) ::        divu(mla%nlevel)

    type(multifab) ::       mtemp(mla%nlevel,mla%dim)
    type(multifab) :: gmres_rhs_v(mla%nlevel,mla%dim)
    type(multifab) :: m_a_fluxdiv_old(mla%nlevel,mla%dim)
    type(multifab) :: m_d_fluxdiv_old(mla%nlevel,mla%dim)
    type(multifab) :: m_a_fluxdiv_new(mla%nlevel,mla%dim)
    type(multifab) :: m_d_fluxdiv_new(mla%nlevel,mla%dim)
    type(multifab) :: m_s_fluxdiv    (mla%nlevel,mla%dim)
    type(multifab) ::        dumac(mla%nlevel,mla%dim)
    type(multifab) ::     umac_old(mla%nlevel,mla%dim)
    type(multifab) ::        gradp(mla%nlevel,mla%dim)
    type(multifab) ::     s_fc_old(mla%nlevel,mla%dim)
    type(multifab) ::   chi_fc_old(mla%nlevel,mla%dim)
    type(multifab) :: vel_bc_n_old(mla%nlevel,mla%dim)
    type(multifab) :: vel_bc_n_delta(mla%nlevel,mla%dim)
    type(multifab), allocatable :: vel_bc_t_old(:,:)
    type(multifab), allocatable :: vel_bc_t_delta(:,:)

    integer :: i,dm,n,nlevs
    logical :: nodal_temp(mla%dim)

    real(kind=dp_t) :: S_fac, theta_fac

    nlevs = mla%nlevel
    dm = mla%dim

    theta_fac = 1.d0/fixed_dt

    allocate(vel_bc_t_old  (nlevs,size(vel_bc_t,dim=2)))
    allocate(vel_bc_t_delta(nlevs,size(vel_bc_t,dim=2)))

    S_fac = (1.d0/rhobar(1) - 1.d0/rhobar(2))
    
    do n=1,nlevs
       call multifab_build(   s_update(n),mla%la(n),nscal,0)
       call multifab_build(gmres_rhs_p(n),mla%la(n),1    ,0)
       call multifab_build(         dp(n),mla%la(n),1    ,1)
       call multifab_build(       divu(n),mla%la(n),1    ,0)
       do i=1,dm
          call multifab_build_edge(          mtemp(n,i),mla%la(n),nscal,1,i)
          call multifab_build_edge(    gmres_rhs_v(n,i),mla%la(n),1    ,0,i)
          call multifab_build_edge(m_a_fluxdiv_old(n,i),mla%la(n),1    ,0,i)
          call multifab_build_edge(m_d_fluxdiv_old(n,i),mla%la(n),1    ,0,i)
          call multifab_build_edge(m_a_fluxdiv_new(n,i),mla%la(n),1    ,0,i)
          call multifab_build_edge(m_d_fluxdiv_new(n,i),mla%la(n),1    ,0,i)
          call multifab_build_edge(m_s_fluxdiv    (n,i),mla%la(n),1    ,0,i)
          call multifab_build_edge(          dumac(n,i),mla%la(n),1    ,1,i)
          call multifab_build_edge(       umac_old(n,i),mla%la(n),1    ,1,i)
          call multifab_build_edge(          gradp(n,i),mla%la(n),1    ,0,i)
          call multifab_build_edge(       s_fc_old(n,i),mla%la(n),nscal,1,i)
          call multifab_build_edge(     chi_fc_old(n,i),mla%la(n),1    ,0,i)
          call multifab_build_edge(vel_bc_n_old(n,i),mla%la(n),1,0,i)
          call multifab_build_edge(vel_bc_n_delta(n,i),mla%la(n),1,0,i)
       end do
       
       if (dm .eq. 2) then
          ! y-velocity bc on x-faces (nodal)
          call multifab_build_nodal(vel_bc_t_old(n,1)  ,mla%la(n),1,0)
          call multifab_build_nodal(vel_bc_t_delta(n,1),mla%la(n),1,0)
          ! x-velocity bc on y-faces (nodal)
          call multifab_build_nodal(vel_bc_t_old(n,2)  ,mla%la(n),1,0)
          call multifab_build_nodal(vel_bc_t_delta(n,2),mla%la(n),1,0)
       else
          ! y-velocity bc on x-faces (nodal in y and x)
          nodal_temp(1) = .true.
          nodal_temp(2) = .true.
          nodal_temp(3) = .false.
          call multifab_build(vel_bc_t_old(n,1)  ,mla%la(n),1,0,nodal_temp)
          call multifab_build(vel_bc_t_delta(n,1),mla%la(n),1,0,nodal_temp)
          ! z-velocity bc on x-faces (nodal in z and x)
          nodal_temp(1) = .true.
          nodal_temp(2) = .false.
          nodal_temp(3) = .true.
          call multifab_build(vel_bc_t_old(n,2)  ,mla%la(n),1,0,nodal_temp)
          call multifab_build(vel_bc_t_delta(n,2),mla%la(n),1,0,nodal_temp)
          ! x-velocity bc on y-faces (nodal in x and y)
          nodal_temp(1) = .true.
          nodal_temp(2) = .true.
          nodal_temp(3) = .false.
          call multifab_build(vel_bc_t_old(n,3)  ,mla%la(n),1,0,nodal_temp)
          call multifab_build(vel_bc_t_delta(n,3),mla%la(n),1,0,nodal_temp)
          ! z-velocity bc on y-faces (nodal in z and y)
          nodal_temp(1) = .false.
          nodal_temp(2) = .true.
          nodal_temp(3) = .true.
          call multifab_build(vel_bc_t_old(n,4)  ,mla%la(n),1,0,nodal_temp)
          call multifab_build(vel_bc_t_delta(n,4),mla%la(n),1,0,nodal_temp)
          ! x-velocity bc on z-faces (nodal in x and z)
          nodal_temp(1) = .true.
          nodal_temp(2) = .false.
          nodal_temp(3) = .true.
          call multifab_build(vel_bc_t_old(n,5)  ,mla%la(n),1,0,nodal_temp)
          call multifab_build(vel_bc_t_delta(n,5),mla%la(n),1,0,nodal_temp)
          ! y-velocity bc on z-faces (nodal in y and z)
          nodal_temp(1) = .false.
          nodal_temp(2) = .true.
          nodal_temp(3) = .true.
          call multifab_build(vel_bc_t_old(n,6)  ,mla%la(n),1,0,nodal_temp)
          call multifab_build(vel_bc_t_delta(n,6),mla%la(n),1,0,nodal_temp)
       end if

    end do

    do n=1,nlevs
       call setval(s_update(n),0.d0,all=.true.)
       do i=1,dm
          call setval(m_a_fluxdiv_old(n,i),0.d0,all=.true.)
          call setval(m_d_fluxdiv_old(n,i),0.d0,all=.true.)
          call setval(m_a_fluxdiv_new(n,i),0.d0,all=.true.)
          call setval(m_d_fluxdiv_new(n,i),0.d0,all=.true.)
          call setval(m_s_fluxdiv    (n,i),0.d0,all=.true.)
          call setval(          dumac(n,i),0.d0,all=.true.)
       end do
    end do

    ! make copies of old quantities
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(  umac_old(n,i),1,  umac(n,i),1,1    ,1)
          call multifab_copy_c(  s_fc_old(n,i),1,  s_fc(n,i),1,nscal,1)
          call multifab_copy_c(chi_fc_old(n,i),1,chi_fc(n,i),1,1    ,0)
          call multifab_copy_c(vel_bc_n_old(n,i),1,vel_bc_n(n,i),1,1,0)
       end do
       do i=1,size(vel_bc_t,dim=2)
          call multifab_copy_c(vel_bc_t_old(n,i),1,vel_bc_t(n,i),1,1,0)
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 1 - Forward-Euler Scalar Predictor
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! set s_update to A^n for scalars
    call mk_advective_s_fluxdiv(mla,umac_old,s_fc,s_update,dx,1,nscal)

    ! add D^n  for rho1 to s_update
    ! add St^n for rho1 to s_update
    ! add baro-diffusion^n to s_update
    do n=1,nlevs
       call multifab_plus_plus_c(s_update(n),2,rhoc_d_fluxdiv(n),1,1,0)
       call multifab_plus_plus_c(s_update(n),2,rhoc_s_fluxdiv(n),1,1,0)
       call multifab_plus_plus_c(s_update(n),2,rhoc_b_fluxdiv(n),1,1,0)
    end do

    ! set snew = s^{*,n+1} = s^n + dt * (A^n + D^n + St^n)
    do n=1,nlevs
       call multifab_mult_mult_s_c(s_update(n),1,fixed_dt,nscal,0)
       call multifab_copy_c(snew(n),1,sold(n),1,nscal,0)
       call multifab_plus_plus_c(snew(n),1,s_update(n),1,nscal,0)
    end do

    ! compute prim^{*,n+1} from s^{*,n+1} in valid region
    call convert_cons_to_prim(mla,snew,prim,.true.)

    ! fill ghost cells for prim^{*,n+1}
    do n=1,nlevs
       call multifab_fill_boundary(prim(n))
       call multifab_physbc(prim(n),1,scal_bc_comp,2,the_bc_tower%bc_tower_array(n),dx(n,:))
    end do

    ! convert prim^{*,n+1} to s^{*,n+1} in valid and ghost region
    ! now snew = s^{*,n+1} has properly filled ghost cells
    call convert_cons_to_prim(mla,snew,prim,.false.)

    ! compute s^{*,n+1} to faces
    do i=1,nscal
       call average_cc_to_face(nlevs,snew,s_fc,i,scal_bc_comp,1,the_bc_tower%bc_tower_array)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 2 - Crank-Nicolson Velocity Predictor
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! build up rhs_v for gmres solve: first set gmres_rhs_v to mold = m^n
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(gmres_rhs_v(n,i),1,mold(n,i),1,1,0)
       end do
    end do

    ! compute mtemp = rho^{*,n+1} * v^n
    call convert_m_to_umac(mla,s_fc,mtemp,umac_old,.false.)

    do n=1,nlevs
       do i=1,dm

          ! subtract rho^{*,n+1} * v^n from gmres_rhs_v
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,mtemp(n,i),1,1,0)

          ! multiply gmres_rhs_v by 1/dt
          call multifab_mult_mult_s_c(gmres_rhs_v(n,i),1,1.d0/fixed_dt,1,0)

       end do
    end do

    ! compute grad p
    call compute_grad(mla,pold,gradp,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

    ! subtract grad p from gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradp(n,i),1,1,0)
       end do
    end do

    ! compute m_a_fluxdiv_old = A^n for momentum
    call mk_advective_m_fluxdiv(mla,umac_old,mold,m_a_fluxdiv_old,dx, &
                                the_bc_tower%bc_tower_array)

    ! add A^n for momentum to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_a_fluxdiv_old(n,i),1,1,0)
       end do
    end do

    ! compute m_d_fluxdiv_old = A_0^n v^n
    call mk_diffusive_m_fluxdiv(mla,m_d_fluxdiv_old,umac_old,eta,eta_ed,kappa,dx, &
                                the_bc_tower%bc_tower_array)

    ! add (1/2) A_0^n v^n to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_s_c(m_d_fluxdiv_old(n,i),1,0.5d0,1,0)
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_d_fluxdiv_old(n,i),1,1,0)
       end do
    end do

    ! compute (chi,eta,kappa)^{*,n+1}
    call compute_chi(mla,chi,chi_fc,prim,dx,the_bc_tower%bc_tower_array)
    call compute_eta(mla,eta,eta_ed,prim,dx,the_bc_tower%bc_tower_array)
    call compute_kappa(mla,kappa,prim,dx)

    ! set m_d_fluxdiv_new = A_0^{*,n+1} v^n
    call mk_diffusive_m_fluxdiv(mla,m_d_fluxdiv_new,umac_old,eta,eta_ed,kappa,dx, &
                                the_bc_tower%bc_tower_array)

    ! add (1/2) A_0^{*,n+1} v^n to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_s_c(m_d_fluxdiv_new(n,i),1,0.5d0,1,0)
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_d_fluxdiv_new(n,i),1,1,0)
       end do
    end do

    ! compute m_s_fluxdiv = div(Sigma^n)
    call mk_stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,m_s_fluxdiv,eta,eta_ed,dx)

    ! add div(Sigma^n) to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_s_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! add gravity term
    if (any(grav(1:dm) .ne. 0.d0)) then
       call mk_grav_force(mla,gmres_rhs_v,s_fc_old,s_fc)
    end if

    ! initialize rhs_p for gmres solve to zero
    do n=1,nlevs
       call setval(gmres_rhs_p(n),0.d0,all=.true.)
    end do

    ! reset inhomogeneous bc condition to deal with reservoirs
    call set_inhomogeneous_vel_bcs(mla,vel_bc_n,vel_bc_t,eta_ed,dx, &
                                   the_bc_tower%bc_tower_array)

    ! add div(rho*chi grad c)^{*,n+1} to rhs_p
    call mk_diffusive_rhoc_fluxdiv(mla,gmres_rhs_p,1,prim,s_fc,chi_fc,dx, &
                                   the_bc_tower%bc_tower_array,vel_bc_n)

    if (use_barodiffusion) then
       ! add baro-diffusion flux diveregnce to rhs_p
       call mk_baro_fluxdiv(mla,gmres_rhs_p,1,s_fc,chi_fc,gp_fc,dx, &
                            the_bc_tower%bc_tower_array,vel_bc_n)
    end if

    ! add div(Psi^n) to rhs_p
    call mk_stochastic_s_fluxdiv(mla,the_bc_tower%bc_tower_array,gmres_rhs_p,s_fc_old, &
                                 chi_fc_old,dx,vel_bc_n)

    do n=1,nlevs
       do i=1,dm
          ! compute change in normal velocity boundary condition over the time step
          ! this deals with time-dependent velocity boundary conditions
          call multifab_copy_c(vel_bc_n_delta(n,i),1,vel_bc_n(n,i),1,1,0)
          call multifab_sub_sub_c(vel_bc_n_delta(n,i),1,vel_bc_n_old(n,i),1,1,0)
       end do
       do i=1,size(vel_bc_t,dim=2)
          call multifab_copy_c(vel_bc_t_delta(n,i),1,vel_bc_t(n,i),1,1,0)
          call multifab_sub_sub_c(vel_bc_t_delta(n,i),1,vel_bc_t_old(n,i),1,1,0)
       end do
    end do

    ! reset s_update for all scalars to zero
    ! then, set s_update for rho1 to F^{*,n+1} = div(rho*chi grad c)^{*,n+1} + div(Psi^n)
    ! it is used in Step 3 below
    do n=1,nlevs
       call multifab_setval_c(s_update(n),0.d0,1,1,all=.true.)
       call multifab_copy_c(s_update(n),2,gmres_rhs_p(n),1,1,0)
    end do

    ! multiply gmres_rhs_p by -S_fac
    ! now gmres_rhs_p = -S^{*,n+1}
    do n=1,nlevs
       call multifab_mult_mult_s_c(gmres_rhs_p(n),1,-S_fac,1,0)
    end do

    ! compute div(v^n)
    call compute_div(mla,umac_old,divu,dx,1,1,1)

    ! add div(v^n) to gmres_rhs_p
    ! now gmres_rhs_p = div(v^n) - S^{*,n+1}
    ! the sign convention is correct since we solve -div(delta v) = gmres_rhs_p
    do n=1,nlevs
       call multifab_plus_plus_c(gmres_rhs_p(n),1,divu(n),1,1,0)
    end do

    ! multiply eta and kappa by 1/2 to put in proper form for gmres solve
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,1.d0/2.d0,1,eta(n)%ng)
       call multifab_mult_mult_s_c(kappa(n),1,1.d0/2.d0,1,kappa(n)%ng)
       if (dm .eq. 2) then
          call multifab_mult_mult_s_c(eta_ed(n,1),1,1.d0/2.d0,1,eta_ed(n,1)%ng)
       else if (dm .eq. 3) then
          call multifab_mult_mult_s_c(eta_ed(n,1),1,1.d0/2.d0,1,eta_ed(n,1)%ng)
          call multifab_mult_mult_s_c(eta_ed(n,2),1,1.d0/2.d0,1,eta_ed(n,2)%ng)
          call multifab_mult_mult_s_c(eta_ed(n,3),1,1.d0/2.d0,1,eta_ed(n,3)%ng)
       end if
    end do

    ! set the initial guess to zero
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(dumac(n,i),0.d0,all=.true.)
          call multifab_setval(   dp(n)  ,0.d0,all=.true.)
       end do
    end do

    ! for inhomogeneous boundary conditions, convert problem to homogeneous by
    ! subtracting from the RHS the result of the operator applied to a solution
    ! vector with zeros everywhere in the problem domain, and ghost cells filled to
    ! respect the boundary conditions
    call convert_to_homogeneous(mla,gmres_rhs_v,gmres_rhs_p,s_fc,eta,eta_ed, &
                                kappa,1.d0/fixed_dt,dx,the_bc_tower, &
                                vel_bc_n_delta,vel_bc_t_delta)

    ! call gmres to compute delta v and delta p
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dp,s_fc, &
               eta,eta_ed,kappa,theta_fac)

    ! restore eta and kappa
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,2.d0,1,eta(n)%ng)
       call multifab_mult_mult_s_c(kappa(n),1,2.d0,1,kappa(n)%ng)
       if (dm .eq. 2) then
          call multifab_mult_mult_s_c(eta_ed(n,1),1,2.d0,1,eta_ed(n,1)%ng)
       else if (dm .eq. 3) then
          call multifab_mult_mult_s_c(eta_ed(n,1),1,2.d0,1,eta_ed(n,1)%ng)
          call multifab_mult_mult_s_c(eta_ed(n,2),1,2.d0,1,eta_ed(n,2)%ng)
          call multifab_mult_mult_s_c(eta_ed(n,3),1,2.d0,1,eta_ed(n,3)%ng)
       end if
    end do

    ! compute v^{*,n+1} = v^n + delta v
    ! compute p^{*,n+1}= p^n + delta p
    ! keep both dp and dumac as initial guess for corrector
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(umac(n,i),1,dumac(n,i),1,1,0)
       end do
       call multifab_copy_c(pnew(n),1,pold(n),1,1,0)
       call multifab_plus_plus_c(pnew(n),1,dp(n),1,1,0)
    end do

    do n=1,nlevs
       ! presure ghost cells
       call multifab_fill_boundary(pnew(n))
       call multifab_physbc(pnew(n),1,pres_bc_comp,1,the_bc_tower%bc_tower_array(n),dx(n,:))
       do i=1,dm
          ! fill periodic and interior ghost cells
          call multifab_fill_boundary(umac(n,i))
          ! set normal velocity on physical domain boundaries
          call multifab_physbc_domainvel(umac(n,i),vel_bc_comp+i-1, &
                                         the_bc_tower%bc_tower_array(n), &
                                         dx(n,:),vel_bc_n(n,:))
          ! set transverse velocity behind physical boundaries
          call multifab_physbc_macvel(umac(n,i),vel_bc_comp+i-1, &
                                      the_bc_tower%bc_tower_array(n), &
                                      dx(n,:),vel_bc_t(n,:))
       end do
    end do

    if (use_barodiffusion) then
       ! compute grad p^{n+1,*}
       call compute_grad(mla,pnew,gp_fc,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)
    end if

    ! convert v^{*,n+1} to m^{*,n+1} in valid and ghost region
    ! now mnew has properly filled ghost cells
    call convert_m_to_umac(mla,s_fc,mnew,umac,.false.)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 3 - Trapezoidal Scalar Corrector
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! s_update already contains D^{*,n+1} + St^{*,n+1} for rho1 from above
    ! add A^{*,n+1} for s to s_update
    call mk_advective_s_fluxdiv(mla,umac,s_fc,s_update,dx,1,nscal)

    ! snew = s^{n+1} 
    !      = (1/2)*s^n + (1/2)*s^{*,n+1} + (dt/2)*(A^{*,n+1} + D^{*,n+1} + St^{*,n+1})
    do n=1,nlevs
       call multifab_plus_plus_c(snew(n),1,sold(n),1,nscal,0)
       call multifab_mult_mult_s_c(snew(n),1,0.5d0,nscal,0)
       call multifab_mult_mult_s_c(s_update(n),1,fixed_dt/2.d0,nscal,0)
       call multifab_plus_plus_c(snew(n),1,s_update(n),1,nscal,0)
    end do

    ! compute prim^{n+1} from s^{n+1} in valid region
    call convert_cons_to_prim(mla,snew,prim,.true.)

    ! fill ghost cells for prim^{n+1}
    do n=1,nlevs
       call multifab_fill_boundary(prim(n))
       call multifab_physbc(prim(n),1,scal_bc_comp,2,the_bc_tower%bc_tower_array(n),dx(n,:))
    end do

    ! convert prim^{n+1} to s^{n+1} in valid and ghost region
    ! now snew = s^{n+1} has properly filled ghost cells
    call convert_cons_to_prim(mla,snew,prim,.false.)

    ! compute s^{n+1} to faces
    do i=1,nscal
       call average_cc_to_face(nlevs,snew,s_fc,i,scal_bc_comp,1,the_bc_tower%bc_tower_array)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 4 - Crank-Nicolson Velocity Corrector
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! build up rhs_v for gmres solve: first set rhs to mold = m^n
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(gmres_rhs_v(n,i),1,mold(n,i),1,1,0)
       end do
    end do

    ! compute mtemp = rho^{n+1} * v^n
    call convert_m_to_umac(mla,s_fc,mtemp,umac_old,.false.)

    do n=1,nlevs
       do i=1,dm

          ! subtract rho^{n+1} * v^n from gmres_rhs_v
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,mtemp(n,i),1,1,0)

          ! multiply gmres_rhs_v by 1/dt
          call multifab_mult_mult_s_c(gmres_rhs_v(n,i),1,1.d0/fixed_dt,1,0)

       end do
    end do

    ! gradp already contains grad p
    ! subtract grad p^n from gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradp(n,i),1,1,0)
       end do
    end do

    ! m_a_fluxdiv_old already contains A^n for momentum
    ! add (1/2) A^n gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_s_c(m_a_fluxdiv_old(n,i),1,0.5d0,1,0)
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_a_fluxdiv_old(n,i),1,1,0)
       end do
    end do

    ! compute m_a_fluxdiv_new = A(s^{n+1},v^{n+1,*}) for momentum
    call mk_advective_m_fluxdiv(mla,umac,mnew,m_a_fluxdiv_new,dx, &
                                the_bc_tower%bc_tower_array)

    ! add (1/2) m_a_fluxdiv_new to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_s_c(m_a_fluxdiv_new(n,i),1,0.5d0,1,0)
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_a_fluxdiv_new(n,i),1,1,0)
       end do
    end do

    ! m_d_fluxdiv_old already contains (1/2) A_0^n v^n
    ! add (1/2) A_0^n v^n to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_d_fluxdiv_old(n,i),1,1,0)
       end do
    end do

    ! compute (chi,eta,kappa)^{n+1}
    call compute_chi(mla,chi,chi_fc,prim,dx,the_bc_tower%bc_tower_array)
    call compute_eta(mla,eta,eta_ed,prim,dx,the_bc_tower%bc_tower_array)
    call compute_kappa(mla,kappa,prim,dx)

    ! reset m_d_fluxdiv_new
    do n=1,nlevs
       do i=1,dm
          call setval(m_d_fluxdiv_new(n,i),0.d0,all=.true.)
       end do
    end do

    ! set m_d_fluxdiv_new = A_0^{n+1} v^n
    call mk_diffusive_m_fluxdiv(mla,m_d_fluxdiv_new,umac_old,eta,eta_ed,kappa,dx, &
                                the_bc_tower%bc_tower_array)

    ! add (1/2) A_0^{n+1} v^n to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_s_c(m_d_fluxdiv_new(n,i),1,0.5d0,1,0)
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_d_fluxdiv_new(n,i),1,1,0)
       end do
    end do

    ! add div(Sigma^n) to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_s_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! add gravity term
    if (any(grav(1:dm) .ne. 0.d0)) then
       call mk_grav_force(mla,gmres_rhs_v,s_fc_old,s_fc)
    end if

    ! initialize rhs_p for gmres solve to zero
    do n=1,nlevs
       call setval(gmres_rhs_p(n),0.d0,all=.true.)
    end do

    ! reset to zero since we only add to them
    ! we store these for use at beginning of next time step
    do n=1,nlevs
       call setval(rhoc_d_fluxdiv(n),0.d0,all=.true.)
       call setval(rhoc_s_fluxdiv(n),0.d0,all=.true.)
       call setval(rhoc_b_fluxdiv(n),0.d0,all=.true.)
    end do

    ! reset inhomogeneous bc condition to deal with reservoirs
    call set_inhomogeneous_vel_bcs(mla,vel_bc_n,vel_bc_t,eta_ed,dx, &
                                   the_bc_tower%bc_tower_array)

    ! set rhoc_d_fluxdiv to div(rho*chi grad c)^{n+1}
    call mk_diffusive_rhoc_fluxdiv(mla,rhoc_d_fluxdiv,1,prim,s_fc,chi_fc,dx, &
                                   the_bc_tower%bc_tower_array,vel_bc_n)

    ! add div(rho*chi grad c)^{n+1} to rhs_p
    do n=1,nlevs
       call multifab_plus_plus_c(gmres_rhs_p(n),1,rhoc_d_fluxdiv(n),1,1,0)
    end do

    if (use_barodiffusion) then
       ! compute baro-diffusion flux divergence
       call mk_baro_fluxdiv(mla,rhoc_b_fluxdiv,1,s_fc,chi_fc,gp_fc,dx, &
                            the_bc_tower%bc_tower_array,vel_bc_n)

       ! add baro-diffusion flux divergence to rhs_p
       do n=1,nlevs
          call multifab_plus_plus_c(gmres_rhs_p(n),1,rhoc_b_fluxdiv(n),1,1,0)
       end do
    end if

    ! fill the stochastic multifabs with a new set of random numbers
    call fill_stochastic(mla)

    ! create div(Psi^{n+1})
    call mk_stochastic_s_fluxdiv(mla,the_bc_tower%bc_tower_array,rhoc_s_fluxdiv, &
                                 s_fc,chi_fc,dx,vel_bc_n)

    ! add div(Psi^{n+1}) to rhs_p
    do n=1,nlevs
       call multifab_plus_plus_c(gmres_rhs_p(n),1,rhoc_s_fluxdiv(n),1,1,0)
    end do    

    ! multiply gmres_rhs_p -S_fac
    do n=1,nlevs
       call multifab_mult_mult_s_c(gmres_rhs_p(n),1,-S_fac,1,0)
    end do

    ! divu already contains div(v^n)
    ! add div(v^n) to gmres_rhs_p
    ! now gmres_rhs_p = div(v^n) - S^{n+1}
    ! the sign convention is correct since we solve -div(delta v) = gmres_rhs_p
    do n=1,nlevs
       call multifab_plus_plus_c(gmres_rhs_p(n),1,divu(n),1,1,0)
    end do

    ! multiply eta and kappa by 1/2 to put in proper form for gmres solve
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,1.d0/2.d0,1,eta(n)%ng)
       call multifab_mult_mult_s_c(kappa(n),1,1.d0/2.d0,1,eta(n)%ng)
       if (dm .eq. 2) then
          call multifab_mult_mult_s_c(eta_ed(n,1),1,1.d0/2.d0,1,eta_ed(n,1)%ng)
       else if (dm .eq. 3) then
          call multifab_mult_mult_s_c(eta_ed(n,1),1,1.d0/2.d0,1,eta_ed(n,1)%ng)
          call multifab_mult_mult_s_c(eta_ed(n,2),1,1.d0/2.d0,1,eta_ed(n,2)%ng)
          call multifab_mult_mult_s_c(eta_ed(n,3),1,1.d0/2.d0,1,eta_ed(n,3)%ng)
       end if
    end do

    do n=1,nlevs
       do i=1,dm
          ! compute change in normal velocity boundary condition over the time step
          ! this deals with time-dependent velocity boundary conditions
          call multifab_copy_c(vel_bc_n_delta(n,i),1,vel_bc_n(n,i),1,1,0)
          call multifab_sub_sub_c(vel_bc_n_delta(n,i),1,vel_bc_n_old(n,i),1,1,0)
       end do
       do i=1,size(vel_bc_t,dim=2)
          call multifab_copy_c(vel_bc_t_delta(n,i),1,vel_bc_t(n,i),1,1,0)
          call multifab_sub_sub_c(vel_bc_t_delta(n,i),1,vel_bc_t_old(n,i),1,1,0)
       end do
    end do

    ! for inhomogeneous boundary conditions, convert problem to homogeneous by
    ! subtracting from the RHS the result of the operator applied to a solution
    ! vector with zeros everywhere in the problem domain, and ghost cells filled to
    ! respect the boundary conditions
    call convert_to_homogeneous(mla,gmres_rhs_v,gmres_rhs_p,s_fc,eta,eta_ed, &
                                kappa,1.d0/fixed_dt,dx,the_bc_tower, &
                                vel_bc_n_delta,vel_bc_t_delta)

    ! call gmres to compute delta v and delta p
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dp,s_fc, &
               eta,eta_ed,kappa,theta_fac)

    ! restore eta and kappa
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,2.d0,1,eta(n)%ng)
       call multifab_mult_mult_s_c(kappa(n),1,2.d0,1,kappa(n)%ng)
       if (dm .eq. 2) then
          call multifab_mult_mult_s_c(eta_ed(n,1),1,2.d0,1,eta_ed(n,1)%ng)
       else if (dm .eq. 3) then
          call multifab_mult_mult_s_c(eta_ed(n,1),1,2.d0,1,eta_ed(n,1)%ng)
          call multifab_mult_mult_s_c(eta_ed(n,2),1,2.d0,1,eta_ed(n,1)%ng)
          call multifab_mult_mult_s_c(eta_ed(n,3),1,2.d0,1,eta_ed(n,1)%ng)
       end if
    end do

    ! compute v^{n+1} = v^n + dumac
    ! compute p^{n+1} = p^n + dp
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(umac(n,i),1,umac_old(n,i),1,1,0)
          call multifab_plus_plus_c(umac(n,i),1,dumac(n,i),1,1,0)
       end do
       call multifab_copy_c(pnew(n),1,pold(n),1,1,0)
       call multifab_plus_plus_c(pnew(n),1,dp(n),1,1,0)
    end do

    do n=1,nlevs
       ! presure ghost cells
       call multifab_fill_boundary(pnew(n))
       call multifab_physbc(pnew(n),1,pres_bc_comp,1,the_bc_tower%bc_tower_array(n),dx(n,:))
       do i=1,dm
          ! fill periodic and interior ghost cells
          call multifab_fill_boundary(umac(n,i))
          ! set normal velocity on physical domain boundaries
          call multifab_physbc_domainvel(umac(n,i),vel_bc_comp+i-1, &
                                         the_bc_tower%bc_tower_array(n), &
                                         dx(n,:),vel_bc_n(n,:))
          ! set transverse velocity behind physical boundaries
          call multifab_physbc_macvel(umac(n,i),vel_bc_comp+i-1, &
                                      the_bc_tower%bc_tower_array(n), &
                                      dx(n,:),vel_bc_t(n,:))
       end do
    end do

    ! convert v^{n+1} to m^{n+1} in valid and ghost region
    ! now mnew has properly filled ghost cells
    call convert_m_to_umac(mla,s_fc,mnew,umac,.false.)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End Time-Advancement
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do n=1,nlevs
       call multifab_destroy(s_update(n))
       call multifab_destroy(gmres_rhs_p(n))
       call multifab_destroy(dp(n))
       call multifab_destroy(divu(n))
       do i=1,dm
          call multifab_destroy(mtemp(n,i))
          call multifab_destroy(gmres_rhs_v(n,i))
          call multifab_destroy(m_a_fluxdiv_old(n,i))
          call multifab_destroy(m_d_fluxdiv_old(n,i))
          call multifab_destroy(m_a_fluxdiv_new(n,i))
          call multifab_destroy(m_d_fluxdiv_new(n,i))
          call multifab_destroy(m_s_fluxdiv(n,i))
          call multifab_destroy(dumac(n,i))
          call multifab_destroy(umac_old(n,i))
          call multifab_destroy(gradp(n,i))
          call multifab_destroy(s_fc_old(n,i))
          call multifab_destroy(chi_fc_old(n,i))
          call multifab_destroy(vel_bc_n_old(n,i))
          call multifab_destroy(vel_bc_n_delta(n,i))
       end do
       do i=1,size(vel_bc_t,dim=2)
          call multifab_destroy(vel_bc_t_old(n,i))
          call multifab_destroy(vel_bc_t_delta(n,i))
       end do
    end do

    deallocate(vel_bc_t_old,vel_bc_t_delta)

  end subroutine advance_timestep

end module advance_timestep_module
