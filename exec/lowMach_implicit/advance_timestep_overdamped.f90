module advance_timestep_overdamped_module

  use ml_layout_module
  use multifab_module
  use define_bc_module
  use convert_stag_module
  use convert_variables_module
  use convert_to_homogeneous_module
  use mk_advective_s_fluxdiv_module
  use mk_advective_m_fluxdiv_module
  use mk_diffusive_fluxdiv_module
  use mk_grav_force_module
  use mk_stochastic_fluxdiv_module
  use bds_module
  use gmres_module
  use init_module
  use div_and_grad_module
  use bc_module
  use multifab_physbc_module
  use multifab_physbc_stag_module
  use probin_lowmach_module, only: nscal, rhobar, grav
  use probin_common_module, only: advection_type

  use analysis_module

  implicit none

  private

  public :: advance_timestep_overdamped

contains

  subroutine advance_timestep_overdamped(mla,mnew,umac,sold,snew,s_fc,prim,pold,pnew, &
                                         chi,chi_fc,eta,eta_ed,kappa,dx,dt,the_bc_tower, &
                                         vel_bc_n,vel_bc_t)

    type(ml_layout), intent(in   ) :: mla
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
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(inout) :: vel_bc_n(:,:)
    type(multifab) , intent(inout) :: vel_bc_t(:,:)

    ! local
    type(multifab) ::    s_update(mla%nlevel)
    type(multifab) ::   bds_force(mla%nlevel)
    type(multifab) :: gmres_rhs_p(mla%nlevel)
    type(multifab) ::          dp(mla%nlevel)
    type(multifab) ::        divu(mla%nlevel)

    type(multifab) ::  gmres_rhs_v(mla%nlevel,mla%dim)
    type(multifab) ::        dumac(mla%nlevel,mla%dim)
    type(multifab) ::        gradp(mla%nlevel,mla%dim)
    type(multifab) :: vel_bc_n_old(mla%nlevel,mla%dim)
    type(multifab) :: vel_bc_n_delta(mla%nlevel,mla%dim)
    type(multifab), allocatable :: vel_bc_t_old(:,:)
    type(multifab), allocatable :: vel_bc_t_delta(:,:)

    integer :: i,dm,n,nlevs
    logical :: nodal_temp(mla%dim)

    real(kind=dp_t) :: S_fac, theta_fac

    nlevs = mla%nlevel
    dm = mla%dim

    theta_fac = 0.d0

    allocate(vel_bc_t_old  (nlevs,size(vel_bc_t,dim=2)))
    allocate(vel_bc_t_delta(nlevs,size(vel_bc_t,dim=2)))

    S_fac = (1.d0/rhobar(1) - 1.d0/rhobar(2))
    
    do n=1,nlevs
       call multifab_build(   s_update(n),mla%la(n),nscal,0)
       call multifab_build(  bds_force(n),mla%la(n),nscal,1)
       call multifab_build(gmres_rhs_p(n),mla%la(n),1    ,0)
       call multifab_build(         dp(n),mla%la(n),1    ,1)
       call multifab_build(       divu(n),mla%la(n),1    ,0)
       do i=1,dm
          call multifab_build_edge(    gmres_rhs_v(n,i),mla%la(n),1    ,0,i)
          call multifab_build_edge(          dumac(n,i),mla%la(n),1    ,1,i)
          call multifab_build_edge(          gradp(n,i),mla%la(n),1    ,0,i)
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
       call setval(bds_force(n),0.d0,all=.true.)
       do i=1,dm
          call setval(dumac(n,i),0.d0,all=.true.)
       end do
    end do

    ! save boundary conditions from corrector
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(vel_bc_n_old(n,i),1,vel_bc_n(n,i),1,1,0)
       end do
       do i=1,size(vel_bc_t,dim=2)
          call multifab_copy_c(vel_bc_t_old(n,i),1,vel_bc_t(n,i),1,1,0)
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Steps 1 and 2 - Predictor Stokes Solve
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! build up rhs_v for gmres solve
    do n=1,nlevs
       do i=1,dm
          call setval(gmres_rhs_v(n,i),0.d0,all=.true.)
       end do
    end do

    ! compute grad p^{n-1/2}
    call compute_grad(mla,pold,gradp,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

    ! subtract grad p^{n-1/2} from gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradp(n,i),1,1,0)
       end do
    end do

    ! add A_0^n v^{n-1/2} to gmres_rhs_v
    call mk_diffusive_m_fluxdiv(mla,gmres_rhs_v,umac,eta,eta_ed,kappa,dx, &
                                the_bc_tower%bc_tower_array)

    ! add div(Sigma^n) to gmres_rhs_v
    call mk_stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,gmres_rhs_v, &
                                 eta,eta_ed,dx,dt)

    ! add gravity term to gmres_rhs_v
    if (any(grav(1:dm) .ne. 0.d0)) then
       call mk_grav_force(mla,gmres_rhs_v,s_fc,s_fc)
    end if

    ! initialize rhs_p for gmres solve to zero
    do n=1,nlevs
       call setval(gmres_rhs_p(n),0.d0,all=.true.)
    end do

    ! reset inhomogeneous bc condition to deal with reservoirs
    call set_inhomogeneous_vel_bcs(mla,vel_bc_n,vel_bc_t,eta_ed,dx, &
                                   the_bc_tower%bc_tower_array)

    ! add div(rho*chi grad c)^n to rhs_p
    call mk_diffusive_rhoc_fluxdiv(mla,gmres_rhs_p,1,prim,s_fc,chi_fc,dx, &
                                   the_bc_tower%bc_tower_array,vel_bc_n)

    ! add div(Psi^n) to rhs_p
    call mk_stochastic_s_fluxdiv(mla,the_bc_tower%bc_tower_array,gmres_rhs_p,s_fc, &
                                 chi_fc,dx,dt,vel_bc_n)

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
    do n=1,nlevs
       call multifab_setval_c(s_update(n),0.d0,1,1,all=.true.)
       call multifab_copy_c(s_update(n),2,gmres_rhs_p(n),1,1,0)
    end do

    ! multiply gmres_rhs_p by -S_fac, so gmres_rhs_p = -S^n
    do n=1,nlevs
       call multifab_mult_mult_s_c(gmres_rhs_p(n),1,-S_fac,1,0)
    end do

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

    ! for inhomogeneous boundary conditions, convert problem to homogeneous by
    ! subtracting from the RHS the result of the operator applied to a solution
    ! vector with zeros everywhere in the problem domain, and ghost cells filled to
    ! respect the boundary conditions
    call convert_to_homogeneous(mla,gmres_rhs_v,gmres_rhs_p,s_fc,eta,eta_ed, &
                                kappa,1.d0/dt,dx,the_bc_tower, &
                                vel_bc_n_delta,vel_bc_t_delta)

    ! call gmres to compute delta v and delta p
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dp,s_fc, &
               eta,eta_ed,kappa,theta_fac)

    ! compute v^* = v^{n-1/2} + delta v
    ! compute p^* = p^{n-1/2} + delta p
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 3 - Forward-Euler Scalar Predictor
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (advection_type .ge. 1) then

       do n=1,nlevs
          ! AJN FIXME - ghost cells will stay set zero
          call multifab_copy_c(bds_force(n),1,s_update(n),1,nscal,0)
          call multifab_fill_boundary(bds_force(n))
       end do

       call bds(mla,umac,sold,s_update,bds_force,s_fc,dx,dt,1,nscal,the_bc_tower)

    else

       ! set s_update to A^n for scalars
       call mk_advective_s_fluxdiv(mla,umac,s_fc,s_update,dx,1,nscal)

    end if

    ! compute s^{*,n+1} = s^n + dt * (A^n + D^n + St^n)
    ! store result in snew (we will later add sold and divide by 2)
    do n=1,nlevs
       call multifab_mult_mult_s_c(s_update(n),1,dt,nscal,0)
       call multifab_copy_c(snew(n),1,sold(n),1,nscal,0)
       call multifab_plus_plus_c(snew(n),1,s_update(n),1,nscal,0)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 4 - Compute Midpoint Estimates
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! add sold to snew and multiply snew by 1/2 so it holds s^{*,n+1/2}
    do n=1,nlevs
       call multifab_plus_plus_c(snew(n),1,sold(n),1,nscal,0)
       call multifab_mult_mult_s_c(snew(n),1,0.5d0,nscal,0)
    end do

    ! convert s^{*,n+1/2} to prim
    call convert_cons_to_prim(mla,snew,prim,.true.)

    ! fill ghost cells for prim
    do n=1,nlevs
       call multifab_fill_boundary(prim(n))
       call multifab_physbc(prim(n),1,scal_bc_comp,2,the_bc_tower%bc_tower_array(n),dx(n,:))
    end do

    ! convert prim to s^{*,n+1/2} in valid and ghost region
    ! now s^{*,n+1/2} properly filled ghost cells
    call convert_cons_to_prim(mla,snew,prim,.false.)

    ! average s^{*,n+1/2} to faces
    do i=1,nscal
       call average_cc_to_face(nlevs,snew,s_fc,i,scal_bc_comp,1,the_bc_tower%bc_tower_array)
    end do

    ! compute (chi,eta,kappa)^{*,n+1/2}
    call compute_chi(mla,chi,chi_fc,prim,dx,the_bc_tower%bc_tower_array)
    call compute_eta(mla,eta,eta_ed,prim,dx,the_bc_tower%bc_tower_array)
    call compute_kappa(mla,kappa,prim,dx)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Steps 5 and 6 - Corrector Stokes Solve
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! build up rhs_v for gmres solve
    do n=1,nlevs
       do i=1,dm
          call setval(gmres_rhs_v(n,i),0.d0,all=.true.)
       end do
    end do

    ! compute grad p^*
    call compute_grad(mla,pnew,gradp,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

    ! subtract grad p^* from gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradp(n,i),1,1,0)
       end do
    end do

    ! add A_0^* v^* to gmres_rhs_v
    call mk_diffusive_m_fluxdiv(mla,gmres_rhs_v,umac,eta,eta_ed,kappa,dx, &
                                the_bc_tower%bc_tower_array)

    ! add div(Sigma^n') to gmres_rhs_v
    call mk_stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,gmres_rhs_v, &
                                 eta,eta_ed,dx,dt)

    ! add gravity term to gmres_rhs_v
    if (any(grav(1:dm) .ne. 0.d0)) then
       call mk_grav_force(mla,gmres_rhs_v,s_fc,s_fc)
    end if

    ! save boundary conditions from predictor
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(vel_bc_n_old(n,i),1,vel_bc_n(n,i),1,1,0)
       end do
       do i=1,size(vel_bc_t,dim=2)
          call multifab_copy_c(vel_bc_t_old(n,i),1,vel_bc_t(n,i),1,1,0)
       end do
    end do

    ! reset inhomogeneous bc condition to deal with reservoirs
    call set_inhomogeneous_vel_bcs(mla,vel_bc_n,vel_bc_t,eta_ed,dx, &
                                   the_bc_tower%bc_tower_array)

    ! add div(rho*chi grad c)^{*,n+1/2} to rhs_p
    call mk_diffusive_rhoc_fluxdiv(mla,gmres_rhs_p,1,prim,s_fc,chi_fc,dx, &
                                   the_bc_tower%bc_tower_array,vel_bc_n)

    ! add div(Psi^n') to rhs_p
    call mk_stochastic_s_fluxdiv(mla,the_bc_tower%bc_tower_array,gmres_rhs_p,s_fc, &
                                 chi_fc,dx,dt,vel_bc_n)

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
    ! then, set s_update for rho1 to F^{*,n+1/2} = div(rho*chi grad c)^{*,n+1/2} + div(Psi^n')
    do n=1,nlevs
       call multifab_setval_c(s_update(n),0.d0,1,1,all=.true.)
       call multifab_copy_c(s_update(n),2,gmres_rhs_p(n),1,1,0)
    end do

    ! multiply gmres_rhs_p by -S_fac, so gmres_rhs_p = -S^{*,n+1/2}
    do n=1,nlevs
       call multifab_mult_mult_s_c(gmres_rhs_p(n),1,-S_fac,1,0)
    end do

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

    ! for inhomogeneous boundary conditions, convert problem to homogeneous by
    ! subtracting from the RHS the result of the operator applied to a solution
    ! vector with zeros everywhere in the problem domain, and ghost cells filled to
    ! respect the boundary conditions
    call convert_to_homogeneous(mla,gmres_rhs_v,gmres_rhs_p,s_fc,eta,eta_ed, &
                                kappa,1.d0/dt,dx,the_bc_tower, &
                                vel_bc_n_delta,vel_bc_t_delta)

    ! call gmres to compute delta v and delta p
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dp,s_fc, &
               eta,eta_ed,kappa,theta_fac)

    ! compute v^{n+1/2} = v^* + delta v
    ! compute p^{n+1/2} = p^* + delta p
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 7 - Trapezoidal Scalar Corrector
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (advection_type .ge. 1) then

       do n=1,nlevs
          ! AJN FIXME - ghost cells will stay set zero
          call multifab_copy_c(bds_force(n),1,s_update(n),1,nscal,0)
          call multifab_fill_boundary(bds_force(n))
       end do

       call bds(mla,umac,sold,s_update,bds_force,s_fc,dx,dt,1,nscal,the_bc_tower)

    else

       ! set s_update to A^{*,n+1/2} for scalars
       call mk_advective_s_fluxdiv(mla,umac,s_fc,s_update,dx,1,nscal)

    end if

    ! compute s^{n+1} = s^n + dt * (A^{*,n+1/2} + D^{*,n+1/2} + St^{*,n+1/2})
    do n=1,nlevs
       call multifab_mult_mult_s_c(s_update(n),1,dt,nscal,0)
       call multifab_copy_c(snew(n),1,sold(n),1,nscal,0)
       call multifab_plus_plus_c(snew(n),1,s_update(n),1,nscal,0)
    end do

    ! convert s^{n+1} to prim
    call convert_cons_to_prim(mla,snew,prim,.true.)

    ! fill ghost cells for prim
    do n=1,nlevs
       call multifab_fill_boundary(prim(n))
       call multifab_physbc(prim(n),1,scal_bc_comp,2,the_bc_tower%bc_tower_array(n),dx(n,:))
    end do

    ! convert prim to s^{n+1} in valid and ghost region
    ! now s^{n+1} properly filled ghost cells
    call convert_cons_to_prim(mla,snew,prim,.false.)

    ! compute s^{n+1} to faces
    do i=1,nscal
       call average_cc_to_face(nlevs,snew,s_fc,i,scal_bc_comp,1,the_bc_tower%bc_tower_array)
    end do


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Compute stuff for plotfile and next time step
    
    ! AJN FIXME
    ! mnew should hold (1/2)(rho^n + rho^{n+1}) v^{*,n+1/2}
    ! this actually does mnew = rho^{n+1} v^{*,n+1/2}
    ! mnew is just a diagnostic, it does not enter the algorithm so for now this is fine
    call convert_m_to_umac(mla,s_fc,mnew,umac,.false.)

    ! compute (chi,eta,kappa)^{n+1}
    call compute_chi(mla,chi,chi_fc,prim,dx,the_bc_tower%bc_tower_array)
    call compute_eta(mla,eta,eta_ed,prim,dx,the_bc_tower%bc_tower_array)
    call compute_kappa(mla,kappa,prim,dx)

    ! fill the stochastic multifabs with a new set of random numbers
    call fill_stochastic(mla)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End Time-Advancement
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
          call multifab_destroy(vel_bc_n_old(n,i))
          call multifab_destroy(vel_bc_n_delta(n,i))
       end do
       do i=1,size(vel_bc_t,dim=2)
          call multifab_destroy(vel_bc_t_old(n,i))
          call multifab_destroy(vel_bc_t_delta(n,i))
       end do
    end do

    deallocate(vel_bc_t_old,vel_bc_t_delta)

  end subroutine advance_timestep_overdamped

end module advance_timestep_overdamped_module
