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
  use probin_lowmach_module, only: nscal, rhobar, diff_coef, visc_coef, grav
  use probin_common_module, only: fixed_dt

  use analysis_module

  implicit none

  private

  public :: advance_timestep_overdamped

contains

  subroutine advance_timestep_overdamped(mla,mnew,umac,sold,snew,s_fc,prim,pold,pnew, &
                                         chi,chi_fc,eta,eta_ed,kappa,dx,the_bc_tower, &
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

    theta_fac = 0.d0

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

    ! compute m_d_fluxdiv_old = A_0^n v^{n-1/2}
    call mk_diffusive_m_fluxdiv(mla,m_d_fluxdiv_old,umac_old,eta,eta_ed,kappa,dx, &
                                the_bc_tower%bc_tower_array)

    ! add A_0^n v^{n-1/2} to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_d_fluxdiv_old(n,i),1,1,0)
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
       call mk_grav_force(mla,gmres_rhs_v,s_fc_old,s_fc_old)
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
    do n=1,nlevs
       call multifab_setval_c(s_update(n),0.d0,1,1,all=.true.)
       call multifab_copy_c(s_update(n),2,gmres_rhs_p(n),1,1,0)
    end do

    ! multiply gmres_rhs_p by -S_fac, so gmres_rhs_p = -S^n
    do n=1,nlevs
       call multifab_mult_mult_s_c(gmres_rhs_p(n),1,-S_fac,1,0)
    end do

    ! compute div v^{n-1/2}
    call compute_div(mla,umac_old,divu,dx,1,1,1)

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
                                kappa,1.d0/fixed_dt,dx,the_bc_tower, &
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

    ! set s_update to A^n for scalars
    call mk_advective_s_fluxdiv(mla,umac,s_fc,s_update,dx,1,nscal)

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
    ! Step 4 - Compute Midpoint Estimates
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Steps 5 and 6 - Corrector Stokes Solve
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 7 - Trapezoidal Scalar Corrector
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    do n=1,nlevs
       call setval(gmres_rhs_p(n),0.d0,all=.true.)
    end do

    ! add div(rho*chi grad c)^{*,n+1} to rhs_p
    call mk_diffusive_rhoc_fluxdiv(mla,gmres_rhs_p,1,prim,s_fc,chi_fc,dx, &
                                   the_bc_tower%bc_tower_array,vel_bc_n)

    ! add div(Psi^n) to rhs_p
    call mk_stochastic_s_fluxdiv(mla,the_bc_tower%bc_tower_array,gmres_rhs_p,s_fc_old, &
                                 chi_fc_old,dx,vel_bc_n)

    ! reset s_update for all scalars to zero
    ! then, set s_update for rho1 to F^{*,n+1} = div(rho*chi grad c)^{*,n+1} + div(Psi^n)
    do n=1,nlevs
       call multifab_setval_c(s_update(n),0.d0,1,1,all=.true.)
       call multifab_copy_c(s_update(n),2,gmres_rhs_p(n),1,1,0)
    end do

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
    ! Compute stuff for plotfile and next time step
    
    ! mnew should hold (1/2)(rho^n + rho^{n+1}) v^*
    !
    !
    !

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

  end subroutine advance_timestep_overdamped

end module advance_timestep_overdamped_module
