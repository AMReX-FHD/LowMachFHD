module advance_timestep_imp_bousq_module

  use ml_layout_module
  use define_bc_module
  use bc_module
  use convert_stag_module
  use mk_advective_s_fluxdiv_module
  use diffusive_m_fluxdiv_module
  use stochastic_m_fluxdiv_module
  use stochastic_mass_fluxdiv_module
  use compute_mass_fluxdiv_charged_module
  use compute_HSE_pres_module
  use convert_m_to_umac_module
  use convert_rhoc_to_c_module
  use eos_check_module
  use mk_advective_m_fluxdiv_module
  use reservoir_bc_fill_module
  use gmres_module
  use div_and_grad_module
  use mk_grav_force_module
  use compute_mixture_properties_module
  use mass_flux_utilities_module
  use multifab_physbc_module
  use multifab_physbc_extrap_module
  use multifab_physbc_stag_module
  use zero_edgeval_module
  use fill_rho_ghost_cells_module
  use fluid_charge_module
  use ml_solve_module
  use bndry_reg_module
  use Epot_mass_fluxdiv_module
  use bl_rng_module
  use bl_random_module
  use probin_common_module, only: advection_type, grav, rhobar, variance_coef_mass, &
                                  variance_coef_mom, barodiffusion_type, project_eos_int, &
                                  use_bl_rng
  use probin_gmres_module, only: gmres_abs_tol, gmres_rel_tol, mg_verbose
  use probin_multispecies_module, only: nspecies
  use probin_charged_module, only: use_charged_fluid, dielectric_type, dielectric_const, &
                                   Epot_wall_bc_type

  use fabio_module

  implicit none

  private

  public :: advance_timestep_imp_bousq

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

  subroutine advance_timestep_imp_bousq(mla,umac,rho_old,rho_new,rhotot_old,rhotot_new, &
                                        gradp_baro,pi,eta,eta_ed,kappa,Temp,Temp_ed, &
                                        Epot_mass_fluxdiv, & ! just a placeholder, can get rid of _old and _tmp versions
                                        diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                        dx,dt,time,the_bc_tower,istep, &
                                        grad_Epot_old,grad_Epot_new,charge_old,charge_new, &
                                        Epot,permittivity,gradPhiApprox)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(inout) :: rho_old(:)
    type(multifab) , intent(inout) :: rho_new(:)
    type(multifab) , intent(inout) :: rhotot_old(:)
    type(multifab) , intent(inout) :: rhotot_new(:)
    type(multifab) , intent(inout) :: gradp_baro(:,:)
    type(multifab) , intent(inout) :: pi(:)
    ! eta and kappa need to enter consistent with old and leave consistent with new
    type(multifab) , intent(inout) :: eta(:)
    type(multifab) , intent(inout) :: eta_ed(:,:) ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(inout) :: kappa(:)
    type(multifab) , intent(inout) :: Temp(:)
    type(multifab) , intent(inout) :: Temp_ed(:,:) ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(inout) :: Epot_mass_fluxdiv(:)
    type(multifab) , intent(inout) :: diff_mass_fluxdiv(:)
    type(multifab) , intent(inout) :: stoch_mass_fluxdiv(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt,time
    type(bc_tower) , intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: istep
    type(multifab) , intent(inout) :: grad_Epot_old(:,:)
    type(multifab) , intent(inout) :: grad_Epot_new(:,:)
    type(multifab) , intent(inout) :: charge_old(:)
    type(multifab) , intent(inout) :: charge_new(:)
    type(multifab) , intent(inout) :: Epot(:)
    ! permittivity enters consistent with old and leaves consistent with new
    type(multifab) , intent(inout) :: permittivity(:)
    ! If present, this will contain 
    ! z^T*F/z^T*A = z^T*W*chi*Gamma*grad(x) / (rho^/(nkT)*z^T*W*chi*W*z)
    ! This is the ambipolar approximation to the gradient of the potential,
    ! which is wrong in general
    type(multifab) , intent(inout), optional :: gradPhiApprox(:,:) 

    ! local
    type(multifab) ::  diff_mass_fluxdiv_old(mla%nlevel)
    type(multifab) :: stoch_mass_fluxdiv_old(mla%nlevel)
    type(multifab) ::  Epot_mass_fluxdiv_old(mla%nlevel)
    type(multifab) ::  Epot_mass_fluxdiv_tmp(mla%nlevel)
    type(multifab) ::   adv_mass_fluxdiv_old(mla%nlevel)
    type(multifab) ::   adv_mass_fluxdiv    (mla%nlevel)

    type(multifab) ::              rho_tmp(mla%nlevel)
    type(multifab) ::           rhotot_tmp(mla%nlevel)
    type(multifab) ::          gmres_rhs_p(mla%nlevel)
    type(multifab) ::                  dpi(mla%nlevel)
    type(multifab) ::                 conc(mla%nlevel)
    type(multifab) ::           charge_np2(mla%nlevel)

    type(multifab) ::            mold(mla%nlevel,mla%dim)
    type(multifab) ::           mtemp(mla%nlevel,mla%dim)
    type(multifab) :: m_a_fluxdiv_new(mla%nlevel,mla%dim)
    type(multifab) :: m_a_fluxdiv_old(mla%nlevel,mla%dim)
    type(multifab) :: m_d_fluxdiv_new(mla%nlevel,mla%dim)
    type(multifab) :: m_d_fluxdiv_old(mla%nlevel,mla%dim)
    type(multifab) :: m_s_fluxdiv_new(mla%nlevel,mla%dim)
    type(multifab) :: m_s_fluxdiv_old(mla%nlevel,mla%dim)
    type(multifab) ::     gmres_rhs_v(mla%nlevel,mla%dim)
    type(multifab) ::           dumac(mla%nlevel,mla%dim)
    type(multifab) ::          gradpi(mla%nlevel,mla%dim)
    type(multifab) ::          rho_fc(mla%nlevel,mla%dim)
    type(multifab) ::       rhotot_fc(mla%nlevel,mla%dim)
    type(multifab) ::      flux_total(mla%nlevel,mla%dim)
    type(multifab) ::       flux_diff(mla%nlevel,mla%dim)
    type(multifab) ::   grad_Epot_np2(mla%nlevel,mla%dim)

    type(multifab) :: m_grav_force_old(mla%nlevel,mla%dim)
    type(multifab) :: m_grav_force_new(mla%nlevel,mla%dim)

    type(multifab) :: Lorentz_force_old(mla%nlevel,mla%dim)
    type(multifab) :: Lorentz_force_new(mla%nlevel,mla%dim)
    type(multifab) :: Lorentz_force_np2(mla%nlevel,mla%dim)

    type(multifab) ::     solver_alpha(mla%nlevel)         ! alpha=0 for Poisson solve
    type(multifab) ::       solver_rhs(mla%nlevel)         ! Poisson solve rhs
    type(multifab) ::            A_Phi(mla%nlevel,mla%dim) ! face-centered A_Phi
    type(multifab) ::      solver_beta(mla%nlevel,mla%dim) ! beta=epsilon+dt*z^T*A_Phi for Poisson solve
    type(multifab) ::  permittivity_fc(mla%nlevel,mla%dim) ! beta=epsilon+dt*z^T*A_Phi for Poisson solve

    type(multifab) :: zdotA(mla%nlevel,mla%dim)
    
    type(bndry_reg) :: fine_flx(mla%nlevel)

    integer :: i,dm,n,nlevs,comp

    real(kind=dp_t) :: theta_alpha, norm_pre_rhs, gmres_abs_tol_in

    real(kind=dp_t) :: weights(1), sum

    weights(1) = 1.d0

    if (any(rhobar(1:nspecies) .ne. rhobar(1))) then
       call bl_error("Implicit Boussinesq algorithm requires all the same rhobar's")
    end if

    nlevs = mla%nlevel
    dm = mla%dim

    theta_alpha = 1.d0/dt
    
    call build_bc_multifabs(mla)
    
    do n=1,nlevs
       call multifab_build( diff_mass_fluxdiv_old(n),mla%la(n),nspecies,0) 
       call multifab_build(stoch_mass_fluxdiv_old(n),mla%la(n),nspecies,0) 
       call multifab_build( Epot_mass_fluxdiv_old(n),mla%la(n),nspecies,0) 
       call multifab_build( Epot_mass_fluxdiv_tmp(n),mla%la(n),nspecies,0) 
       call multifab_build(adv_mass_fluxdiv_old(n),mla%la(n),nspecies,0)
       call multifab_build(adv_mass_fluxdiv    (n),mla%la(n),nspecies,0)

       call multifab_build(    rho_tmp(n),mla%la(n),nspecies,rho_old(n)%ng)
       call multifab_build( rhotot_tmp(n),mla%la(n),1       ,rhotot_old(n)%ng)
       call multifab_build(gmres_rhs_p(n),mla%la(n),1       ,0)
       call multifab_build(        dpi(n),mla%la(n),1       ,1)
       call multifab_build(       conc(n),mla%la(n),nspecies,rho_old(n)%ng)
       call multifab_build( charge_np2(n),mla%la(n),1       ,1)
       do i=1,dm
          call multifab_build_edge(              mold(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(             mtemp(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(   m_a_fluxdiv_new(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(   m_a_fluxdiv_old(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(   m_d_fluxdiv_new(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(   m_d_fluxdiv_old(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(   m_s_fluxdiv_new(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(   m_s_fluxdiv_old(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(       gmres_rhs_v(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(             dumac(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(            gradpi(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(         rhotot_fc(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(            rho_fc(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge(        flux_total(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge(         flux_diff(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge(     grad_Epot_np2(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(  m_grav_force_old(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(  m_grav_force_new(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge( Lorentz_force_old(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge( Lorentz_force_new(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge( Lorentz_force_np2(n,i),mla%la(n),1       ,0,i)
       end do
       call multifab_build(solver_alpha(n),mla%la(n),1,0)
       call multifab_build(solver_rhs(n),mla%la(n),1,0)
       do i=1,dm
          call multifab_build_edge(A_Phi(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge(solver_beta(n,i),mla%la(n),1,0,i)
          call multifab_build_edge(permittivity_fc(n,i),mla%la(n),1,0,i)
          call multifab_build_edge(zdotA(n,i),mla%la(n),1,0,i)
       end do
       call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
    end do

    do n=1,nlevs
       ! alpha=0
       ! this is used for the electric potential Poisson solve
       ! and for the evaluation of div A_Phi grad Phi
       call multifab_setval(solver_alpha(n),0.d0)
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
    do n=1,nlevs
       call multifab_setval(adv_mass_fluxdiv_old(n),0.d0)
    end do
    call mk_advective_s_fluxdiv(mla,umac,rho_fc,adv_mass_fluxdiv_old,dx,1,nspecies)

    ! compute m_d_fluxdiv_old = L_0^n v^n
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(m_d_fluxdiv_old(n,i),0.d0,all=.true.)
       end do
    end do
    call diffusive_m_fluxdiv(mla,m_d_fluxdiv_old,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    ! compute m_s_fluxdiv_old = div (sqrt() (W + W^T)^{n:n+1})
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(m_s_fluxdiv_old(n,i),0.d0,all=.true.)
       end do
    end do
    if (variance_coef_mom .ne. 0.d0) then   
       ! fill the stochastic momentum multifabs with new sets of random numbers
       call fill_m_stochastic(mla)
       call stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,m_s_fluxdiv_old, &
                                 eta,eta_ed,Temp,Temp_ed,dx,dt,weights)
    end if

    ! compute "old" gravity force
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(m_grav_force_old(n,i),0.d0,all=.true.)
       end do
    end do
    if (any(grav(1:dm) .ne. 0.d0)) then
       call mk_grav_force(mla,m_grav_force_old,rhotot_fc,rhotot_fc,the_bc_tower)
    end if

    ! compute "old" Lorentz force (NOT USED RIGHT NOW)
    call compute_Lorentz_force(mla,Lorentz_force_old,grad_Epot_old,permittivity, &
                               charge_old,dx,the_bc_tower)

    ! compute momentum, mold = rho^n v^n
    call convert_m_to_umac(mla,rhotot_fc,mold,umac,.false.)

    ! compute advective flux divergence, m_a_fluxdiv_old = div(-rho v v)^n
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(m_a_fluxdiv_old(n,i),0.d0,all=.true.)
       end do
    end do
    call mk_advective_m_fluxdiv(mla,umac,mold,m_a_fluxdiv_old,dx,the_bc_tower%bc_tower_array)

    ! compute A_Phi^n for Poisson solve (does not have z^T)
    call implicit_potential_coef(mla,rho_old,Temp,A_Phi,the_bc_tower)

    ! build RHS for Poisson solve (store in rho_tmp)
    ! (rho w)^n - dt div (F_a^n + F_d^n)
    do n=1,nlevs
       call multifab_copy_c(rho_tmp(n),1,rho_old(n),1,nspecies,0)
       call multifab_saxpy_3_cc(rho_tmp(n),1,dt, adv_mass_fluxdiv_old(n),1,nspecies)
       call multifab_saxpy_3_cc(rho_tmp(n),1,dt,diff_mass_fluxdiv_old(n),1,nspecies)
       if (variance_coef_mass .ne. 0.d0) then
          call multifab_saxpy_3_cc(rho_tmp(n),1,dt,stoch_mass_fluxdiv_old(n),1,nspecies)
       end if
    end do

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

    ! right-hand-side for Poisson solve is z^T [ (rho w)^n - dt div (F_a^n + F_d^n) ]
    call dot_with_z(mla,rho_tmp,solver_rhs)

    ! compute z^T A_Phi^n, store in solver_beta
    call dot_with_z_face(mla,A_Phi,solver_beta)

    ! permittivity on faces
    call average_cc_to_face(nlevs,permittivity,permittivity_fc,1,scal_bc_comp,1, &
                            the_bc_tower%bc_tower_array)

    ! compute solver_beta = epsilon + dt theta z^T A_Phi^n
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_s_c(solver_beta(n,i),1,dt,1,0)
          call multifab_plus_plus_c(solver_beta(n,i),1,permittivity_fc(n,i),1,1,0)
       end do
    end do

    ! initial guess for Phi
    do n=1,nlevs
       call multifab_setval(Epot(n),0.d0,all=.true.)
       ! fill ghost cells for Epot at walls using Dirichlet value
       call multifab_physbc(Epot(n),1,Epot_bc_comp,1,the_bc_tower%bc_tower_array(n), &
                            dx_in=dx(n,:))
    end do

    ! for inhomogeneous Neumann bc's for electric potential, put in homogeneous form
    if (Epot_wall_bc_type .eq. 2) then
       call inhomogeneous_neumann_fix(mla,solver_rhs,permittivity,dx,the_bc_tower)
    end if

    ! solve -div (epsilon + dt theta z^T A_Phi^n) grad Phi^{n+1,*} = z^T RHS
    call ml_cc_solve(mla,solver_rhs,Epot,fine_flx,solver_alpha,solver_beta,dx, &
                     the_bc_tower,Epot_bc_comp,verbose=mg_verbose)

    ! for periodic problems subtract off the average of Epot
    ! we can generalize this later for walls
    if ( all(mla%pmask(1:dm)) ) then
       sum = multifab_sum(Epot(1)) / multifab_volume(Epot(1))
       call multifab_sub_sub_s(Epot(1),sum)
       do n=1,nlevs
          call multifab_physbc(Epot(n),1,Epot_bc_comp,1,the_bc_tower%bc_tower_array(n), &
                               dx_in=dx(n,:))
          call multifab_fill_boundary(Epot(n))
       end do
    end if
       
    ! compute grad Phi^{n+1,*} for use in Lorentz force
    call compute_grad(mla,Epot,grad_Epot_new,dx,1,Epot_bc_comp,1,1,the_bc_tower%bc_tower_array)

    do comp=1,nspecies

       ! copy component of A_Phi^n into solver_beta and multiply by grad_Epot_new
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(solver_beta(n,i),1,A_Phi(n,i),comp,1,0)
             call multifab_mult_mult_c(solver_beta(n,i),1,grad_Epot_new(n,i),1,1,0)
          end do
          ! zero mass flux on walls
          call zero_edgeval_walls(solver_beta(n,:),1,1,the_bc_tower%bc_tower_array(n))
       end do

       ! compute Epot_mass_fluxdiv_old = div A_Phi^n grad Epot^{n+1,*}
       call compute_div(mla,solver_beta,Epot_mass_fluxdiv_old,dx,1,comp,1)
    end do

    ! add dt*Epot_mass_fluxdiv_old to RHS to get rho_tmp = rho^{n+1,*}
    do n=1,nlevs
       call multifab_saxpy_3_cc(rho_tmp(n),1,dt,Epot_mass_fluxdiv_old(n),1,nspecies)
    end do

    ! add (dt/2)*Epot_mass_fluxdiv_old to continue building rho_new
    do n=1,nlevs
       call multifab_saxpy_3_cc(rho_new(n),1,0.5d0*dt,Epot_mass_fluxdiv_old(n),1,nspecies)
    end do

    ! compute rhotot^{n+1,*} from rho^{n+1,*} in VALID REGION
    call compute_rhotot(mla,rho_tmp,rhotot_tmp)

    ! rho to conc - NO GHOST CELLS
    call convert_rhoc_to_c(mla,rho_tmp,rhotot_tmp,conc,.true.)
    call fill_c_ghost_cells(mla,conc,dx,the_bc_tower)

    ! fill conc ghost cells
    do n=1,nlevs
       call fill_rho_ghost_cells(conc(n),rhotot_tmp(n),the_bc_tower%bc_tower_array(n))
    end do

    ! conc to rho - INCLUDING GHOST CELLS
    call convert_rhoc_to_c(mla,rho_tmp,rhotot_tmp,conc,.false.)

    ! print out EOS drift
    call eos_check(mla,rho_tmp)

    ! average (rho,rhotot)^{n+1,*} to faces
    call average_cc_to_face(nlevs,   rho_tmp,   rho_fc,1,   c_bc_comp,nspecies,the_bc_tower%bc_tower_array)
    call average_cc_to_face(nlevs,rhotot_tmp,rhotot_fc,1,scal_bc_comp,       1,the_bc_tower%bc_tower_array)

    ! compute charge^{n+1,*}
    call dot_with_z(mla,rho_tmp,charge_new)

    ! compute permittivity^{n+1,*}
    if (dielectric_type .ne. 0) then
       call compute_permittivity(mla,permittivity,rho_tmp,rhotot_tmp,the_bc_tower)
    end if

    ! compute rho^{n+1,*}*g
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(m_grav_force_new(n,i),0.d0,all=.true.)
       end do
    end do
    if (any(grav(1:dm) .ne. 0.d0)) then
       call mk_grav_force(mla,m_grav_force_new,rhotot_fc,rhotot_fc,the_bc_tower)
    end if

    ! compute (eta,kappa)^{n+1,*}
    call compute_eta_kappa(mla,eta,eta_ed,kappa,rho_tmp,rhotot_tmp,Temp,dx, &
                           the_bc_tower%bc_tower_array)

    ! compute "new" Lorentz_force^{n+1,*}
    call compute_Lorentz_force(mla,Lorentz_force_new,grad_Epot_new,permittivity, &
                               charge_new,dx,the_bc_tower)

    ! diff_mass_fluxdiv = -div (rho W chi Gamma grad x + ...)^{n+1,*}
    ! stoch_mass_fluxdiv = -div sqrt(...) B^{n+1,*} Z^{n:n+1}
    ! and flux_total for reservoir boundary conditions on velocity
    call compute_mass_fluxdiv_charged(mla,rho_tmp,gradp_baro, &
                                      diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                      Temp,flux_total,flux_diff, &
                                      dt,time,dx,weights,the_bc_tower)

    ! now fluxdivs contain "-div(F) = div (rho W chi Gamma grad x + ...)", etc.
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

    if (barodiffusion_type .ne. 0) then
       call bl_error("advance_timestep_imp_bousq: barodiffusion not supported yet")
    end if
       
    ! subtract grad pi^n from gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradpi(n,i),1,1,0)
       end do
    end do

    ! add advective momentum flux divergence, m_a_fluxdiv_old = div(-rho v v)^n
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,1.d0,m_a_fluxdiv_old(n,i),1,1)
       end do
    end do

    ! compute m_d_fluxdiv_new = L_0^{n+1,*} vbar, where vbar = v^n with t^{n+1} bc's
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(m_d_fluxdiv_new(n,i),0.d0,all=.true.)
       end do
    end do
    call diffusive_m_fluxdiv(mla,m_d_fluxdiv_new,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    ! add (1/2) (L_0^n v^n + L_0^{n+1,*} vbar^n)
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_d_fluxdiv_old(n,i),1,1)
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_d_fluxdiv_new(n,i),1,1)
       end do
    end do

    ! compute m_s_fluxdiv_new = div (sqrt() (W + W^T)^{n:n+1})
    ! (these should only differ from the t^n stochastic fluxdiv because of eta^{n+1,*})
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(m_s_fluxdiv_new(n,i),0.d0,all=.true.)
       end do
    end do
    if (variance_coef_mom .ne. 0.d0) then
       call stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,m_s_fluxdiv_new, &
                                 eta,eta_ed,Temp,Temp_ed,dx,dt,weights)
    end if

    ! add (1/2) (t^n + t^{n+1,*}) stochastic fluxes
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_s_fluxdiv_old(n,i),1,1)
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_s_fluxdiv_new(n,i),1,1)
       end do
    end do

    ! add (1/2) (t^n + t^{n+1,*}) gravitational force
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_grav_force_old(n,i),1,1)
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_grav_force_new(n,i),1,1)
       end do
    end do

    ! subtract t^{n+1,*} Lorentz force
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,-1.d0,Lorentz_force_new(n,i),1,1)
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
    do n=1,nlevs
       call multifab_setval(adv_mass_fluxdiv(n),0.d0)
    end do
    call mk_advective_s_fluxdiv(mla,umac,rho_fc,adv_mass_fluxdiv,dx,1,nspecies)

    ! compute advective flux divergence, m_a_fluxdiv_new = div(-rho v v)^{n+1,*}
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(m_a_fluxdiv_new(n,i),0.d0,all=.true.)
       end do
    end do
    call mk_advective_m_fluxdiv(mla,umac,mtemp,m_a_fluxdiv_new,dx,the_bc_tower%bc_tower_array)

    ! compute A_Phi^{n+1,*} for Poisson solve (does not have z^T)
    call implicit_potential_coef(mla,rho_tmp,Temp,A_Phi,the_bc_tower)

    ! build RHS for Poisson solve (store in rho_tmp)
    ! (rho w)^{n+1,*} - dt div (F_a^{n+1,*} + F_d^{n+1,*})
    do n=1,nlevs
       call multifab_saxpy_3_cc(rho_tmp(n),1,dt,adv_mass_fluxdiv(n),1,nspecies)
       call multifab_saxpy_3_cc(rho_tmp(n),1,dt,diff_mass_fluxdiv(n),1,nspecies)
       if (variance_coef_mass .ne. 0.d0) then
          call multifab_saxpy_3_cc(rho_tmp(n),1,dt,stoch_mass_fluxdiv(n),1,nspecies)
       end if
    end do

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

    ! right-hand-side for Poisson solve is z^T [ (rho w)^{n+1,*} - dt div (F_a^{n+1,*} + F_d^{n+1,*}) ]
    call dot_with_z(mla,rho_tmp,solver_rhs)

    ! compute z^T A_Phi^{n+1,*}, store in solver_beta
    call dot_with_z_face(mla,A_Phi,solver_beta)

    ! permittivity^{n+1,*} on faces
    call average_cc_to_face(nlevs,permittivity,permittivity_fc,1,scal_bc_comp,1, &
                            the_bc_tower%bc_tower_array)

    ! compute solver_beta = epsilon + dt theta z^T A_Phi^{n+1,*}
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_s_c(solver_beta(n,i),1,dt,1,0)
          call multifab_plus_plus_c(solver_beta(n,i),1,permittivity_fc(n,i),1,1,0)
       end do
    end do

    ! initial guess for Phi
    do n=1,nlevs
       call multifab_setval(Epot(n),0.d0,all=.true.)
       ! fill ghost cells for Epot at walls using Dirichlet value
       call multifab_physbc(Epot(n),1,Epot_bc_comp,1,the_bc_tower%bc_tower_array(n), &
                            dx_in=dx(n,:))
    end do

    ! for inhomogeneous Neumann bc's for electric potential, put in homogeneous form
    if (Epot_wall_bc_type .eq. 2) then
       call inhomogeneous_neumann_fix(mla,solver_rhs,permittivity,dx,the_bc_tower)
    end if

    ! solve -div (epsilon + dt theta z^T A_Phi^{n+1,*}) grad Phi^{n+2,*} = z^T RHS
    call ml_cc_solve(mla,solver_rhs,Epot,fine_flx,solver_alpha,solver_beta,dx, &
                     the_bc_tower,Epot_bc_comp,verbose=mg_verbose)

    ! for periodic problems subtract off the average of Epot
    ! we can generalize this later for walls
    if ( all(mla%pmask(1:dm)) ) then
       sum = multifab_sum(Epot(1)) / multifab_volume(Epot(1))
       call multifab_sub_sub_s(Epot(1),sum)
       do n=1,nlevs
          call multifab_physbc(Epot(n),1,Epot_bc_comp,1,the_bc_tower%bc_tower_array(n), &
                               dx_in=dx(n,:))
          call multifab_fill_boundary(Epot(n))
       end do
    end if
       
    ! compute grad Phi^{n+2,*} for use in Lorentz force
    call compute_grad(mla,Epot,grad_Epot_np2,dx,1,Epot_bc_comp,1,1,the_bc_tower%bc_tower_array)

    do comp=1,nspecies

       ! copy component of A_Phi^{n+1,*} into solver_beta and multiply by grad_Epot_np2
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(solver_beta(n,i),1,A_Phi(n,i),comp,1,0)
             call multifab_mult_mult_c(solver_beta(n,i),1,grad_Epot_np2(n,i),1,1,0)
          end do
          ! zero mass flux on walls
          call zero_edgeval_walls(solver_beta(n,:),1,1,the_bc_tower%bc_tower_array(n))
       end do

       ! compute Epot_mass_fluxdiv_tmp = div A_Phi^{n+1,*} grad Epot^{n+2,*}
       call compute_div(mla,solver_beta,Epot_mass_fluxdiv_tmp,dx,1,comp,1)
    end do

    ! add dt*Epot_mass_fluxdiv_tmp to RHS to get rho_tmp = rho^{n+2,*}
    do n=1,nlevs
       call multifab_saxpy_3_cc(rho_tmp(n),1,dt,Epot_mass_fluxdiv_tmp(n),1,nspecies)
    end do

    ! compute rhotot^{n+2,*} from rho^{n+2,*} in VALID REGION
    call compute_rhotot(mla,rho_tmp,rhotot_tmp)

    ! rho to conc - NO GHOST CELLS
    call convert_rhoc_to_c(mla,rho_tmp,rhotot_tmp,conc,.true.)
    call fill_c_ghost_cells(mla,conc,dx,the_bc_tower)

    ! fill conc ghost cells
    do n=1,nlevs
       call fill_rho_ghost_cells(conc(n),rhotot_tmp(n),the_bc_tower%bc_tower_array(n))
    end do

    ! conc to rho - INCLUDING GHOST CELLS
    call convert_rhoc_to_c(mla,rho_tmp,rhotot_tmp,conc,.false.)

    ! compute charge^{n+2,*}
    call dot_with_z(mla,rho_tmp,charge_np2)

    ! compute permittivity^{n+2,*}
    if (dielectric_type .ne. 0) then
       call compute_permittivity(mla,permittivity,rho_tmp,rhotot_tmp,the_bc_tower)
    end if

    ! compute Lorentz_force^{n+2,*} (NOT USED RIGHT NOW)
    call compute_Lorentz_force(mla,Lorentz_force_np2,grad_Epot_np2,permittivity, &
                               charge_np2,dx,the_bc_tower)

    ! add (dt/2)*Epot_mass_fluxdiv_tmp to finishg building rho_new = rho^{n+1}
    do n=1,nlevs
       call multifab_saxpy_3_cc(rho_new(n),1,0.5d0*dt,Epot_mass_fluxdiv_tmp(n),1,nspecies)
    end do

    ! compute rhotot^{n+1} from rho^{n+1} in VALID REGION
    call compute_rhotot(mla,rho_new,rhotot_new)

    ! rho to conc - NO GHOST CELLS
    call convert_rhoc_to_c(mla,rho_new,rhotot_new,conc,.true.)
    call fill_c_ghost_cells(mla,conc,dx,the_bc_tower)

    ! fill conc ghost cells
    do n=1,nlevs
       call fill_rho_ghost_cells(conc(n),rhotot_new(n),the_bc_tower%bc_tower_array(n))
    end do

    ! conc to rho - INCLUDING GHOST CELLS
    call convert_rhoc_to_c(mla,rho_new,rhotot_new,conc,.false.)

    ! print out EOS drift
    call eos_check(mla,rho_new)

    ! average (rho,rhotot)^{n+1} to faces
    call average_cc_to_face(nlevs,   rho_new,   rho_fc,1,   c_bc_comp,nspecies,the_bc_tower%bc_tower_array)
    call average_cc_to_face(nlevs,rhotot_new,rhotot_fc,1,scal_bc_comp,       1,the_bc_tower%bc_tower_array)

    ! compute charge^{n+1}
    call dot_with_z(mla,rho_new,charge_new)

    ! compute permittivity^{n+1}
    if (dielectric_type .ne. 0) then
       call compute_permittivity(mla,permittivity,rho_new,rhotot_new,the_bc_tower)
    end if

    ! compute rho^{n+1}*g
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(m_grav_force_new(n,i),0.d0,all=.true.)
       end do
    end do
    if (any(grav(1:dm) .ne. 0.d0)) then
       call mk_grav_force(mla,m_grav_force_new,rhotot_fc,rhotot_fc,the_bc_tower)
    end if

    ! compute (eta,kappa)^{n+1}
    call compute_eta_kappa(mla,eta,eta_ed,kappa,rho_new,rhotot_new,Temp,dx, &
                           the_bc_tower%bc_tower_array)

    ! AJN HACK - compute a new grad Epot^{n+1} here?
    ! this will become grad Epot^n at the beginning of the next time step
!    if (dielectric_const .ne. 0.d0) then
!       call average_cc_to_face(nlevs,permittivity,permittivity_fc,1,scal_bc_comp,1, &
!                               the_bc_tower%bc_tower_array)
!       do n=1,nlevs
!          call multifab_setval(Epot(n),0.d0,all=.true.)
!          call multifab_copy_c(solver_rhs(n),1,charge_new(n),1,1,0)
!          do i=1,dm
!             call multifab_copy_c(solver_beta(n,i),1,permittivity_fc(n,i),1,1,0)
!          end do
!       end do
!       call ml_cc_solve(mla,solver_rhs,Epot,fine_flx,solver_alpha,solver_beta,dx, &
!                        the_bc_tower,Epot_bc_comp,verbose=mg_verbose)
!       call compute_grad(mla,Epot,grad_Epot_new,dx,1,Epot_bc_comp,1,1,the_bc_tower%bc_tower_array)
!       ! compute "new" Lorentz_force^{n+1}
!       call compute_Lorentz_force(mla,Lorentz_force_new,grad_Epot_new,permittivity, &
!                                  charge_new,dx,the_bc_tower)
!    else
!       do n=1,nlevs
!          do i=1,dm
!             call multifab_copy_c(Lorentz_force_new(n,i),1,Lorentz_force_np2(n,i),1,1,0)
!          end do
!       end do
!    end if

    ! fill the stochastic mass multifabs with new sets of random numbers
    if (variance_coef_mass .ne. 0.d0) then
       if (use_bl_rng) then
          call bl_rng_copy_engine(rng_eng_diffusion_old,rng_eng_diffusion)
       end if
       call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)
    end if

    ! mass fluxes for the next time step
    ! diff_mass_fluxdiv = -div (rho W chi Gamma grad x + ...)^{n+1}
    ! stoch_mass_fluxdiv = -div sqrt(...) B^{n+1} Z^{n+1:n+2}
    ! and flux_total for reservoir boundary conditions on velocity
    call compute_mass_fluxdiv_charged(mla,rho_new,gradp_baro, &
                                      diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                      Temp,flux_total,flux_diff, &
                                      dt,time,dx,weights,the_bc_tower)

    ! now fluxes contain "-div(F) = div (rho W chi Gamma grad x + ...)", etc.
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
       
    if(present(gradPhiApprox)) then 
       ! compute grad(phi) approximation, z^T Fmass / z^T A_Phi
       call dot_with_z_face(mla,flux_diff,gradPhiApprox)
       call dot_with_z_face(mla,A_Phi,zdotA)
       do n=1,nlevs
          do i=1,dm
             call multifab_div_div_c(gradPhiApprox(n,i),1,zdotA(n,i),1,1,0)
          end do
       end do
    end if

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

    if (barodiffusion_type .ne. 0) then
       call bl_error("advance_timestep_imp_bousq: barodiffusion not supported yet")
    end if
       
    ! subtract grad pi^{n+1,*} from gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradpi(n,i),1,1,0)
       end do
    end do

    ! add advective momentum flux divergence, (1/2) (div(-rho v v)^n + div(-rho v v)^{n+1,*})
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_a_fluxdiv_old(n,i),1,1)
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_a_fluxdiv_new(n,i),1,1)
       end do
    end do

    ! compute m_d_fluxdiv_new = L_0^{n+1} vbar, where vbar = v^{n+1,*} with t^{n+1} bc's
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(m_d_fluxdiv_new(n,i),0.d0,all=.true.)
       end do
    end do
    call diffusive_m_fluxdiv(mla,m_d_fluxdiv_new,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    ! add (1/2) (L_0^n v^n + L_0^{n+1} vbar^{n+1,*})
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_d_fluxdiv_old(n,i),1,1)
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_d_fluxdiv_new(n,i),1,1)
       end do
    end do

    ! compute m_s_fluxdiv_new = div (sqrt() (W + W^T)^{n:n+1})
    ! (these should only differ from the t^n stochastic fluxdiv because of eta^{n+1})
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(m_s_fluxdiv_new(n,i),0.d0,all=.true.)
       end do
    end do
    if (variance_coef_mom .ne. 0.d0) then
       call stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,m_s_fluxdiv_new, &
                                 eta,eta_ed,Temp,Temp_ed,dx,dt,weights)
    end if

    ! add (1/2) (t^n + t^{n+1}) stochastic fluxes
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_s_fluxdiv_old(n,i),1,1)
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_s_fluxdiv_new(n,i),1,1)
       end do
    end do

    ! add (1/2) (t^n + t^{n+1}) gravitational force
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_grav_force_old(n,i),1,1)
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_grav_force_new(n,i),1,1)
       end do
    end do

    ! subtract t^{n+1,*} Lorentz force
    do n=1,nlevs
       do i=1,dm
         call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,-1.d0,Lorentz_force_new(n,i),1,1)
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

    call destroy_bc_multifabs(mla)

    do n=1,nlevs
       call multifab_destroy(diff_mass_fluxdiv_old(n))
       call multifab_destroy(stoch_mass_fluxdiv_old(n))
       call multifab_destroy(Epot_mass_fluxdiv_old(n))
       call multifab_destroy(Epot_mass_fluxdiv_tmp(n))
       call multifab_destroy(adv_mass_fluxdiv_old(n))
       call multifab_destroy(adv_mass_fluxdiv(n))
       call multifab_destroy(rho_tmp(n))
       call multifab_destroy(rhotot_tmp(n))
       call multifab_destroy(gmres_rhs_p(n))
       call multifab_destroy(dpi(n))
       call multifab_destroy(conc(n))
       call multifab_destroy(charge_np2(n))
       do i=1,dm
          call multifab_destroy(mold(n,i))
          call multifab_destroy(mtemp(n,i))
          call multifab_destroy(m_a_fluxdiv_new(n,i))
          call multifab_destroy(m_a_fluxdiv_old(n,i))
          call multifab_destroy(m_d_fluxdiv_new(n,i))
          call multifab_destroy(m_d_fluxdiv_old(n,i))
          call multifab_destroy(m_s_fluxdiv_new(n,i))
          call multifab_destroy(m_s_fluxdiv_old(n,i))
          call multifab_destroy(gmres_rhs_v(n,i))
          call multifab_destroy(dumac(n,i))
          call multifab_destroy(rhotot_fc(n,i))
          call multifab_destroy(gradpi(n,i))
          call multifab_destroy(rho_fc(n,i))
          call multifab_destroy(flux_total(n,i))
          call multifab_destroy(flux_diff(n,i))
          call multifab_destroy(grad_Epot_np2(n,i))
          call multifab_destroy(m_grav_force_old(n,i))
          call multifab_destroy(m_grav_force_new(n,i))
          call multifab_destroy(Lorentz_force_old(n,i))
          call multifab_destroy(Lorentz_force_new(n,i))
          call multifab_destroy(Lorentz_force_np2(n,i))
       end do
       call multifab_destroy(solver_alpha(n))
       call multifab_destroy(solver_rhs(n))
       do i=1,dm
          call multifab_destroy(A_Phi(n,i))
          call multifab_destroy(solver_beta(n,i))
          call multifab_destroy(permittivity_fc(n,i))
          call multifab_destroy(zdotA(n,i))
       end do
       call bndry_reg_destroy(fine_flx(n))
    end do

  end subroutine advance_timestep_imp_bousq

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

end module advance_timestep_imp_bousq_module
