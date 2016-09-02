module advance_timestep_iterative_module

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
  use probin_common_module, only: advection_type, grav, rhobar, variance_coef_mass, &
                                  variance_coef_mom, barodiffusion_type, project_eos_int
  use probin_gmres_module, only: gmres_abs_tol, gmres_rel_tol, mg_verbose
  use probin_multispecies_module, only: nspecies
  use probin_charged_module, only: use_charged_fluid, theta_pot, dielectric_type, &
                                   num_pot_iters, dpdt_factor, Epot_wall_bc_type

  implicit none

  private

  public :: advance_timestep_iterative

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

  subroutine advance_timestep_iterative(mla,umac,rho_old,rho_new,rhotot_old,rhotot_new, &
                                        gradp_baro,pi,eta,eta_ed,kappa,Temp,Temp_ed, &
                                        Epot_mass_fluxdiv, &
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
    ! If present, this will contain z^T*F/z^T*A = z^T*W*chi*Gamma*grad(x) / (rho^/(nkT)*z^T*W*chi*W*z)
    ! This is the ambipolar approximation to the gradient of the potential, which is wrong in general
    type(multifab) , intent(inout), optional :: gradPhiApprox(:,:) 

    ! local
    type(multifab) ::  diff_mass_fluxdiv_old(mla%nlevel)
    type(multifab) :: stoch_mass_fluxdiv_old(mla%nlevel)
    type(multifab) ::  Epot_mass_fluxdiv_old(mla%nlevel)

    type(multifab) :: adv_mass_fluxdiv_old(mla%nlevel)
    type(multifab) :: adv_mass_fluxdiv_new(mla%nlevel)
    type(multifab) ::          gmres_rhs_p(mla%nlevel)
    type(multifab) ::                  dpi(mla%nlevel)
    type(multifab) ::                 divu(mla%nlevel)
    type(multifab) ::                 conc(mla%nlevel)
    type(multifab) ::                S_inc(mla%nlevel)

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

    type(multifab) :: m_grav_force_old(mla%nlevel,mla%dim)
    type(multifab) :: m_grav_force_new(mla%nlevel,mla%dim)

    type(multifab) :: Lorentz_force_old(mla%nlevel,mla%dim)
    type(multifab) :: Lorentz_force_new(mla%nlevel,mla%dim)

    type(multifab) ::     solver_alpha(mla%nlevel)         ! alpha=0 for Poisson solve
    type(multifab) ::       solver_rhs(mla%nlevel)         ! Poisson solve rhs
    type(multifab) ::            A_Phi(mla%nlevel,mla%dim) ! face-centered A_Phi
    type(multifab) ::      solver_beta(mla%nlevel,mla%dim) ! beta=epsilon+dt*z^T*A_Phi for Poisson solve
    type(multifab) ::  permittivity_fc(mla%nlevel,mla%dim) ! beta=epsilon+dt*z^T*A_Phi for Poisson solve

    type(multifab) :: zdotA(mla%nlevel,mla%dim)
    
    type(bndry_reg) :: fine_flx(mla%nlevel)

    integer :: i,dm,n,nlevs,comp,l

    real(kind=dp_t) :: theta_alpha, norm_pre_rhs, gmres_abs_tol_in

    real(kind=dp_t) :: weights(1), sum

    weights(1) = 1.d0

    nlevs = mla%nlevel
    dm = mla%dim

    theta_alpha = 1.d0/dt
    
    call build_bc_multifabs(mla)
    
    do n=1,nlevs
       call multifab_build( diff_mass_fluxdiv_old(n),mla%la(n),nspecies,0) 
       call multifab_build(stoch_mass_fluxdiv_old(n),mla%la(n),nspecies,0) 
       call multifab_build( Epot_mass_fluxdiv_old(n),mla%la(n),nspecies,0) 

       call multifab_build(adv_mass_fluxdiv_old(n),mla%la(n),nspecies,0)
       call multifab_build(adv_mass_fluxdiv_new(n),mla%la(n),nspecies,0)
       call multifab_build(      gmres_rhs_p(n),mla%la(n),1       ,0)
       call multifab_build(              dpi(n),mla%la(n),1       ,1)
       call multifab_build(             divu(n),mla%la(n),1       ,0)
       call multifab_build(             conc(n),mla%la(n),nspecies,rho_old(n)%ng)
       call multifab_build(            S_inc(n),mla%la(n),1       ,0)
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
          call multifab_build_edge(  m_grav_force_old(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(  m_grav_force_new(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge( Lorentz_force_old(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge( Lorentz_force_new(n,i),mla%la(n),1       ,0,i)
       end do
       call multifab_build(solver_alpha(n),mla%la(n),1,0)
       call multifab_build(solver_rhs(n),mla%la(n),1,0)
       do i=1,dm
          call multifab_build_edge(A_Phi(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge(solver_beta(n,i),mla%la(n),1,0,i)
          call multifab_build_edge(permittivity_fc(n,i),mla%la(n),1,0,i)
          call multifab_build_edge(zdotA(n,i),mla%la(n),1,0,i)
       end do
    end do

    do n = 1,nlevs
       call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
    end do
    do n=1,nlevs
       ! alpha=0
       ! this is used for the electric potential Poisson solve
       ! and for the evaluation of div A_Phi grad Phi
       call multifab_setval(solver_alpha(n),0.d0)
       ! this is the iterative low Mach constraint correction
       call multifab_setval(S_inc(n),0.d0)
    end do

    !!!!!!!!!!!!!!!!!!!!!!
    ! copy old densities into new, including ghost cells
    do n=1,nlevs
       call multifab_copy_c(   rho_new(n),1,   rho_old(n),1,nspecies,   rho_new(n)%ng)
       call multifab_copy_c(rhotot_new(n),1,rhotot_old(n),1,       1,rhotot_new(n)%ng)
    end do
    
    ! average rho_old and rhotot_old to faces
    call average_cc_to_face(nlevs,   rho_old,   rho_fc,1,   c_bc_comp,nspecies,the_bc_tower%bc_tower_array)
    call average_cc_to_face(nlevs,rhotot_old,rhotot_fc,1,scal_bc_comp,       1,the_bc_tower%bc_tower_array)
    
    !!!!!!!!!!!!!!!!!!!!!!
    ! t^n terms for the density update that do not change with iteration

    ! make copies of the "old" diffusive, stochastic, and potential mass fluxes
    do n=1,nlevs
       call multifab_copy_c( diff_mass_fluxdiv_old(n),1, diff_mass_fluxdiv(n),1,nspecies,0)
       call multifab_copy_c(stoch_mass_fluxdiv_old(n),1,stoch_mass_fluxdiv(n),1,nspecies,0)
       call multifab_copy_c( Epot_mass_fluxdiv_old(n),1, Epot_mass_fluxdiv(n),1,nspecies,0)
    end do

    ! compute "old" advective mass fluxes and copy into new
    do n=1,nlevs
       call multifab_setval(adv_mass_fluxdiv_old(n),0.d0)
    end do
    call mk_advective_s_fluxdiv(mla,umac,rho_fc,adv_mass_fluxdiv_old,dx,1,nspecies)
    do n=1,nlevs
       call multifab_copy_c(adv_mass_fluxdiv_new(n),1,adv_mass_fluxdiv_old(n),1,nspecies,0)
    end do
    
    !!!!!!!!!!!!!!!!!!!!!!
    ! t^n terms for gmres_rhs_v that do not change with iteration

    ! compute A_0^n v^n
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(m_d_fluxdiv_old(n,i),0.d0,all=.true.)
       end do
    end do
    call diffusive_m_fluxdiv(mla,m_d_fluxdiv_old,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    ! compute "old" stochastic momentum fluxes
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(m_s_fluxdiv_old(n,i),0.d0,all=.true.)
       end do
    end do
    if (variance_coef_mom .ne. 0.d0) then
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

    ! compute "old" Lorentz force
    call compute_Lorentz_force(mla,Lorentz_force_old,grad_Epot_old,permittivity, &
                               charge_old,dx,the_bc_tower)

    ! compute mold = rho^n v^n
    call convert_m_to_umac(mla,rhotot_fc,mold,umac,.false.)

    ! compute m_a_fluxdiv_old = div(-rho^n v^n v^n) and copy into new
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(m_a_fluxdiv_old(n,i),0.d0,all=.true.)
       end do
    end do
    call mk_advective_m_fluxdiv(mla,umac,mold,m_a_fluxdiv_old,dx,the_bc_tower%bc_tower_array)
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(m_a_fluxdiv_new(n,i),1,m_a_fluxdiv_old(n,i),1,1,0)
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!
    ! iterative loop over l
    do l=1,num_pot_iters

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Step 1 - Density Update
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if (l .gt. 1) then

          ! compute adv_mass_fluxdiv_new
          do n=1,nlevs
             call multifab_setval(adv_mass_fluxdiv_new(n),0.d0)
          end do
          call mk_advective_s_fluxdiv(mla,umac,rho_fc,adv_mass_fluxdiv_new,dx,1,nspecies)

       end if

       ! compute A_Phi^{n+1,l} to solve for Epot_mass_fluxdiv_new via Poisson solve
       call implicit_potential_coef(mla,rho_new,Temp,A_Phi,the_bc_tower)

       ! build RHS for Poisson solve (store in rho_new)
       ! first set RHS = (1-theta)*(A^n + D^n + St^n) + theta*(A^{n+1,l} + D^{n+1,l} + St^{n+1,l})
       do n=1,nlevs
          call multifab_setval(rho_new(n),0.d0,all=.true.)
          call multifab_saxpy_3_cc(rho_new(n),1,1.d0-theta_pot,diff_mass_fluxdiv_old(n),1,nspecies)
          call multifab_saxpy_3_cc(rho_new(n),1,     theta_pot,diff_mass_fluxdiv    (n),1,nspecies)
          if (variance_coef_mass .ne. 0.d0) then
             call multifab_saxpy_3_cc(rho_new(n),1,1.d0-theta_pot,stoch_mass_fluxdiv_old(n),1,nspecies)
             call multifab_saxpy_3_cc(rho_new(n),1,     theta_pot,stoch_mass_fluxdiv    (n),1,nspecies)
          end if
          call multifab_saxpy_3_cc(rho_new(n),1,1.d0-theta_pot,adv_mass_fluxdiv_old(n),1,nspecies)
          call multifab_saxpy_3_cc(rho_new(n),1,     theta_pot,adv_mass_fluxdiv_new(n),1,nspecies)
       end do

       ! add (1-theta) Epot_mass_fluxdiv_old to RHS
       do n=1,nlevs
          call multifab_saxpy_3_cc(rho_new(n),1,1.d0-theta_pot,Epot_mass_fluxdiv_old(n),1,nspecies)
       end do

       ! multiply by dt and add rho_old
       do n=1,nlevs
          call multifab_mult_mult_s_c(rho_new(n),1,dt,nspecies,0)
          call multifab_plus_plus_c(rho_new(n),1,rho_old(n),1,nspecies,0)
       end do

       ! right-hand-side for Poisson solve is z^T RHS
       call dot_with_z(mla,rho_new,solver_rhs)

       ! compute z^T A_Phi^{n+1,l}, store in solver_beta
       call dot_with_z_face(mla,A_Phi,solver_beta)

       ! permittivity on faces
       call average_cc_to_face(nlevs,permittivity,permittivity_fc,1,scal_bc_comp,1, &
                               the_bc_tower%bc_tower_array)

       ! compute solver_beta = epsilon + dt theta z^T A_Phi^{n+1,l}
       do n=1,nlevs
          do i=1,dm
             call multifab_mult_mult_s_c(solver_beta(n,i),1,dt*theta_pot,1,0)
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

       ! solve -div (epsilon + dt theta z^T A_Phi^{n+1,l}) grad Phi^{n+1,l+1} = z^T RHS
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
       
       ! compute the gradient of the electric potential for use in momentum force
       call compute_grad(mla,Epot,grad_Epot_new,dx,1,Epot_bc_comp,1,1,the_bc_tower%bc_tower_array)

       do comp=1,nspecies

          ! copy component of A_Phi^{n+1,l} into beta and multiply by grad_Epot_new
          do n=1,nlevs
             do i=1,dm
                call multifab_copy_c(solver_beta(n,i),1,A_Phi(n,i),comp,1,0)
                call multifab_mult_mult_c(solver_beta(n,i),1,grad_Epot_new(n,i),1,1,0)
             end do
             ! zero mass flux on walls
             call zero_edgeval_walls(solver_beta(n,:),1,1,the_bc_tower%bc_tower_array(n))
          end do

          ! compute Epot_mass_fluxdiv = div A_Phi^{n+1,l} grad Epot
          call compute_div(mla,solver_beta,Epot_mass_fluxdiv,dx,1,comp,1)
       end do

       ! add dt*theta*Epot_mass_fluxdiv to RHS to get rho^{n+1,l+1}
       do n=1,nlevs
          call multifab_saxpy_3_cc(rho_new(n),1,dt*theta_pot,Epot_mass_fluxdiv(n),1,nspecies)
       end do

       ! compute rhotot from rho in VALID REGION
       call compute_rhotot(mla,rho_new,rhotot_new)

       ! rho to conc - NO GHOST CELLS
       call convert_rhoc_to_c(mla,rho_new,rhotot_new,conc,.true.)
       call fill_c_ghost_cells(mla,conc,dx,the_bc_tower)

       do n=1,nlevs
          call fill_rho_ghost_cells(conc(n),rhotot_new(n),the_bc_tower%bc_tower_array(n))
       end do

       ! conc to rho - INCLUDING GHOST CELLS
       call convert_rhoc_to_c(mla,rho_new,rhotot_new,conc,.false.)

       ! print out EOS drift
       call eos_check(mla,rho_new)

       ! average rho_new and rhotot_new to faces
       call average_cc_to_face(nlevs,   rho_new,   rho_fc,1,   c_bc_comp,nspecies,the_bc_tower%bc_tower_array)
       call average_cc_to_face(nlevs,rhotot_new,rhotot_fc,1,scal_bc_comp,       1,the_bc_tower%bc_tower_array)

       ! compute total charge
       call dot_with_z(mla,rho_new,charge_new)

       ! compute new permittivity
       if (dielectric_type .ne. 0) then
          call compute_permittivity(mla,permittivity,rho_new,rhotot_new, &
                                    the_bc_tower)
       end if

       ! compute mtemp = rho^{n+1,l+1} v^{n+1,l}
       call convert_m_to_umac(mla,rhotot_fc,mtemp,umac,.false.)

       ! compute (eta,kappa)^{n+1,l+1}
       call compute_eta_kappa(mla,eta,eta_ed,kappa,rho_new,rhotot_new,Temp,dx, &
                              the_bc_tower%bc_tower_array)

       ! compute A_Phi^{n+1,l+1} to compute updated Epot_mass_fluxdiv_new
       call implicit_potential_coef(mla,rho_new,Temp,A_Phi,the_bc_tower)

       do comp=1,nspecies

          ! copy component of A_Phi^{n+1,l+1} into beta and multiply by grad_Epot_new
          do n=1,nlevs
             do i=1,dm
                call multifab_copy_c(solver_beta(n,i),1,A_Phi(n,i),comp,1,0)
                call multifab_mult_mult_c(solver_beta(n,i),1,grad_Epot_new(n,i),1,1,0)
             end do
             ! zero mass flux on walls
             call zero_edgeval_walls(solver_beta(n,:),1,1,the_bc_tower%bc_tower_array(n))
          end do

          ! compute Epot_mass_fluxdiv = div A_Phi^{n+1,l+1} grad Epot
          call compute_div(mla,solver_beta,Epot_mass_fluxdiv,dx,1,comp,1)
       end do

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Step 3 - Velocity Update
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! compute new gravity
       do n=1,nlevs
          do i=1,dm
             call multifab_setval(m_grav_force_new(n,i),0.d0,all=.true.)
          end do
       end do
       if (any(grav(1:dm) .ne. 0.d0)) then
          call mk_grav_force(mla,m_grav_force_new,rhotot_fc,rhotot_fc,the_bc_tower)
       end if

       ! compute "new" Lorentz force
       call compute_Lorentz_force(mla,Lorentz_force_new,grad_Epot_new,permittivity, &
                                  charge_new,dx,the_bc_tower)

       if (l .gt. 1) then

          ! compute m_a_fluxdiv_new = div(-rho^{n+1,l} v^{n+1,l} v^{n+1,l})
          do n=1,nlevs
             do i=1,dm
                call multifab_setval(m_a_fluxdiv_new(n,i),0.d0,all=.true.)
             end do
          end do
          call mk_advective_m_fluxdiv(mla,umac,mtemp,m_a_fluxdiv_new,dx,the_bc_tower%bc_tower_array)

       end if

       ! build gmres_rhs_v
       ! first set gmres_rhs_v = (mold - mtemp) / dt
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(gmres_rhs_v(n,i),1,mold(n,i),1,1,0)
             call multifab_sub_sub_c(gmres_rhs_v(n,i),1,mtemp(n,i),1,1,0)
             call multifab_mult_mult_s_c(gmres_rhs_v(n,i),1,1.d0/dt,1,0)
          end do
       end do

       ! compute grad pi^{n+1,l}
       call compute_grad(mla,pi,gradpi,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

       if (barodiffusion_type .ne. 0) then
          call bl_error("advance_timestep_iterative: barodiffusion not supported yet")
       end if

       ! subtract grad pi^{n+1,l} from gmres_rhs_v
       do n=1,nlevs
          do i=1,dm
             call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradpi(n,i),1,1,0)
          end do
       end do

       ! add (1/2) (A^n + A^{n+1,l})
       do n=1,nlevs
          do i=1,dm
             call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_a_fluxdiv_old(n,i),1,1)
             call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_a_fluxdiv_new(n,i),1,1)
          end do
       end do

       ! compute A_0^{n+1,l+1} vbar^{n+1,l}
       do n=1,nlevs
          do i=1,dm
             call multifab_setval(m_d_fluxdiv_new(n,i),0.d0,all=.true.)
          end do
       end do
       call diffusive_m_fluxdiv(mla,m_d_fluxdiv_new,umac,eta,eta_ed,kappa,dx, &
                                the_bc_tower%bc_tower_array)

       ! add (1/2) (A_0^n v^n + A_0^{n+1,l+1} vbar^{n+1,l})
       do n=1,nlevs
          do i=1,dm
             call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_d_fluxdiv_old(n,i),1,1)
             call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_d_fluxdiv_new(n,i),1,1)
          end do
       end do

       ! compute "new" stochastic momentum fluxes
       do n=1,nlevs
          do i=1,dm
             call multifab_setval(m_s_fluxdiv_new(n,i),0.d0,all=.true.)
          end do
       end do
       if (variance_coef_mom .ne. 0.d0) then
          call stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,m_s_fluxdiv_new, &
                                    eta,eta_ed,Temp,Temp_ed,dx,dt,weights)
       end if

       ! add (1/2) (old + new) stochastic fluxes
       do n=1,nlevs
          do i=1,dm
             call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_s_fluxdiv_old(n,i),1,1)
             call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_s_fluxdiv_new(n,i),1,1)
          end do
       end do

       ! add (1/2) (old + new) gravitational force
       do n=1,nlevs
          do i=1,dm
             call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_grav_force_old(n,i),1,1)
             call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,m_grav_force_new(n,i),1,1)
          end do
       end do

       ! subtract momentum charge force
       do n=1,nlevs
          do i=1,dm
             call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,-(1.d0-theta_pot),Lorentz_force_old(n,i),1,1)
             call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,      -theta_pot ,Lorentz_force_new(n,i),1,1)
          end do
       end do

       ! new random fluxes? - is this right?
       if (l .eq. num_pot_iters .and. variance_coef_mass .ne. 0.d0) then
          call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)
       end if

       ! compute diff_mass_fluxdiv_new and stoch_mass_fluxdiv_new for gmres_rhs_p
       call compute_mass_fluxdiv_charged(mla,rho_new,gradp_baro, &
                                         diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                         Temp,flux_total,flux_diff, &
                                         dt,time,dx,weights,the_bc_tower)

       if(present(gradPhiApprox) .and. (l==num_pot_iters)) then 
          ! compute grad(phi) approximation, z^T Fmass / z^T A_Phi
          call dot_with_z_face(mla,flux_diff,gradPhiApprox)
          call dot_with_z_face(mla,A_Phi,zdotA)
          do n=1,nlevs
             do i=1,dm
                call multifab_div_div_c(gradPhiApprox(n,i),1,zdotA(n,i),1,1,0)
             end do
          end do
       end if

       ! now fluxes contain "-F = rho*W*chi*Gamma*grad(x) + ..."
       do n=1,nlevs
          call multifab_mult_mult_s_c(diff_mass_fluxdiv(n),1,-1.d0,nspecies,0)
          if (variance_coef_mass .ne. 0) then
             call multifab_mult_mult_s_c(stoch_mass_fluxdiv    (n),1,-1.d0,nspecies,0)
          end if
          do i=1,dm
             call multifab_mult_mult_s_c(flux_total(n,i),1,-1.d0,nspecies,0)
          end do
       end do

       ! set the Dirichlet velocity value on reservoir faces
       call reservoir_bc_fill(mla,flux_total,vel_bc_n,the_bc_tower%bc_tower_array)

       ! put "-S = div(F_i/rho_i)" into gmres_rhs_p (we will later add divu)
       do n=1,nlevs
          call multifab_setval(gmres_rhs_p(n),0.d0,all=.true.)
          do i=1,nspecies
             call multifab_saxpy_3_cc(gmres_rhs_p(n),1,-1.d0/rhobar(i), diff_mass_fluxdiv(n),i,1)
             call multifab_saxpy_3_cc(gmres_rhs_p(n),1,-1.d0/rhobar(i), Epot_mass_fluxdiv(n),i,1)
             if (variance_coef_mass .ne. 0.d0) then
                call multifab_saxpy_3_cc(gmres_rhs_p(n),1,-1.d0/rhobar(i),stoch_mass_fluxdiv(n),i,1)
             end if
          end do
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

       ! compute div vbar^{n+1,l}
       call compute_div(mla,umac,divu,dx,1,1,1)

       ! add div vbar^{n+1,l} to gmres_rhs_p
       ! now gmres_rhs_p = div vbar^{n+1,l} - S^{n+1,l+1}
       ! the sign convention is correct since we solve -div(delta v) = gmres_rhs_p
       do n=1,nlevs
          call multifab_plus_plus_c(gmres_rhs_p(n),1,divu(n),1,1,0)
       end do

       ! Compute Velocity Constraint Correction
       if (dpdt_factor .ne. 0.d0) then
          call modify_S(mla,rho_new,S_inc,dt)
          do n=1,nlevs
             call multifab_plus_plus_c(gmres_rhs_p(n),1,S_inc(n),1,1,0)
          end do
       end if

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

       if (l .eq. 1) then
          gmres_abs_tol_in = gmres_abs_tol ! Save this 
       end if

       ! call gmres to compute delta v and delta pi
       call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dpi,rhotot_fc, &
                  eta,eta_ed,kappa,theta_alpha,norm_pre_rhs)

       ! for the corrector gmres solve we want the stopping criteria based on the
       ! norm of the preconditioned rhs from the predictor gmres solve.  otherwise
       ! for cases where du in the corrector should be small the gmres stalls
       if (l .eq. 1) then
          gmres_abs_tol = max(gmres_abs_tol_in, norm_pre_rhs*gmres_rel_tol)
       end if

       ! restore eta and kappa
       do n=1,nlevs
          call multifab_mult_mult_s_c(eta(n)  ,1,2.d0,1,eta(n)%ng)
          call multifab_mult_mult_s_c(kappa(n),1,2.d0,1,kappa(n)%ng)
          do i=1,size(eta_ed,dim=2)
             call multifab_mult_mult_s_c(eta_ed(n,i),1,2.d0,1,eta_ed(n,i)%ng)
          end do
       end do
       
       ! compute v^{n+1,l+1} = v^{n+1,l} + dumac
       ! compute pi^{n+1,l+1}= pi^{n+1,l} + dpi
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

    end do ! end loop l=1,num_pot_iters

    gmres_abs_tol = gmres_abs_tol_in ! Restore the desired tolerance   

    ! fill the stochastic momentum multifabs with a new set of random numbers
    if (variance_coef_mom .ne. 0.d0) then
       call fill_m_stochastic(mla)
    end if

    call destroy_bc_multifabs(mla)

    do n=1,nlevs
       call multifab_destroy(diff_mass_fluxdiv_old(n))
       call multifab_destroy(stoch_mass_fluxdiv_old(n))
       call multifab_destroy(Epot_mass_fluxdiv_old(n))
       call multifab_destroy(adv_mass_fluxdiv_old(n))
       call multifab_destroy(adv_mass_fluxdiv_new(n))
       call multifab_destroy(gmres_rhs_p(n))
       call multifab_destroy(dpi(n))
       call multifab_destroy(divu(n))
       call multifab_destroy(conc(n))
       call multifab_destroy(S_inc(n))
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
          call multifab_destroy(m_grav_force_old(n,i))
          call multifab_destroy(m_grav_force_new(n,i))
          call multifab_destroy(Lorentz_force_old(n,i))
          call multifab_destroy(Lorentz_force_new(n,i))
       end do
       call multifab_destroy(solver_alpha(n))
       call multifab_destroy(solver_rhs(n))
       do i=1,dm
          call multifab_destroy(A_Phi(n,i))
          call multifab_destroy(solver_beta(n,i))
          call multifab_destroy(permittivity_fc(n,i))
          call multifab_destroy(zdotA(n,i))
       end do
    end do
    do n = 1,nlevs
       call bndry_reg_destroy(fine_flx(n))
    end do

  end subroutine advance_timestep_iterative

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

end module advance_timestep_iterative_module
