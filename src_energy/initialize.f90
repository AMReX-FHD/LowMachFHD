module initialize_module

  use ml_layout_module
  use bndry_reg_module
  use ml_solve_module
  use multifab_physbc_module
  use bc_module
  use define_bc_module
  use energy_eos_module
  use energy_eos_wrapper_module
  use convert_rhoc_to_c_module
  use convert_rhoh_to_h_module
  use mass_fluxdiv_energy_module
  use rhoh_fluxdiv_energy_module
  use macproject_module
  use convert_stag_module
  use div_and_grad_module
  use mk_advective_s_fluxdiv_module
  use mass_flux_utilities_module
  use probin_multispecies_module, only: nspecies
  use probin_common_module, only: n_cells

  use fabio_module

  implicit none

  private

  public :: initialize

contains

  ! this routine performs "Step 0" of the algorithm (see exec/energy/doc/)
  ! -Compute P_0^{*,n+1}
  ! -Compute v^n with a projection
  ! -Compute rho_i^{*,n+1} explicitly
  ! -Compute (rho h)^{*,n+1} implicitly
  ! -If necessary, compute volume discrepancy correction and return to beginning of step
  subroutine initialize(mla,umac_old,rho_old,rho_new,rhotot_old,rhotot_new, &
                        rhoh_old,rhoh_new,p0_old,p0_new, &
                        gradp_baro,Temp_old,Temp_new, &
                        dx,dt,time,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: umac_old(:,:)
    type(multifab) , intent(inout) :: rho_old(:)
    type(multifab) , intent(inout) :: rho_new(:)
    type(multifab) , intent(inout) :: rhotot_old(:)
    type(multifab) , intent(inout) :: rhotot_new(:)
    type(multifab) , intent(inout) :: rhoh_old(:)
    type(multifab) , intent(inout) :: rhoh_new(:)
    real(kind=dp_t), intent(in   ) :: p0_old
    real(kind=dp_t), intent(inout) :: p0_new
    type(multifab) , intent(inout) :: gradp_baro(:,:)
    type(multifab) , intent(inout) :: Temp_old(:)
    type(multifab) , intent(inout) :: Temp_new(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt,time
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables

    ! temporary copy of initial umac
    type(multifab) :: umac_tmp(mla%nlevel,mla%dim)

    ! this will hold F_k and div(F_k)
    type(multifab) :: mass_flux_old(mla%nlevel,mla%dim)
    type(multifab) :: mass_flux_new(mla%nlevel,mla%dim)
    type(multifab) :: mass_fluxdiv_old(mla%nlevel)
    type(multifab) :: mass_fluxdiv_new(mla%nlevel)

    ! this will hold div(Q) + sum(div(hk*Fk)) + rho*Hext
    type(multifab) :: rhoh_fluxdiv_old(mla%nlevel)
    type(multifab) :: rhoh_fluxdiv_new(mla%nlevel)

    ! this holds h
    type(multifab) :: enth(mla%nlevel)

    ! this will hold rho and rhoh on faces
    type(multifab) :: rho_fc(mla%nlevel,mla%dim)
    type(multifab) :: rhoh_fc(mla%nlevel,mla%dim)

    ! This will hold (rhoh)^n/dt - div(rhoh*v)^n + (Sbar^n+Sbarcorr^n)/alphabar^n 
    !                + (1/2)(div(Q^n) + sum(div(h_k^n F_k^n)) + (rho Hext)^n)
    ! for the RHS of the temperature diffusion solve.
    ! Each of these terms stays fixed over all l iterations.
    type(multifab) :: deltaT_rhs1(mla%nlevel)

    ! This will hold -(rho^{*,n+1}h^{*,n+1,l})/dt
    !                + (1/2)(div(Q^{*,n+1,l}) + sum(div(h_k^{*,n+1,l}F_k^{*,n+1,l}))
    !                + (rho Hext)^(*,n+1))
    ! for the RHS of the temperature diffusion solve.
    ! Each of these terms may change for each l iteration.
    type(multifab) :: deltaT_rhs2(mla%nlevel)

    ! temporary storage for concentrations and mole fractions
    type(multifab) :: conc    (mla%nlevel)
    type(multifab) :: molefrac(mla%nlevel)

    ! the implicit energy solve computes deltaT
    type(multifab) :: deltaT(mla%nlevel)

    ! coefficients for deltaT (energy) solve
    type(multifab) :: cc_solver_alpha(mla%nlevel)
    type(multifab) :: cc_solver_beta (mla%nlevel,mla%dim)

    ! shear viscosity
    type(multifab) :: eta_old(mla%nlevel)
    type(multifab) :: eta_new(mla%nlevel)

    ! bulk viscosity
    type(multifab) :: kappa_old(mla%nlevel)
    type(multifab) :: kappa_new(mla%nlevel)

    ! thermal diffusivity
    type(multifab) :: lambda_old(mla%nlevel)
    type(multifab) :: lambda_new(mla%nlevel)

    ! diffusion matrix
    type(multifab) :: chi_old(mla%nlevel)
    type(multifab) :: chi_new(mla%nlevel)

    ! thermodiffusion coefficients
    type(multifab) :: zeta_old(mla%nlevel)
    type(multifab) :: zeta_new(mla%nlevel)

    ! div(u) + alpha dP_0/dt = S_old, where
    ! S_old = Sbar_old + deltaS_old
    type(multifab)  :: deltaS_old(mla%nlevel)
    real(kind=dp_t) :: Sbar_old

    ! volume discrepancy correction
    ! Scorr = Scorrbar + deltaScorr
    type(multifab)  :: Scorr_old(mla%nlevel)
    type(multifab)  :: deltaScorr_old(mla%nlevel)
    real(kind=dp_t) :: Scorrbar_old

    ! coefficient multiplying dP_0/dt in constraint at old-time
    ! alpha = alphabar + deltaalpha
    type(multifab)  :: deltaalpha_old(mla%nlevel)
    real(kind=dp_t) :: alphabar_old

    ! coefficient multiplying dP_0/dt in volume discrepancy correction at new-time
    type(multifab)  :: alpha_new(mla%nlevel)

    ! the RHS for the projection
    type(multifab) :: Sproj(mla%nlevel)
    ! solution of the pressure-projection solve
    type(multifab) :: phi(mla%nlevel)
    ! coefficient for projection
    type(multifab) :: rhotot_fc(mla%nlevel,mla%dim)
    type(multifab) :: rhototinv_fc(mla%nlevel,mla%dim)

    ! this holds the thermodynamic pressure
    type(multifab) :: Peos(mla%nlevel)

    ! for energy implicit solve
    ! doesn't actually do anything for single-level solves
    type(bndry_reg) :: fine_flx(2:mla%nlevel)

    integer :: n,nlevs,i,dm,n_cell,k,l

    nlevs = mla%nlevel
    dm = mla%dim

    if (dm .eq. 2) then
       n_cell = n_cells(1)*n_cells(2)
    else
       n_cell = n_cells(1)*n_cells(2)*n_cells(3)
    end if

    do n=1,nlevs
       do i=1,dm
          call multifab_build_edge(umac_tmp(n,i),mla%la(n),1,1,i)
          call multifab_build_edge(mass_flux_old(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge(mass_flux_new(n,i),mla%la(n),nspecies,0,i)
       end do
       call multifab_build(mass_fluxdiv_old(n),mla%la(n),nspecies,0)
       call multifab_build(mass_fluxdiv_new(n),mla%la(n),nspecies,0)

       call multifab_build(enth(n),mla%la(n),1,rho_old(n)%ng)

       do i=1,dm
          call multifab_build_edge(rho_fc(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge(rhoh_fc(n,i),mla%la(n),1,0,i)
       end do

       call multifab_build(rhoh_fluxdiv_old(n),mla%la(n),1,0)
       call multifab_build(rhoh_fluxdiv_new(n),mla%la(n),1,0)

       call multifab_build(deltaT_rhs1(n),mla%la(n),1,0)
       call multifab_build(deltaT_rhs2(n),mla%la(n),1,0)

       call multifab_build(    conc(n),mla%la(n),nspecies,rho_old(n)%ng)
       call multifab_build(molefrac(n),mla%la(n),nspecies,rho_old(n)%ng)

       call multifab_build(deltaT(n),mla%la(n),1,1)

       call multifab_build(cc_solver_alpha(n),mla%la(n),1,0)
       do i=1,dm
          call multifab_build_edge(cc_solver_beta(n,i),mla%la(n),1,0,i)
       end do

       call multifab_build(eta_old(n),mla%la(n),1,1)
       call multifab_build(eta_new(n),mla%la(n),1,1)

       call multifab_build(kappa_old(n),mla%la(n),1,1)
       call multifab_build(kappa_new(n),mla%la(n),1,1)

       call multifab_build(lambda_old(n),mla%la(n),1,2)
       call multifab_build(lambda_new(n),mla%la(n),1,2)

       call multifab_build(chi_old(n),mla%la(n),nspecies**2,1)
       call multifab_build(chi_new(n),mla%la(n),nspecies**2,1)

       call multifab_build(zeta_old(n),mla%la(n),nspecies,1)
       call multifab_build(zeta_new(n),mla%la(n),nspecies,1)

       call multifab_build(deltaS_old(n),mla%la(n),1,0)

       call multifab_build(     Scorr_old(n),mla%la(n),1,0)
       call multifab_build(deltaScorr_old(n),mla%la(n),1,0)

       call multifab_build(     alpha_new(n),mla%la(n),1,0)
       call multifab_build(deltaalpha_old(n),mla%la(n),1,0)

       call multifab_build(Sproj(n),mla%la(n),1,0)
       call multifab_build(phi(n),mla%la(n),1,1)
       do i=1,dm
          call multifab_build_edge(   rhotot_fc(n,i),mla%la(n),1,0,i)
          call multifab_build_edge(rhototinv_fc(n,i),mla%la(n),1,0,i)
       end do

       call multifab_build(Peos(n),mla%la(n),1,0)

    end do

    ! make a copy of the initial velocity
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(umac_tmp(n,i),1,umac_old(n,i),1,1,1)
       end do
    end do

    ! compute mass fractions in valid region and then fill ghost cells
    call convert_rhoc_to_c(mla,rho_old,rhotot_old,conc,.true.)
    call fill_c_ghost_cells(mla,conc,dx,the_bc_tower)

    ! compute mole fractions in VALID + GHOST regions
    call convert_conc_to_molefrac(mla,conc,molefrac,.true.)

    ! compute initial transport properties
    call ideal_mixture_transport_wrapper(mla,rhotot_old,Temp_old,p0_old,conc,molefrac, &
                                         eta_old,lambda_old,kappa_old,chi_old,zeta_old)

    ! compute mass_fluxdiv_old = div(F^n)
    call mass_fluxdiv_energy(mla,rho_old,rhotot_old,molefrac,chi_old,zeta_old, &
                             gradp_baro,Temp_old,mass_fluxdiv_old,mass_flux_old,dx,the_bc_tower)

    ! compute rhoh_fluxdiv_old = div(Q)^n + sum(div(hk*Fk))^n + rho*Hext^n
    call rhoh_fluxdiv_energy(mla,lambda_old,Temp_old,mass_flux_old,rhotot_old,rhoh_fluxdiv_old, &
                             dx,time,the_bc_tower)

    ! compute S_old and alpha_old (store them in deltaS_old and deltaalpha_old)
    call compute_S_alpha(mla,deltaS_old,deltaalpha_old,mass_fluxdiv_old,rhoh_fluxdiv_old,conc, &
                         Temp_old,rhotot_old,p0_old)

    ! split S_old and alpha_old into average and perturbational pieces
    ! S_old = Sbar_old + deltaS_old
    ! alpha_old = alphabar_old + deltaalpha_old
    do n=1,nlevs
       Sbar_old     = multifab_sum_c(deltaS_old(n)    ,1,1) / dble(n_cell)
       alphabar_old = multifab_sum_c(deltaalpha_old(n),1,1) / dble(n_cell)
       call multifab_sub_sub_s_c(deltaS_old(n)    ,1,Sbar_old    ,1,0)
       call multifab_sub_sub_s_c(deltaalpha_old(n),1,alphabar_old,1,0)
    end do

    ! zero out volume discrepancy correction and its decomposition
    Scorrbar_old = 0.d0
    do n=1,nlevs
       call multifab_setval(Scorr_old(n)     ,0.d0,all=.true.)
       call multifab_setval(deltaScorr_old(n),0.d0,all=.true.)
    end do

    ! average rhotot to faces
    call average_cc_to_face(nlevs,rhotot_old,rhotot_fc,1,scal_bc_comp,1, &
                            the_bc_tower%bc_tower_array)

    ! compute (1/rhotot) on faces)
    do n=1,nlevs
       do i=1,dm
          call setval(rhototinv_fc(n,i),1.d0,all=.true.)
          call multifab_div_div_c(rhototinv_fc(n,i),1,rhotot_fc(n,i),1,1,0)
       end do
    end do

    ! begin loop here over Steps 0a-0e
    do k=1,dpdt_iters

       ! hack - need to set umac_old back to its initial value
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(umac_old(n,i),1,umac_tmp(n,i),1,1,1)
          end do
       end do

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Step 0a: Compute a pressure update
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! update pressure
       p0_new = p0_old + dt*(Sbar_old + Scorrbar_old)/alphabar_old

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Step 0b: Compute the velocity field using a projection
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! compute Sproj = deltaS_old + deltaScorr_old 
       !                 - deltaalpha_old * (Sbar_old + Scorrbar_old)/alphabar_old
       do n=1,nlevs
          call multifab_copy_c(Sproj(n),1,deltaalpha_old(n),1,1,0)
          call multifab_mult_mult_s_c(Sproj(n),1,-(Sbar_old+Scorrbar_old)/alphabar_old,1,0)
          call multifab_plus_plus_c(Sproj(n),1,deltaS_old(n),1,1,0)
          call multifab_plus_plus_c(Sproj(n),1,deltaScorr_old(n),1,1,0)
       end do

       ! build rhs for projection, div(v^init) - Sproj
       ! first multiply Sproj by -1
       do n=1,nlevs
          call multifab_mult_mult_s_c(Sproj(n),1,-1.d0,1,0)
       end do

       ! add div(v^init) to Sproj
       call compute_div(mla,umac_old,Sproj,dx,1,1,1,increment_in=.true.)

       ! solve div (1/rhotot) grad phi = div(v^init) - S^0
       ! solve to completion, i.e., use the 'full' solver
       call macproject(mla,phi,umac_old,rhototinv_fc,Sproj,dx,the_bc_tower,.true.)

       ! v^0 = v^init - (1/rho^0) grad phi
       call subtract_weighted_gradp(mla,umac_old,rhototinv_fc,phi,dx,the_bc_tower)

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Step 0c: Advance the densities using forward-Euler advective fluxes and 
       !          explicit mass diffusion
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       call average_cc_to_face(nlevs,rho_old,rho_fc,1,c_bc_comp,nspecies, &
                               the_bc_tower%bc_tower_array)

       ! compute -div(rho*v) and store it in rho_new
       call mk_advective_s_fluxdiv(mla,umac_old,rho_fc,rho_new,dx,1,nspecies)

       ! rho_new = rho_old + dt(-div(rho*v) + div(F))
       do n=1,nlevs
          call multifab_plus_plus_c(rho_new(n),1,mass_fluxdiv_old(n),1,nspecies,0)
          call multifab_mult_mult_s_c(rho_new(n),1,dt,nspecies,0)
          call multifab_plus_plus_c(rho_new(n),1,rho_old(n),1,nspecies,0)
       end do

       ! compute rhotot_new = sum(rho_new)
       call compute_rhotot(mla,rho_new,rhotot_new)

       ! fill ghost cells
       do n=1,nlevs
          call multifab_fill_boundary(rhotot_new(n))
          call multifab_physbc(rhotot_new(n),1,scal_bc_comp,1, &
                               the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
       end do

       ! compute mass fractions in valid region and then fill ghost cells
       call convert_rhoc_to_c(mla,rho_new,rhotot_new,conc,.true.)
       call fill_c_ghost_cells(mla,conc,dx,the_bc_tower)

       ! compute mole fractions in VALID + GHOST regions
       call convert_conc_to_molefrac(mla,conc,molefrac,.true.)

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Step 0d: Advance the enthalpy by iteratively looping over an energy solve
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! The portion of the RHS that stays fixed over all l iterations is:
       !   (rhoh)^n/dt - div(rhoh*v)^n + (Sbar^n+Scorrbar^n)/alphabar^n  
       !               + (1/2)(div(Q^n) + sum(div(h_k^n F_k^n)) + (rho Hext)^n)
       ! store this in deltaT_rhs1

       ! set deltaT_rhs1 = (rhoh)^n/dt + (Sbar^n+Scorrbar^n)/alphabar^n  
       do n=1,nlevs
          call multifab_copy_c(deltaT_rhs1(n),1,rhoh_old(n),1,1,0)
          call multifab_mult_mult_s_c(deltaT_rhs1(n),1,1.d0/dt,1,0)
          call multifab_plus_plus_s_c(deltaT_rhs1(n),1,(Sbar_old+Scorrbar_old)/alphabar_old,1,0)
       end do

       ! compute h
       call convert_rhoh_to_h(mla,rhoh_old,rhotot_old,enth,.true.)
       call fill_h_ghost_cells(mla,enth,dx,the_bc_tower)

       ! average h to faces
       call average_cc_to_face(nlevs,enth,rhoh_fc,1,h_bc_comp,1, &
                               the_bc_tower%bc_tower_array)

       ! multiply h on faces by rho on faces
       do n=1,nlevs
          do i=1,dm
             call multifab_mult_mult_c(rhoh_fc(n,i),1,rhotot_fc(n,i),1,1,0)
          end do
       end do

       ! add -div(rhoh*v)^n to deltaT_rhs1
       call mk_advective_s_fluxdiv(mla,umac_old,rhoh_fc,deltaT_rhs1,dx,1,1)

       ! add (1/2)(div(Q^n) + sum(div(h_k^n F_k^n)) + (rho Hext)^n) to deltaT_rhs1
       do n=1,nlevs
          call multifab_saxpy_3(deltaT_rhs1(n),0.5d0,rhoh_fluxdiv_old(n))
       end do

       ! set (rhoh)^{*,n+1,l} = rho^{*,n+1} * h^n
       do n=1,nlevs
          call multifab_copy_c(rhoh_new(n),1,enth(n),1,1,rhoh_new(n)%ng)
          call multifab_mult_mult_c(rhoh_new(n),1,rhotot_new(n),1,1,rhoh_new(n)%ng)
       end do

       ! set Temp^{*,n+1,l} = T^n
       do n=1,nlevs
          call multifab_copy_c(Temp_new(n),1,Temp_old(n),1,1,Temp_new(n)%ng)
       end do

       do l=1,deltaT_iters

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Step 0d-1: Compute (lambda,cp,F)^{*,n+1,l}
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! compute time-advanced transport properties
          call ideal_mixture_transport_wrapper(mla,rhotot_new,Temp_new,p0_new,conc,molefrac, &
                                         eta_new,lambda_new,kappa_new,chi_new,zeta_new)

          ! compute mass_fluxdiv = div(F^{*,n+1,l}))
          call mass_fluxdiv_energy(mla,rho_new,rhotot_new,molefrac,chi_new,zeta_new, &
                                   gradp_baro,Temp_new,mass_fluxdiv_new,mass_flux_new,dx,the_bc_tower)

          ! compute rhoh_fluxdiv_new = div(Q)^{*,n+1,l} + sum(div(hk*Fk))^{*,n+1,l} + rho*Hext^{*,n+1,l}
          call rhoh_fluxdiv_energy(mla,lambda_new,Temp_new,mass_flux_new,rhotot_new,rhoh_fluxdiv_new, &
                                   dx,time,the_bc_tower)

          ! The portion of the RHS that changes for each l iteration is:
          !  -(rho^{*,n+1}h^{*,n+1,l})/dt 
          !               + (1/2)(div(Q^{*,n+1,l}) + sum(div(h_k^{*,n+1,l}F_k^{*,n+1,l})) + (rho Hext)^(*,n+1))
          ! store this in deltaT_rhs2

          ! set deltaT_rhs2 = -(rho^{*,n+1}h^{*,n+1,l})/dt 
          do n=1,nlevs
             call multifab_copy_c(deltaT_rhs2(n),1,rhoh_new(n),1,1,0)
             call multifab_mult_mult_s_c(deltaT_rhs2(n),1,-1.d0/dt,1,0)
          end do

          ! add (1/2)(div(Q) + sum(div(h_k F_k)) + (rho Hext))^{*,n+1,l} to deltaT_rhs2
          do n=1,nlevs
             call multifab_saxpy_3(deltaT_rhs2(n),0.5d0,rhoh_fluxdiv_new(n))
          end do

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Step 0d-2: Solve for deltaT implicitly
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! cc_solver_alpha = rho^{*,n+1} c_p^{*,n+1,l} / dt
          call compute_cp(mla,cc_solver_alpha,conc,Temp_new)
          do n=1,nlevs
             call multifab_mult_mult_c(cc_solver_alpha(n),1,rhotot_new(n),1,1,0)
             call multifab_mult_mult_s_c(cc_solver_alpha(n),1,1.d0/dt,1,0)
          end do

          ! cc_solver_beta = (1/2) lambda^{*,n+1,l}
          call average_cc_to_face(nlevs,lambda_new,cc_solver_beta,1,tran_bc_comp,1, &
                                  the_bc_tower%bc_tower_array)
          
          do n=1,nlevs
             do i=1,dm
                call multifab_mult_mult_s_c(cc_solver_beta(n,i),1,0.5d0,1,0)
             end do
          end do

          ! cc_solver_rhs = deltaT_rhs1 + deltaT_rhs2 (store this in deltaT_rhs2)
          do n=1,nlevs
             call multifab_plus_plus_c(deltaT_rhs2(n),1,deltaT_rhs1(n),1,1,0)
          end do

          ! solve for deltaT
          do n = 2,nlevs
             call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
          end do

          do n=1,nlevs
             call setval(deltaT(n),0.d0,all=.true.)
          end do

          call ml_cc_solve(mla,deltaT_rhs2,deltaT,fine_flx,cc_solver_alpha,cc_solver_beta,dx, &
                           the_bc_tower,temp_bc_comp)

          do n = 2,nlevs
             call bndry_reg_destroy(fine_flx(n))
          end do

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Step 0d-3: Update the temperature and enthalpy
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! T^{*,n+1,l+1} = T^{*,n+1,l} + deltaT
          do n=1,nlevs
             call multifab_plus_plus_c(Temp_new(n),1,deltaT(n),1,1,0)
          end do

          ! fill T ghost cells
          do n=1,nlevs
             call multifab_fill_boundary(Temp_new(n))
             call multifab_physbc(Temp_new(n),1,temp_bc_comp,1,the_bc_tower%bc_tower_array(n), &
                                  dx_in=dx(n,:))
          end do

          ! h^{*,n+1,l+1} = h(rho^{*,n+1},T^{*,n+1,l+1})
          call compute_h(mla,Temp_new,rhoh_new)
          
          do n=1,nlevs
             call multifab_mult_mult_c(rhoh_new(n),1,rho_new(n),1,1,1)
          end do

       end do ! end loop l over deltaT iterations

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Step 0e: If the thermodynamic drift is unacceptable, update the volume
       !          discrepancy correction
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! compute thermodynamic pressure
       call compute_p(mla,rhotot_new,Temp_new,conc,Peos)

       ! compute alpha^{*,n+1}
       call compute_alpha(mla,alpha_new,conc,Temp_new,p0_new)

       ! Scorr = Scorr + alpha^{*,n+1} * [(Peos^{*,n+1} - P0^{*,n+1})/dt]
       do n=1,nlevs
          call multifab_sub_sub_s_c(Peos(n),1,p0_new,1,0)

!          if (k .eq. 1) then
!             call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift1")
!          else if (k .eq. 2) then
!             call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift2")
!          else if (k .eq. 3) then
!             call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift3")
!          else if (k .eq. 4) then
!             call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift4")
!          else if (k .eq. 5) then
!             call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift5")
!          else if (k .eq. 6) then
!             call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift6")
!          else if (k .eq. 7) then
!             call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift7")
!          else if (k .eq. 8) then
!             call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift8")
!          else if (k .eq. 9) then
!             call fabio_ml_multifab_write_d(Peos,mla%mba%rr(:,1),"a_drift9")
!          end if

          call multifab_mult_mult_s_c(Peos(n),1,1.d0/dt,1,0)
          call multifab_mult_mult_c(Peos(n),1,alpha_new(n),1,1,0)
          call multifab_plus_plus_c(Scorr_old(n),1,Peos(n),1,1,0)
       end do

       ! split Scorr into average and perturbational components
       do n=1,nlevs
          Scorrbar_old = multifab_sum_c(Scorr_old(n),1,1) / dble(n_cell)
          call multifab_copy_c(deltaScorr_old(n),1,Scorr_old(n),1,1,0)
          call multifab_sub_sub_s_c(deltaScorr_old(n),1,Scorrbar_old,1,0)
       end do

    end do  ! end loop k over dpdt iterations

    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(umac_tmp(n,i))
          call multifab_destroy(mass_flux_old(n,i))
          call multifab_destroy(mass_flux_new(n,i))
       end do
       call multifab_destroy(mass_fluxdiv_old(n))
       call multifab_destroy(mass_fluxdiv_new(n))
       call multifab_destroy(enth(n))
       do i=1,dm
          call multifab_destroy(rho_fc(n,i))
          call multifab_destroy(rhoh_fc(n,i))
       end do
       call multifab_destroy(rhoh_fluxdiv_old(n))
       call multifab_destroy(rhoh_fluxdiv_new(n))
       call multifab_destroy(deltaT_rhs1(n))
       call multifab_destroy(deltaT_rhs2(n))
       call multifab_destroy(conc(n))
       call multifab_destroy(molefrac(n))
       call multifab_destroy(deltaT(n))
       call multifab_destroy(cc_solver_alpha(n))
       do i=1,dm
          call multifab_destroy(cc_solver_beta(n,i))
       end do
       call multifab_destroy(eta_old(n))
       call multifab_destroy(eta_new(n))
       call multifab_destroy(kappa_old(n))
       call multifab_destroy(kappa_new(n))
       call multifab_destroy(lambda_old(n))
       call multifab_destroy(lambda_new(n))
       call multifab_destroy(chi_old(n))
       call multifab_destroy(chi_new(n))
       call multifab_destroy(zeta_old(n))
       call multifab_destroy(zeta_new(n))
       call multifab_destroy(deltaS_old(n))
       call multifab_destroy(Scorr_old(n))
       call multifab_destroy(deltaScorr_old(n))
       call multifab_destroy(alpha_new(n))
       call multifab_destroy(deltaalpha_old(n))
       call multifab_destroy(Sproj(n))
       call multifab_destroy(phi(n))
       do i=1,dm
          call multifab_destroy(rhotot_fc(n,i))
          call multifab_destroy(rhototinv_fc(n,i))
       end do
       call multifab_destroy(Peos(n))
    end do

    stop

  end subroutine initialize

end module initialize_module
