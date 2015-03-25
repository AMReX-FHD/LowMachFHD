module initialize_module

  use ml_layout_module
  use multifab_physbc_module
  use bc_module
  use define_bc_module
  use energy_eos_module
  use energy_eos_wrapper_module
  use convert_variables_module
  use mass_fluxdiv_energy_module
  use rhoh_fluxdiv_energy_module
  use probin_multispecies_module, only: nspecies
  use probin_common_module, only: n_cells

  use fabio_module

  implicit none

  private

  public :: initialize

contains

  ! this routine performs "Step 0" of the algorithm:
  ! -Update P_0 and compute v^n with a projection
  ! -Advance rho_i and (rho h).
  ! -If necessary, compute volume discrepancy correction and return to 
  !  projection part of this step
  subroutine initialize(mla,umac,rho_old,rho_new,rhotot_old,rhotot_new, &
                        rhoh_old,rhoh_new,p0_old,p0_new, &
                        gradp_baro,pi,Temp, &
                        dx,dt,time,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(inout) :: rho_old(:)
    type(multifab) , intent(inout) :: rho_new(:)
    type(multifab) , intent(in   ) :: rhotot_old(:)
    type(multifab) , intent(inout) :: rhotot_new(:)
    type(multifab) , intent(in   ) :: rhoh_old(:)
    type(multifab) , intent(inout) :: rhoh_new(:)
    real(kind=dp_t), intent(in   ) :: p0_old
    real(kind=dp_t), intent(inout) :: p0_new
    type(multifab) , intent(inout) :: gradp_baro(:,:)
    type(multifab) , intent(inout) :: pi(:)
    type(multifab) , intent(inout) :: Temp(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt,time
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables

    ! this will hold -div(rho*v)^n + div(F^n)
    type(multifab) :: rho_update(mla%nlevel)

    ! this will hold F_k and div(F_k)
    type(multifab) :: mass_flux(mla%nlevel,mla%dim)
    type(multifab) :: mass_fluxdiv(mla%nlevel)

    ! this will hold div(Q) + sum(div(hk*Fk)) + rho*Hext
    type(multifab) :: rhoh_fluxdiv(mla%nlevel)

    ! This will hold (rhoh)^n/dt - div(rhoh*v)^n + Sbar^n/alphabar^n 
    !                + (1/2)(div(Q^n) + sum(div(h_k^n F_k^n)) + (rho Hext)^n)
    ! for the RHS of the temperature diffusion solve.
    ! Each of these terms stays fixed over all l iterations.
    type(multifab) :: deltaT_rhs1(mla%nlevel)

    ! This will hold -(rho^{*,n+1}h^{*,n+1,l})/dt + Sbarcorr^n/alphabar^n
    !                + (1/2)(div(Q^{*,n+1,l}) + sum(div(h_k^{*,n+1,l}F_k^{*,n+1,l})
    !                + (rho Hext)^(*,n+1)
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

    ! this is the div(u)=S
    ! S = Sbar + deltaS
    type(multifab)  :: deltaS(mla%nlevel)
    real(kind=dp_t) :: Sbar

    ! volume discrepancy correction
    ! Scorr = Scorrbar + deltaScorr
    type(multifab)  :: Scorr     (mla%nlevel)
    type(multifab)  :: deltaScorr(mla%nlevel)
    real(kind=dp_t) :: Scorrbar

    ! coefficient multiplying dP_0/dt in constraint
    ! alpha = alphabar + deltaalpha
    type(multifab)  :: deltaalpha(mla%nlevel)
    real(kind=dp_t) :: alphabar

    ! the RHS for the projection
    type(multifab) :: Sproj(mla%nlevel)
    ! solution of the pressure-projection solve
    type(multifab) :: phi(mla%nlevel)
    ! coefficient for projection
    type(multifab) :: rhoinv_fc(mla%nlevel,mla%dim)

    integer :: n,nlevs,i,dm,n_cell

    nlevs = mla%nlevel
    dm = mla%dim

    if (dm .eq. 2) then
       n_cell = n_cells(1)*n_cells(2)
    else
       n_cell = n_cells(1)*n_cells(2)*n_cells(3)
    end if

    do n=1,nlevs
       call multifab_build(rho_update(n),mla%la(n),nspecies,0)

       do i=1,dm
          call multifab_build_edge(mass_flux(n,i),mla%la(n),nspecies,0,i)
       end do
       call multifab_build(mass_fluxdiv(n),mla%la(n),nspecies,0)
       call multifab_build(rhoh_fluxdiv(n),mla%la(n),1,0)

       call multifab_build(deltaT_rhs1(n),mla%la(n),1,0)
       call multifab_build(deltaT_rhs2(n),mla%la(n),1,0)

       call multifab_build(    conc(n),mla%la(n),nspecies,rho_old(n)%ng)
       call multifab_build(molefrac(n),mla%la(n),nspecies,rho_old(n)%ng)

       call multifab_build(deltaT(n),mla%la(n),1,1)

       call multifab_build(cc_solver_alpha(n),mla%la(n),1,1)
       do i=1,dm
          call multifab_build_edge(cc_solver_beta(n,i),mla%la(n),1,1,i)
       end do

       call multifab_build(eta_old(n),mla%la(n),1,1)
       call multifab_build(eta_new(n),mla%la(n),1,1)

       call multifab_build(kappa_old(n),mla%la(n),1,1)
       call multifab_build(kappa_new(n),mla%la(n),1,1)

       call multifab_build(lambda_old(n),mla%la(n),1,1)
       call multifab_build(lambda_new(n),mla%la(n),1,1)

       call multifab_build(chi_old(n),mla%la(n),nspecies**2,1)
       call multifab_build(chi_new(n),mla%la(n),nspecies**2,1)

       call multifab_build(zeta_old(n),mla%la(n),nspecies,1)
       call multifab_build(zeta_new(n),mla%la(n),nspecies,1)

       call multifab_build(deltaS(n),mla%la(n),1,0)

       call multifab_build(     Scorr(n),mla%la(n),1,0)
       call multifab_build(deltaScorr(n),mla%la(n),1,0)

       call multifab_build(deltaalpha(n),mla%la(n),1,0)

       call multifab_build(Sproj(n),mla%la(n),1,0)
       call multifab_build(phi(n),mla%la(n),1,1)
       do i=1,dm
          call multifab_build_edge(rhoinv_fc(n,i),mla%la(n),1,1,i)
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 0a: Compute (S,alpha)^n, decompose (S,alpha,Scorr)^n, and 
    !          compute a pressure update
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute mass fractions in valid region and then fill ghost cells
    call convert_rho_to_conc(mla,rho_old,rhotot_old,conc,.true.)
    do n=1,nlevs
       ! fill ghost cells for two adjacent grids including periodic boundary ghost cells
       call multifab_fill_boundary(conc(n))
       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(conc(n),1,c_bc_comp,nspecies,the_bc_tower%bc_tower_array(n), &
                            dx_in=dx(n,:))
    end do

    ! compute mole fractions
    call convert_conc_to_molefrac(mla,conc,molefrac,.true.)

    ! compute initial transport properties
    call ideal_mixture_transport_wrapper(mla,rhotot_old,Temp,p0_old,conc,molefrac, &
                                         eta_old,lambda_old,kappa_old,chi_old,zeta_old)

    ! compute mass_fluxdiv = div(F^n)
    call mass_fluxdiv_energy(mla,rho_old,rhotot_old,molefrac,chi_old,zeta_old, &
                             gradp_baro,Temp,mass_fluxdiv,mass_flux,dx,the_bc_tower)

    ! compute rhoh_fluxdiv = div(Q)^n + sum(div(hk*Fk))^n + rho*Hext^n
    call rhoh_fluxdiv_energy(mla,lambda_old,Temp,mass_flux,rhotot_old,rhoh_fluxdiv, &
                             dx,time,the_bc_tower)

    ! compute S and alpha (store them in deltaS and deltaalpha)
    call compute_S_alpha(mla,deltaS,deltaalpha,mass_fluxdiv,rhoh_fluxdiv,conc, &
                         Temp,rhotot_old,p0_old)

    ! split S and alpha into average and perturbational pieces, e.g., (Sbar + deltaS)
    do n=1,nlevs
       Sbar     = multifab_sum_c(deltaS(n)    ,1,1)
       alphabar = multifab_sum_c(deltaalpha(n),1,1)
       call multifab_sub_sub_s_c(deltaS(n)    ,1,Sbar    ,1,0)
       call multifab_sub_sub_s_c(deltaalpha(n),1,alphabar,1,0)
    end do

    ! zero out volume discrepancy correction
    Scorrbar = 0.d0
    do n=1,nlevs
       call multifab_setval(Scorr(n)     ,0.d0,all=.true.)
       call multifab_setval(deltaScorr(n),0.d0,all=.true.)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 0b: Compute the velocity field using a projection
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 0c: Advance the densities using forward-Euler advective fluxes and 
    !          explicit mass diffusion
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 0d: Advance the enthalpy by iteratively looping over an energy solve
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 0d-1: Compute (lambda,cp,F)^{*,n+1,l}
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 0d-2: Solve for deltaT implicitly
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 0d-3: Update the temperature and enthalpy
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 0e: If the thermodynamic drift is unacceptable, update the volume
    !          discrepancy correction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    do n=1,nlevs
       call multifab_destroy(rho_update(n))
       do i=1,dm
          call multifab_destroy(mass_flux(n,i))
       end do
       call multifab_destroy(mass_fluxdiv(n))
       call multifab_destroy(rhoh_fluxdiv(n))
       call multifab_destroy(deltaT_rhs1(n))
       call multifab_destroy(deltaT_rhs2(n))
       call multifab_destroy(    conc(n))
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
       call multifab_destroy(deltaS(n))
       call multifab_destroy(     Scorr(n))
       call multifab_destroy(deltaScorr(n))
       call multifab_destroy(deltaalpha(n))
       call multifab_destroy(Sproj(n))
       call multifab_destroy(phi(n))
       do i=1,dm
          call multifab_destroy(rhoinv_fc(n,i))
       end do
    end do

  end subroutine initialize

end module initialize_module
