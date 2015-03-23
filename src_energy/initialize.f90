module initialize_module

  use ml_layout_module
  use multifab_physbc_module
  use bc_module
  use define_bc_module
  use energy_eos_module
  use energy_eos_wrapper_module
  use convert_variables_module
  use probin_multispecies_module, only: nspecies

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
                        gradp_baro,pi,eta,eta_ed,kappa,Temp,Temp_ed, &
                        diff_mass_fluxdiv,stoch_mass_fluxdiv, &
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

    ! local variables

    ! this will hold -div(rho*v)^n + div(F^n)
    type(multifab) :: rho_update(mla%nlevel)

    ! This will hold (rhoh)^n/dt - div(rhoh*v)^n + Sbar^n/alphabar^n 
    !                + (1/2)(div(Q^n) + sum(div(h_k^n F_k^n)) + (rho Hext)^n)
    ! for the RHS of the temperature diffusion solve.
    ! Each of these terms stays fixed over all l iterations.
    type(multifab) :: rhoh_update1(mla%nlevel)

    ! This will hold -(rho^{*,n+1}h^{*,n+1,l})/dt + Sbarcorr^n/alphabar^n
    !                + (1/2)(div(Q^{*,n+1,l}) + sum(div(h_k^{*,n+1,l}F_k^{*,n+1,l})
    !                + (rho Hext)^(*,n+1)
    ! for the RHS of the temperature diffusion solve.
    ! Each of these terms may change for each l iteration.
    type(multifab) :: rhoh_update2(mla%nlevel)

    type(multifab) :: conc_old(mla%nlevel)
    type(multifab) :: molefrac_old(mla%nlevel)

    type(multifab) :: deltaT(mla%nlevel)

    type(multifab) :: solver_alpha(mla%nlevel)
    type(multifab) :: solver_beta (mla%nlevel,mla%dim)

    type(multifab) :: lambda_old(mla%nlevel)
    type(multifab) :: lambda_new(mla%nlevel)

    type(multifab) :: S(mla%nlevel)
    type(multifab) :: deltaS(mla%nlevel)

    type(multifab) :: Scorr(mla%nlevel)
    type(multifab) :: deltaScorr(mla%nlevel)

    type(multifab) :: alpha(mla%nlevel)
    type(multifab) :: deltaalpha(mla%nlevel)

    integer :: Sbar, Scorrbar, alphabar

    integer :: n,nlevs,i,dm

    nlevs = mla%nlevel
    dm = mla%dim

    do n=1,nlevs
       call multifab_build(    conc_old(n),mla%la(n),nspecies,rho_old(n)%ng)
       call multifab_build(molefrac_old(n),mla%la(n),nspecies,rho_old(n)%ng)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 0a: Compute a pressure update
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute mass fractions in valid region and then fill ghost cells
    call convert_rho_to_conc(mla,rho_old,rhotot_old,conc_old,.true.)
    do n=1,nlevs
       ! fill ghost cells for two adjacent grids including periodic boundary ghost cells
       call multifab_fill_boundary(conc_old(n))
       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(conc_old(n),1,c_bc_comp,nspecies,the_bc_tower%bc_tower_array(n), &
                            dx_in=dx(n,:))
    end do

    ! compute mole fractions
    call convert_conc_to_molefrac(mla,conc_old,molefrac_old,.true.)

    ! compute initial transport properties
!    call ideal_mixture_transport(rhotot_old,Temp,p0_old,)

    ! Construct S.  Many pieces of S are used in later parts of the algorithm,
    ! e.g., density update or enthalpy solve, but with different scalings

    ! compute mass flux

    ! set enthalpy_update to div(Q)

    ! increment enthalpy_update by sum_k div (h_k F_k)

    ! increment enthalpy_update by rho*H_ext

    ! set rho_update to div(F_i) - deterministic part

    ! compute S

    ! compute alpha

    ! split S and alpha into average and perturbational pieces, e.g., (Sbar + deltaS)

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

  end subroutine initialize

end module initialize_module
