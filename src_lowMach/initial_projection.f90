module initial_projection_module

  use multifab_module
  use ml_layout_module
  use convert_stag_module
  use define_bc_module
  use macproject_module
  use div_and_grad_module
  use bc_module
  use multifab_physbc_stag_module
  use compute_mass_fluxdiv_module
  use stochastic_mass_fluxdiv_module
  use reservoir_bc_fill_module
  use fluid_charge_module
  use probin_common_module, only: rhobar, variance_coef_mass, algorithm_type, &
                                  molmass, nspecies, restart
  use probin_charged_module, only: dielectric_const, use_charged_fluid
  use probin_chemistry_module, only: nreactions, use_Poisson_rng
  use chemical_rates_module

  implicit none

  private

  public :: initial_projection

contains

  ! this routine is only called for all inertial simulations (both restart and non-restart)
  ! it does the following:
  ! 1. fill mass random numbers
  ! 2. computes mass fluxes and flux divergences
  ! if restarting, the subroutine ends; otherwise
  ! 3. perform an initial projection
  !
  ! overdamped schemes need to do 1. and 2. within the advance_timestep routine
  ! in principle, performing an initial projection for overdamped will change
  ! the reference state for the GMRES solver
  ! For overdamped the first ever solve cannot have a good reference state
  ! so in general there is the danger it will be less accurate than subsequent solves
  ! but I do not see how one can avoid that
  ! From this perspective it may be useful to keep initial_projection even in overdamped
  ! because different gmres tolerances may be needed in the first step than in the rest
  subroutine initial_projection(mla,umac,rho,rhotot,gradp_baro, &
                                Epot_mass_fluxdiv, &
                                diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                stoch_mass_flux,chem_rate, &
                                Temp,eta,eta_ed,dt,dx,the_bc_tower, &
                                charge_old,grad_Epot_old,Epot,permittivity)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(inout) :: rho(:)
    type(multifab) , intent(inout) :: rhotot(:)
    type(multifab) , intent(in   ) :: gradp_baro(:,:)
    type(multifab) , intent(inout) :: Epot_mass_fluxdiv(:)
    type(multifab) , intent(inout) :: diff_mass_fluxdiv(:)
    type(multifab) , intent(inout) :: stoch_mass_fluxdiv(:)
    type(multifab) , intent(inout) :: stoch_mass_flux(:,:)
    type(multifab) , intent(inout) :: chem_rate(:)
    type(multifab) , intent(in   ) :: Temp(:)
    type(multifab) , intent(in   ) :: eta(:)
    type(multifab) , intent(in   ) :: eta_ed(:,:)  ! nodal (2d); edge-centered (3d)
    real(kind=dp_t), intent(in   ) :: dt
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(inout) :: charge_old(:)
    type(multifab) , intent(inout) :: grad_Epot_old(:,:)
    type(multifab) , intent(inout) :: Epot(:)
    type(multifab) , intent(in   ) :: permittivity(:)

    ! local
    integer :: i,dm,n,nlevs
    real(kind=dp_t) :: dt_eff

    type(multifab) ::            mac_rhs(mla%nlevel)
    type(multifab) ::               divu(mla%nlevel)
    type(multifab) ::                phi(mla%nlevel)

    type(multifab) ::       rhotot_fc(mla%nlevel,mla%dim)
    type(multifab) ::    rhototinv_fc(mla%nlevel,mla%dim)
    type(multifab) ::  diff_mass_flux(mla%nlevel,mla%dim)
    type(multifab) :: total_mass_flux(mla%nlevel,mla%dim)


    type(multifab) :: n_cc(mla%nlevel)

    real(kind=dp_t), allocatable :: weights(:)
    
    type(bl_prof_timer), save :: bpt

    call build(bpt,"initial_projection")

    if (algorithm_type .eq. 2) then
       call bl_error("Should not call initial_projection for overdamped schemes")
    end if
    
    if (algorithm_type .eq. 5 .or. algorithm_type .eq. 6) then
       allocate(weights(2))
       weights(1) = 1.d0
       weights(2) = 0.d0
       ! for midpoint scheme where predictor goes to t^{n+1/2}
       dt_eff = 0.5d0*dt
    else
       allocate(weights(1))
       weights(1) = 1.d0
       ! predictor integrates over full time step
       dt_eff = dt
    end if

    if (nreactions > 0) then
       if (algorithm_type .ne. 5 .and. algorithm_type .ne. 6) then
          call bl_error('Error: only algorithm_type=(5 or 6) allowed for nreactions>0')
       else if (use_Poisson_rng .eq. 2) then
          call bl_error('Error: currently use_Poisson_rng=2 not allowed for algorithm_type=(5 or 6) and nreactions>0')
       end if
    end if

    dm = mla%dim
    nlevs = mla%nlevel
 
    do n=1,nlevs
       call multifab_build(mac_rhs(n),mla%la(n),1,0)
       call multifab_build(divu(n),mla%la(n),1,0)
       call multifab_build(phi(n),mla%la(n),1,1)
       do i=1,dm
          call multifab_build_edge(      rhotot_fc(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(   rhototinv_fc(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge( diff_mass_flux(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge(total_mass_flux(n,i),mla%la(n),nspecies,0,i)
       end do       
    end do

    if (nreactions > 0) then
       do n=1,nlevs
          call multifab_build(n_cc(n),mla%la(n),nspecies,0)
       end do
    end if

    ! reset inhomogeneous bc condition to deal with reservoirs
    call set_inhomogeneous_vel_bcs(mla,vel_bc_n,vel_bc_t,eta_ed,dx,0.d0, &
                                   the_bc_tower%bc_tower_array)

    ! compute diffusive, stochastic, and potential mass fluxes
    ! with barodiffusion and thermodiffusion
    ! this computes "-F = rho W chi [Gamma grad x... ]"
    if (use_charged_fluid) then
       do n=1,nlevs
          call multifab_setval(Epot_mass_fluxdiv(n),0.d0,all=.true.)
          call multifab_setval(Epot(n),0.d0,all=.true.)
       end do

    end if

    if (variance_coef_mass .ne. 0.d0) then
       call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)
    end if

    call compute_mass_fluxdiv(mla,rho,rhotot,gradp_baro,Temp, &
                              diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                              diff_mass_flux,stoch_mass_flux, &
                              dt_eff,0.d0,dx,weights,the_bc_tower, &
                              charge_old,grad_Epot_old,Epot,permittivity)

    ! assemble total fluxes to be used in reservoirs
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(total_mass_flux(n,i),1,diff_mass_flux(n,i),1,nspecies,0)
          if (variance_coef_mass .ne. 0.d0) then
             call multifab_plus_plus_c(total_mass_flux(n,i),1,stoch_mass_flux(n,i),1,nspecies,0)
          end if
       end do
    end do

    ! compute chemical rates m_i*R_i
    if (nreactions > 0) then
       ! convert rho (mass densities rho_i) into n_cc (number densities n_i=rho_i/m_i)
       do n=1,nlevs
          call multifab_copy_c(n_cc(n),1,rho(n),1,nspecies,0)
          do i=1,nspecies
             call multifab_div_div_s_c(n_cc(n),i,molmass(i),1,0)
          end do
       end do

       ! compute chemical rates R_i (units=[number density]/[time])
       call chemical_rates(mla,n_cc,chem_rate,dx,dt_eff)

       ! convert chemical rates R_i into m_i*R_i (units=[mass density]/[time])
       do n=1,nlevs
          do i=1,nspecies
             call multifab_mult_mult_s_c(chem_rate(n),i,molmass(i),1,0)
          end do
       end do
    end if

    ! set the Dirichlet velocity value on reservoir faces
    call reservoir_bc_fill(mla,total_mass_flux,vel_bc_n,the_bc_tower%bc_tower_array)

    if (restart .lt. 0) then

       ! project the velocities
       ! only for non-restarting runs
       call setval(mac_rhs(n),0.d0)

       if (algorithm_type .ne. 6) then

          ! set mac_rhs to -S = sum_i div(F_i)/rhobar_i
          do n=1,nlevs

             do i=1,nspecies
                call multifab_saxpy_3_cc(mac_rhs(n),1,-1.d0/rhobar(i),diff_mass_fluxdiv(n),i,1)
                if (use_charged_fluid) then
                   call multifab_saxpy_3_cc(mac_rhs(n),1,-1.d0/rhobar(i),Epot_mass_fluxdiv(n),i,1)
                end if
                if (variance_coef_mass .ne. 0.d0) then
                   call multifab_saxpy_3_cc(mac_rhs(n),1,-1.d0/rhobar(i),stoch_mass_fluxdiv(n),i,1)
                end if
                if (nreactions > 0) then
                   ! if nreactions>0, also add sum_i -(m_i*R_i)/rhobar_i
                   call multifab_saxpy_3_cc(mac_rhs(n),1,-1.d0/rhobar(i),chem_rate(n),i,1)
                end if
             end do
          end do

       end if

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! build rhs = div(v^init) - S^0
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       do n=1,nlevs
          do i=1,dm
             ! to deal with reservoirs
             ! set normal velocity on physical domain boundaries
             call multifab_physbc_domainvel(umac(n,i),vel_bc_comp+i-1, &
                                            the_bc_tower%bc_tower_array(n), &
                                            dx(n,:),vel_bc_n(n,:))
             ! fill periodic and interior ghost cells
             call multifab_fill_boundary(umac(n,i))
          end do
       end do

       ! set divu = div(v^init)
       call compute_div(mla,umac,divu,dx,1,1,1)

       ! add div(v^init) to mac_rhs
       ! now mac_rhs = div(v^init) - S
       do n=1,nlevs
          call multifab_plus_plus_c(mac_rhs(n),1,divu(n),1,1,0)
       end do

       ! average rhotot to faces
       call average_cc_to_face(nlevs,rhotot,rhotot_fc,1,scal_bc_comp,1,the_bc_tower%bc_tower_array)

       ! compute (1/rhotot) on faces
       do n=1,nlevs
          do i=1,dm
             call setval(rhototinv_fc(n,i),1.d0,all=.true.)
             call multifab_div_div_c(rhototinv_fc(n,i),1,rhotot_fc(n,i),1,1,0)
          end do
       end do

       ! solve div (1/rhotot) grad phi = div(v^init) - S^0
       ! solve to completion, i.e., use the 'full' solver
       do n=1,nlevs
          call setval(phi(n),0.d0)
       end do
       call macproject(mla,phi,umac,rhototinv_fc,mac_rhs,dx,the_bc_tower,.true.)

       ! v^0 = v^init - (1/rho^0) grad phi
       call subtract_weighted_gradp(mla,umac,rhototinv_fc,phi,dx,the_bc_tower)

       ! fill ghost cells
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

    end if

    deallocate(weights)

    do n=1,nlevs
       call multifab_destroy(mac_rhs(n))
       call multifab_destroy(divu(n))
       call multifab_destroy(phi(n))
       do i=1,dm
          call multifab_destroy(rhotot_fc(n,i))
          call multifab_destroy(rhototinv_fc(n,i))
          call multifab_destroy(diff_mass_flux(n,i))
          call multifab_destroy(total_mass_flux(n,i))
       end do
    end do

    if (nreactions > 0) then
       do n=1,nlevs
          call multifab_destroy(n_cc(n))
       end do
    end if

    call destroy(bpt)

  end subroutine initial_projection

end module initial_projection_module
