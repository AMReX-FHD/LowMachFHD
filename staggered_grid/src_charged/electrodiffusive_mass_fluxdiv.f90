module electrodiffusive_mass_fluxdiv_module

  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use div_and_grad_module
  use convert_stag_module
  use zero_edgeval_module
  use fluid_charge_module
  use bndry_reg_module
  use ml_solve_module
  use multifab_physbc_module
  use matvec_mul_module
  use probin_common_module, only: nspecies, variance_coef_mass, shift_cc_to_boundary, bc_lo, bc_hi
  use probin_charged_module, only: Epot_wall_bc_type, Epot_wall, E_ext_type, electroneutral, &
                                   zero_eps_on_wall_type, epot_mg_verbose, epot_mg_abs_tol, &
                                   epot_mg_rel_tol, charge_per_mass, relxn_param_charge
  use probin_multispecies_module, only: is_nonisothermal, use_multiphase
  
  use fabio_module

  implicit none

  private

  public :: electrodiffusive_mass_flux, electrodiffusive_mass_fluxdiv, &
       inhomogeneous_neumann_fix

contains

  subroutine electrodiffusive_mass_fluxdiv(mla,rho,Temp,rhoWchi, &
                                           diff_mass_flux,diff_mass_fluxdiv, &
                                           stoch_mass_flux, &
                                           dx,the_bc_tower,charge, &
                                           grad_Epot,Epot,permittivity,dt, &
                                           zero_initial_Epot)

    ! this adds -div(F) = div(A_Phi grad Phi) to diff_mass_fluxdiv
    ! grad_Epot = grad Phi

    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: rho(:)
    type(multifab) , intent(in   )  :: Temp(:)
    type(multifab) , intent(in   )  :: rhoWchi(:)
    type(multifab) , intent(inout)  :: diff_mass_flux(:,:)
    type(multifab) , intent(inout)  :: diff_mass_fluxdiv(:)
    type(multifab) , intent(in   )  :: stoch_mass_flux(:,:)
    real(kind=dp_t), intent(in   )  :: dx(:,:)
    type(bc_tower) , intent(in   )  :: the_bc_tower
    type(multifab) , intent(inout)  :: charge(:)
    type(multifab) , intent(inout)  :: grad_Epot(:,:)
    type(multifab) , intent(inout)  :: Epot(:)
    type(multifab) , intent(in   )  :: permittivity(:)
    real(kind=dp_t), intent(in   )  :: dt
    logical        , intent(in   )  :: zero_initial_Epot

    ! local variables
    integer i,dm,n,nlevs

    ! local array of multifabs for grad and div; one for each direction
    type(multifab) :: electro_mass_flux(mla%nlevel,mla%dim)
    
    type(bl_prof_timer), save :: bpt

    call build(bpt, "electrodiffusive_mass_fluxdiv")
    
    nlevs = mla%nlevel  ! number of levels 
    dm    = mla%dim     ! dimensionality
 
    ! build the local multifabs
    do n=1,nlevs
       do i=1,dm
          ! electro_mass_flux(i) is face-centered, has nspecies component, zero ghost 
          ! cells & nodal in direction i
          call multifab_build_edge(electro_mass_flux(n,i),mla%la(n),nspecies,0,i)
       end do
    end do   
    
    ! compute the face-centered electro_mass_flux (each direction: cells+1 faces while 
    ! cells contain interior+2 ghost cells) 
    call electrodiffusive_mass_flux(mla,rho,Temp,rhoWchi,electro_mass_flux, &
                                    diff_mass_flux,stoch_mass_flux,dx,the_bc_tower, &
                                    charge,grad_Epot,Epot,permittivity,dt,zero_initial_Epot)

    ! add fluxes to diff_mass_flux
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(diff_mass_flux(n,i),1,electro_mass_flux(n,i),1,nspecies,0)
       end do
    end do

    ! add flux divergence to diff_mass_fluxdiv
    call compute_div(mla,electro_mass_flux,diff_mass_fluxdiv,dx,1,1,nspecies, &
                     increment_in=.true.)
    
    ! destroy the multifab to free the memory
    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(electro_mass_flux(n,i))
       end do
    end do

    call destroy(bpt)

  end subroutine electrodiffusive_mass_fluxdiv
 
  subroutine electrodiffusive_mass_flux(mla,rho,Temp,rhoWchi,electro_mass_flux, &
                                        diff_mass_flux,stoch_mass_flux, &
                                        dx,the_bc_tower,charge,grad_Epot,Epot, &
                                        permittivity,dt,zero_initial_Epot)

    ! this computes "-F = A_Phi grad Phi"

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: rho(:) 
    type(multifab) , intent(in   ) :: Temp(:)  
    type(multifab) , intent(in   ) :: rhoWchi(:)
    type(multifab) , intent(inout) :: electro_mass_flux(:,:)
    type(multifab) , intent(in   ) ::    diff_mass_flux(:,:)
    type(multifab) , intent(in   ) ::   stoch_mass_flux(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(inout) :: charge(:)
    type(multifab) , intent(inout) :: grad_Epot(:,:)
    type(multifab) , intent(inout) :: Epot(:)
    type(multifab) , intent(in   ) :: permittivity(:)
    real(kind=dp_t), intent(in   ) :: dt
    logical        , intent(in   ) :: zero_initial_Epot

    ! local variables
    integer :: n,i,comp,dm,nlevs
 
    ! local face-centered multifabs 
    type(multifab)  :: rhoWchi_face(mla%nlevel,mla%dim)

    type(multifab)  ::       alpha(mla%nlevel)
    type(multifab)  :: charge_coef(mla%nlevel)
    type(multifab)  ::         rhs(mla%nlevel)
    type(multifab)  ::      rhsvec(mla%nlevel)

    type(multifab)  ::             beta(mla%nlevel,mla%dim)
    type(multifab)  ::  permittivity_fc(mla%nlevel,mla%dim)
    type(multifab)  ::            E_ext(mla%nlevel,mla%dim)
    type(multifab)  ::             A_Phi(mla%nlevel,mla%dim)

    type(multifab)  :: diffstoch_mass_flux(mla%nlevel,mla%dim)

    type(bndry_reg) :: fine_flx(mla%nlevel)
  
    real(kind=dp_t) :: sum

    real(kind=dp_t) :: Epot_wall_save(2,mla%dim)

    real(kind=dp_t) :: norm, epot_mg_abs_tol_temp

    type(bl_prof_timer), save :: bpt
    
    call build(bpt,"electrodiffusive_mass_flux")

    dm    = mla%dim     ! dimensionality
    nlevs = mla%nlevel  ! number of levels 
    
    do n=1,nlevs
       call multifab_build(      alpha(n),mla%la(n),       1,0)
       call multifab_build(charge_coef(n),mla%la(n),nspecies,1)
       call multifab_build(        rhs(n),mla%la(n),       1,0)
       call multifab_build( rhsvec(n),mla%la(n),nspecies,0)
       do i=1,dm
          call multifab_build_edge(            beta(n,i),mla%la(n),          1,0,i)
          call multifab_build_edge(    rhoWchi_face(n,i),mla%la(n),nspecies**2,0,i)
          call multifab_build_edge( permittivity_fc(n,i),mla%la(n),          1,0,i)
          call multifab_build_edge(           E_Ext(n,i),mla%la(n),          1,0,i)
          call multifab_build_edge(           A_Phi(n,i),mla%la(n),   nspecies,0,i)
          call multifab_build_edge( diffstoch_mass_flux(n,i),mla%la(n),   nspecies,0,i)
       end do
    end do

    ! build the boundary flux register
    do n=1,nlevs
       call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
    end do

    ! if periodic, ensure charge sums to zero by subtracting off the average
    if (all(mla%pmask(1:dm))) then
       sum = multifab_sum_c(charge(1),1,1) / multifab_volume(charge(1))
       call multifab_sub_sub_s_c(charge(1),1,sum,1,0)
       if (abs(sum) .gt. 1.d-12) then
          if (parallel_IOProcessor()) then
             print*,'average charge =',sum
          end if
          call bl_warn("Warning: electrodiffusive_mass_flux - average charge is not zero")
       end if
    end if

    ! compute face-centered rhoWchi from cell-centered values 
    if (any(shift_cc_to_boundary(:,:) .eq. 1)) then
       call shift_cc_to_boundary_face(nlevs, rhoWchi, rhoWchi_face, 1, tran_bc_comp, &
                                      nspecies**2, the_bc_tower%bc_tower_array, .false.) 
    else
       call average_cc_to_face(nlevs, rhoWchi, rhoWchi_face, 1, tran_bc_comp, &
                               nspecies**2, the_bc_tower%bc_tower_array, .false.) 
    end if

    ! solve poisson equation for phi (the electric potential)
    ! -del dot epsilon grad Phi = charge
    do n=1,nlevs

       if (zero_initial_Epot) then
          call setval(Epot(n),0.d0,all=.true.)
       end if

       ! fill ghost cells for Epot at walls using Dirichlet value
       call multifab_physbc(Epot(n),1,Epot_bc_comp,1,the_bc_tower%bc_tower_array(n), &
                            dx_in=dx(n,:))

       ! set alpha=0
       call setval(alpha(n),0.d0,all=.true.)

    end do

    if (.not. electroneutral .or. E_ext_type .ne. 0) then
       ! permittivity on faces
       if (any(shift_cc_to_boundary(:,:) .eq. 1)) then
          call shift_cc_to_boundary_face(nlevs,permittivity,permittivity_fc,1,scal_bc_comp,1, &
                                         the_bc_tower%bc_tower_array)
       else
          call average_cc_to_face(nlevs,permittivity,permittivity_fc,1,scal_bc_comp,1, &
                                  the_bc_tower%bc_tower_array)
       end if
    end if

    if (electroneutral) then

       ! For electroneutral we only support homogeneous Neumann BCs for potential
       ! This is the correct Poisson BC for impermeable walls
       ! For reservoirs, the BCs are actually inhomogeneous but computed on-the-fly by the code later on
       ! Here we setup just the homogeneous Poisson problem -- this is all that the multigrid solver can handle     
       ! Reactive walls are not yet supported
       
       ! check to make sure physical boundary use homogeneous Neumann conditions on electric potential
       do i=1,dm
          if (bc_lo(i) .ne. PERIODIC) then
             if (Epot_wall_bc_type(1,i) .eq. 1 .or. Epot_wall(1,i) .ne. 0.d0) then
                call bl_error("electroneutral algorithm requires homogeneous Neumann potential bc's on physical boundaries")
             end if
          end if
          if (bc_hi(i) .ne. PERIODIC) then
             if (Epot_wall_bc_type(2,i) .eq. 1 .or. Epot_wall(2,i) .ne. 0.d0) then
                call bl_error("electroneutral algorithm requires homogeneous Neumann potential bc's on physical boundaries")
             end if
          end if
       end do

       ! compute A_Phi for Poisson solve (does not have z^T)
       call implicit_potential_coef(mla,rho,Temp,A_Phi,the_bc_tower,rhoWchi_face)
              
       ! compute z^T A_Phi^n, store in solver_beta
       call dot_with_z_face(mla,A_Phi,beta)

       ! combine F_diffstoch = F_d + F_s
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(diffstoch_mass_flux(n,i),1,diff_mass_flux(n,i),1,nspecies,0)
             if (variance_coef_mass .ne. 0.d0) then
                call multifab_plus_plus_c(diffstoch_mass_flux(n,i),1,stoch_mass_flux(n,i),1,nspecies,0)
             end if
          end do
       end do
       
       ! zero F_diffstoch on all walls -- this is required for reservoirs
       ! This should not really be necessary except if there are thermodiffusive fluxes
       ! since the deterministic and stochastic diffusive fluxes at boundaries should be zero for walls
       do n=1,nlevs
          call zero_edgeval_physical(diffstoch_mass_flux(n,:),1,nspecies,the_bc_tower%bc_tower_array(n))
       end do

       ! compute RHS = div (z^T (F_d + F_s))
       ! first, set rhsvec = div (F_d + F_s)
       call compute_div(mla,diffstoch_mass_flux,rhsvec,dx,1,1,nspecies,increment_in=.false.)

       ! In order to prevent slow charge buildup, we include the charge density in the rhs
       ! We allow for a relaxation factor here, though the default value of 1 is fine
       ! increment rhsvec by rho; we will dot with z below
       do n=1,nlevs

           ! OPTION 2: add (rho w) / dt to RHS (dotted with z below)
          call multifab_saxpy_3(rhsvec(n),relxn_param_charge/dt,rho(n))  ! crashes

       end do

       !!!!!!!!!!!!!!!!!!!!!!
       ! change solver tolerance based on scales of the problem

       ! dot abs(z) with (div F)
       call dot_with_z(mla,rhsvec,rhs,abs_z=.true.)
       ! compute norm
       norm = multifab_norm_inf(rhs(1))

       ! set absolute tolerance to be the norm*epot_mg_rel_tol
       epot_mg_abs_tol_temp = epot_mg_abs_tol
       epot_mg_abs_tol = norm*epot_mg_rel_tol
       !!!!!!!!!!!!!!!!!!!!!!

       ! compute rhs for Poisson zolve, z^T (div F)
       call dot_with_z(mla,rhsvec,rhs)

       ! When including charge on the rhs, the sum may not be zero anymore though it should be very close to zero
       ! subtract off average of rhs to make system solvable
       sum = multifab_sum_c(rhs(1),1,1) / multifab_volume(rhs(1))
       call multifab_sub_sub_s_c(rhs(1),1,sum,1,0)

    else
       
       ! non-electroneutral

       ! set beta=permittivity (epsilon)
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(beta(n,i),1,permittivity_fc(n,i),1,1,0)
          end do
       end do

       if (zero_eps_on_wall_type .gt. 0) then
          ! set beta to set to zero on certain boundary faces
          call zero_eps_on_wall(mla,beta,dx)
       end if

       ! set rhs equal to charge
       do n=1,nlevs
          call multifab_copy_c(rhs(n),1,charge(n),1,1,0)
       end do

       ! for inhomogeneous Neumann bc's for electric potential, put in homogeneous form
       if ( any((Epot_wall_bc_type(1:2,1:dm) .eq. 2) .and. &
                (Epot_wall(1:2,1:dm) .ne. 0.d0     )) ) then

          ! save the numerical values for the Dirichlet and Neumann conditions
          Epot_wall_save(1:2,1:dm) = Epot_wall(1:2,1:dm)

          ! for Dirichlet conditions, temporarily set the numerical values to zero
          ! so we can put the Neumann boundaries into homogeneous form
          do comp=1,dm
             if (Epot_wall_bc_type(1,comp) .eq. 1) then
                Epot_wall(1,comp) = 0.d0
             end if
             if (Epot_wall_bc_type(2,comp) .eq. 1) then
                Epot_wall(2,comp) = 0.d0
             end if
          end do

          call inhomogeneous_neumann_fix(mla,rhs,permittivity,dx,the_bc_tower)

          ! restore the numerical values for the Dirichlet and Neumann conditions
          Epot_wall(1:2,1:dm) = Epot_wall_save(1:2,1:dm)

       end if

    end if

    if (E_ext_type .ne. 0) then

       ! compute external electric field on edges
       call compute_E_ext(mla,E_ext)

       ! compute epsilon*E_ext
       do n=1,nlevs
          do i=1,dm
             call multifab_mult_mult_c(permittivity_fc(n,i),1,E_ext(n,i),1,1,0)
          end do
       end do

       ! compute div (epsilon*E_ext) and SUBTRACT it to solver rhs
       ! this needs to be tested with spatially-varying E_ext OR epsilon
       do n=1,nlevs
          do i=1,dm
             call multifab_mult_mult_s(permittivity_fc(n,i),-1.d0,0)
          end do
       end do
       call compute_div(mla,permittivity_fc,rhs,dx,1,1,1,increment_in=.true.)
       do n=1,nlevs
          do i=1,dm
             call multifab_mult_mult_s(permittivity_fc(n,i),-1.d0,0)
          end do
       end do

    end if

    ! solve (alpha - del dot beta grad) Epot = charge (for electro-explicit)
    !   Inhomogeneous Dirichlet or homogeneous Neumann is OK
    ! solve (alpha - del dot beta grad) Epot = z^T F (for electro-neutral)
    !   Only homogeneous Neumann BCs supported
    call ml_cc_solve(mla,rhs,Epot,fine_flx,alpha,beta,dx(:,1:dm),the_bc_tower,Epot_bc_comp, &
                     eps=epot_mg_rel_tol, &
                     abs_eps=epot_mg_abs_tol, &
                     verbose=epot_mg_verbose, &
                     ok_to_fix_singular=.false.)

    ! restore original solver tolerance
    if (electroneutral) then
       epot_mg_abs_tol = epot_mg_abs_tol_temp
    end if

    ! for periodic problems subtract off the average of Epot
    ! we can generalize this later for walls
    if ( all(mla%pmask(1:dm)) ) then
       sum = multifab_sum(Epot(1)) / boxarray_dvolume(get_boxarray(Epot(1)))
       call multifab_sub_sub_s(Epot(1),sum)
    end if  

    ! fill ghost cells for electric potential
    if (electroneutral .and. &
         ( any(bc_lo(1:dm) .eq. NO_SLIP_RESERVOIR) .or. &
           any(bc_hi(1:dm) .eq. NO_SLIP_RESERVOIR) .or. &
           any(bc_lo(1:dm) .eq. SLIP_RESERVOIR)    .or. &
           any(bc_hi(1:dm) .eq. SLIP_RESERVOIR) ) ) then

       ! for electroneutral problems with reservoirs,
       ! the inhomogeneous BC for phi must be computed here for each face:
       ! grad(phi) = -z^T*F_diffstoch/(z^T*A_Phi)
       ! We only need this to fill in BCs for phi

       ! combine F_diffstoch = F_d + F_s
       ! Recompute this so the value on the boundary is the correct one
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(diffstoch_mass_flux(n,i),1,diff_mass_flux(n,i),1,nspecies,0)
             if (variance_coef_mass .ne. 0.d0) then
                call multifab_plus_plus_c(diffstoch_mass_flux(n,i),1,stoch_mass_flux(n,i),1,nspecies,0)
             end if
          end do
       end do

       ! fill ghost cells for phi using
       ! grad(phi) = -z^T*F_diffstoch/(z^T*A_Phi)
       do n=1,nlevs
          call fill_phi_bc_eln_reservoir(Epot(n), diffstoch_mass_flux(n,:), &
                                        A_Phi(n,:), beta(n,:), &
                                        1,nspecies,the_bc_tower%bc_tower_array(n),dx(n,:))
       end do


    else

       ! for all other problems, use the Dirichlet or Neumann values supplied by
       ! Epot_wall_bc_type and Epot_wall
       ! note that for the inhomogeneous Neumann phi case, since the solver assumed 
       ! homogeneous BC's, this routine will properly fill the ghost cells
       do n=1,nlevs
          call multifab_physbc(Epot(n),1,Epot_bc_comp,1,the_bc_tower%bc_tower_array(n), &
                               dx_in=dx(n,:))
          call multifab_fill_boundary(Epot(n))
       end do

    end if

    ! compute the gradient of the electric potential
    call compute_grad(mla,Epot,grad_Epot,dx,1,Epot_bc_comp,1,1,the_bc_tower%bc_tower_array)

    if (E_ext_type .ne. 0) then
       ! add external electric field
       ! since E = -grad(Epot), we have to subtract the external field from grad_Epot
       do n=1,nlevs
          do i=1,dm
             call multifab_sub_sub_c(grad_Epot(n,i),1,E_ext(n,i),1,1,0)
          end do
       end do
    end if

    if (zero_eps_on_wall_type .gt. 0) then
       ! Set E-field ie grad_Epot to be zero on certain boundary faces.
       ! This enforces dphi/dn = 0 on the parts of the (Dirichlet) wall we want. 
       call zero_eps_on_wall(mla,grad_Epot,dx)
    end if

    do n=1,nlevs
       do i=1,dm
          call multifab_fill_boundary(grad_Epot(n,i))
       end do
    end do

    ! compute the charge flux coefficient
    call compute_charge_coef(mla,rho,Temp,charge_coef)

    ! average charge flux coefficient to faces, store in flux
    if (any(shift_cc_to_boundary(:,:) .eq. 1)) then
       call shift_cc_to_boundary_face(nlevs,charge_coef,electro_mass_flux,1,c_bc_comp,nspecies, &
                                      the_bc_tower%bc_tower_array,.true.)
    else
       call average_cc_to_face(nlevs,charge_coef,electro_mass_flux,1,c_bc_comp,nspecies, &
                               the_bc_tower%bc_tower_array,.true.)
    end if

    ! multiply flux coefficient by gradient of electric potential
    do n=1,nlevs
       do i=1,dm
          do comp=1,nspecies
             call multifab_mult_mult_c(electro_mass_flux(n,i), comp, grad_Epot(n,i), 1, 1)
          end do
       end do
    end do

    if (use_multiphase) then
       call limit_emf(rho, electro_mass_flux, grad_Epot)
    end if

    ! compute -rhoWchi * (... ) on faces
    do n=1,nlevs
       do i=1,dm
          call matvec_mul(mla, electro_mass_flux(n,i), rhoWchi_face(n,i), nspecies)
       end do
    end do

    if (.false..and.electroneutral) then
       ! Print some fluxes for debugging
       
       do n=1,nlevs
          call print_flux_reservoir(electro_mass_flux(n,:), &
                                    grad_Epot(n,:), &
                                    A_Phi(n,:), &
                                    beta(n,:), &
                                    1,nspecies,the_bc_tower%bc_tower_array(n))
       end do

    end if

    ! for walls we need to zero the electro_mass_flux since we have already zero'd the diff and stoch
    ! mass fluxes.  For inhomogeneous Neumann conditions on Epot, the physically correct thing would
    ! have been to compute species gradients to exactly counterbalance the Neumann conditions on Epot
    ! so the total mass flux (diff + stoch + Epot) is zero, but the numerical remedy here is to simply
    ! zero them individually.
    do n=1,nlevs
       call zero_edgeval_walls(electro_mass_flux(n,:),1,nspecies, the_bc_tower%bc_tower_array(n))
    end do

    do n=1,nlevs
       call bndry_reg_destroy(fine_flx(n))
    end do

    do n=1,nlevs
       call multifab_destroy(alpha(n))
       call multifab_destroy(charge_coef(n))
       call multifab_destroy(rhs(n))
       call multifab_destroy(rhsvec(n))
       do i=1,dm
          call multifab_destroy(rhoWchi_face(n,i))
          call multifab_destroy(beta(n,i))
          call multifab_destroy(permittivity_fc(n,i))
          call multifab_destroy(E_ext(n,i))
          call multifab_destroy(A_Phi(n,i))
          call multifab_destroy(diffstoch_mass_flux(n,i))
       end do
    end do

    call destroy(bpt)

  end subroutine electrodiffusive_mass_flux
  
  ! We would like to solve A x = b with inhomogeneous bc's
  ! Here, "A" is -div epsilon grad
  ! This is equivalent to A_H x = b - A x_H, where
  !   A   is the inhomogeneous operator
  !   A_H is the homogeneous operator
  !   x_H is a multifab filled with zeros, but ghost cells filled to respect bc's
  ! We use this for walls with inhomogeneous Neumann conditions on the electric potential
  subroutine inhomogeneous_neumann_fix(mla,rhs,permittivity,dx,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rhs(:)
    type(multifab) , intent(in   ) :: permittivity(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    type(multifab) :: zerofab(mla%nlevel)
    type(multifab) :: gradphi(mla%nlevel,mla%dim)
    integer :: i,dm,n,nlevs

    dm    = mla%dim     ! dimensionality
    nlevs = mla%nlevel  ! number of levels 

    do n=1,nlevs
       call multifab_build(zerofab(n),mla%la(n),1,1)
       do i=1,dm
          call multifab_build_edge(gradphi(n,i),mla%la(n),1,0,i)
       end do
    end do

    do n=1,nlevs
       call multifab_setval(zerofab(n),0.d0,all=.true.)
    end do

    do n=1,nlevs

       ! fill ghost cells for zerofab
       call multifab_physbc(zerofab(n),1,Epot_bc_comp,1,the_bc_tower%bc_tower_array(n), &
                            dx_in=dx(n,:))

       call multifab_fill_boundary(zerofab(n))

       ! multiply zerofab everywhere (including ghost cells) by the permittivity since
       ! we are incrementing the RHS (rhs) by -A x_H
       call multifab_mult_mult_c(zerofab(n),1,permittivity(n),1,1,1)

    end do

    ! compute gradient of zerofab
    call compute_grad(mla,zerofab,gradphi,dx,1,Epot_bc_comp,1,1,the_bc_tower%bc_tower_array)

    ! increment rhs with negative divergence
    call compute_div(mla,gradphi,rhs,dx,1,1,1,increment_in=.true.)

    do n=1,nlevs
       call multifab_fill_boundary(rhs(n))
    end do

    do n=1,nlevs
       call multifab_destroy(zerofab(n))
       do i=1,dm
          call multifab_destroy(gradphi(n,i))
       end do
    end do

  end subroutine inhomogeneous_neumann_fix


  ! The inhomogeneous BC for phi must be computed for electroneutral reservoirs
  ! grad(phi) = -z^T*F_diffstoch/(z^T*A_Phi)
  subroutine fill_phi_bc_eln_reservoir(Epot,diffstoch_mass_flux,A_phi,z_dot_A, &
                                       start_comp,num_comp,the_bc_level,dx)


    type(multifab) , intent(inout) :: Epot
    type(multifab) , intent(in)    :: diffstoch_mass_flux(:)
    type(multifab) , intent(in)    :: A_phi(:)
    type(multifab) , intent(in)    :: z_dot_A(:)
    integer        , intent(in   ) :: start_comp,num_comp
    type(bc_level) , intent(in   ) :: the_bc_level
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! Local
    integer                  :: lo(get_dim(Epot)),hi(get_dim(Epot))
    integer                  :: ng_e,ng_d,ng_a,ng_z,i,dm
    real(kind=dp_t), pointer :: ep(:,:,:,:)
    real(kind=dp_t), pointer :: dpx(:,:,:,:), dpy(:,:,:,:), dpz(:,:,:,:)
    real(kind=dp_t), pointer :: apx(:,:,:,:), apy(:,:,:,:), apz(:,:,:,:)
    real(kind=dp_t), pointer :: zpx(:,:,:,:), zpy(:,:,:,:), zpz(:,:,:,:)

    dm = get_dim(diffstoch_mass_flux(1))
    ng_e = nghost(Epot)
    ng_d = nghost(diffstoch_mass_flux(1))
    ng_a = nghost(A_phi(1))
    ng_z = nghost(z_dot_A(1))
    
    do i=1,nfabs(Epot)
       ep  => dataptr(Epot,i)
       dpx => dataptr(diffstoch_mass_flux(1),i)
       dpy => dataptr(diffstoch_mass_flux(2),i)
       apx => dataptr(A_Phi(1),i)
       apy => dataptr(A_Phi(2),i)
       zpx => dataptr(z_dot_A(1),i)
       zpy => dataptr(z_dot_A(2),i)

       lo = lwb(get_box(Epot,i))
       hi = upb(get_box(Epot,i))
       select case (dm)
       case (2)
          call fill_phi_bc_eln_reservoir_2d(ep(:,:,1,1), ng_e, &
                                            dpx(:,:,1,start_comp:start_comp+num_comp-1), &
                                            dpy(:,:,1,start_comp:start_comp+num_comp-1), ng_d, &
                                            apx(:,:,1,start_comp:start_comp+num_comp-1), &
                                            apy(:,:,1,start_comp:start_comp+num_comp-1), ng_a, &
                                            zpx(:,:,1,1), zpy(:,:,1,1), ng_z, &
                                            lo, hi, the_bc_level%phys_bc_level_array(i,:,:), dx)
       case (3)
          dpz => dataptr(diffstoch_mass_flux(3),i)
          apz => dataptr(A_Phi(3),i)
          zpz => dataptr(z_dot_A(3),i)
          call fill_phi_bc_eln_reservoir_3d(ep(:,:,:,1), ng_e, &
                                            dpx(:,:,:,start_comp:start_comp+num_comp-1), &
                                            dpy(:,:,:,start_comp:start_comp+num_comp-1), &
                                            dpz(:,:,:,start_comp:start_comp+num_comp-1), ng_d, &
                                            apx(:,:,:,start_comp:start_comp+num_comp-1), &
                                            apy(:,:,:,start_comp:start_comp+num_comp-1), &
                                            apz(:,:,:,start_comp:start_comp+num_comp-1), ng_a, &
                                            zpx(:,:,:,1), zpy(:,:,:,1), zpz(:,:,:,1), ng_z, &
                                            lo, hi, the_bc_level%phys_bc_level_array(i,:,:), dx)
       end select
    end do
 
  end subroutine fill_phi_bc_eln_reservoir

  subroutine fill_phi_bc_eln_reservoir_2d(Epot,ng_e, &
                                          diffstoch_mass_fluxx,diffstoch_mass_fluxy,ng_d, &
                                          A_Phix,A_Phiy,ng_a, &
                                          z_dot_Ax,z_dot_Ay,ng_z, &
                                          lo,hi,bc,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng_e,ng_d,ng_a,ng_z
    real(kind=dp_t), intent(inout) ::                 Epot(lo(1)-ng_e:,lo(2)-ng_e:)
    real(kind=dp_t), intent(inout) :: diffstoch_mass_fluxx(lo(1)-ng_d:,lo(2)-ng_d:,:)
    real(kind=dp_t), intent(inout) :: diffstoch_mass_fluxy(lo(1)-ng_d:,lo(2)-ng_d:,:)
    real(kind=dp_t), intent(inout) ::               A_Phix(lo(1)-ng_a:,lo(2)-ng_a:,:)
    real(kind=dp_t), intent(inout) ::               A_Phiy(lo(1)-ng_a:,lo(2)-ng_a:,:)
    real(kind=dp_t), intent(inout) ::             z_dot_Ax(lo(1)-ng_z:,lo(2)-ng_z:)
    real(kind=dp_t), intent(inout) ::             z_dot_Ay(lo(1)-ng_z:,lo(2)-ng_z:)
    integer        , intent(in   ) :: bc(:,:)
    real(kind=dp_t)                :: dx(:)    
    
    integer :: i,j
    real(kind=dp_t) :: zTF, gradphi

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. NO_SLIP_RESERVOIR .or. bc(1,1) .eq. SLIP_RESERVOIR) then
       do j=lo(2),hi(2)
          zTF=dot_product(charge_per_mass(1:nspecies), diffstoch_mass_fluxx(lo(1),j,1:nspecies))
          gradphi= -zTF / dot_product(charge_per_mass(1:nspecies), A_Phix(lo(1),j,1:nspecies))
          Epot(lo(1)-1,j) = Epot(lo(1),j) - gradphi*dx(1)/2.d0
       end do
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. NO_SLIP_RESERVOIR .or. bc(1,2) .eq. SLIP_RESERVOIR) then
       do j=lo(2),hi(2)
          zTF=dot_product(charge_per_mass(1:nspecies), diffstoch_mass_fluxx(hi(1)+1,j,1:nspecies))
          gradphi= -zTF / dot_product(charge_per_mass(1:nspecies), A_Phix(hi(1)+1,j,1:nspecies))
          Epot(hi(1)+1,j) = Epot(hi(1),j) + gradphi*dx(1)/2.d0
       end do
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. NO_SLIP_RESERVOIR .or. bc(2,1) .eq. SLIP_RESERVOIR) then
       do i=lo(1),hi(1)
          zTF=dot_product(charge_per_mass(1:nspecies), diffstoch_mass_fluxy(i,lo(2),1:nspecies))
          gradphi= -zTF / dot_product(charge_per_mass(1:nspecies), A_Phiy(i,lo(2),1:nspecies))
          Epot(i,lo(2)-1) = Epot(i,lo(2)) - gradphi*dx(2)/2.d0
       end do
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. NO_SLIP_RESERVOIR .or. bc(2,2) .eq. SLIP_RESERVOIR) then
       do i=lo(1),hi(1)
          zTF=dot_product(charge_per_mass(1:nspecies), diffstoch_mass_fluxy(i,hi(2)+1,1:nspecies))
          gradphi= -zTF / dot_product(charge_per_mass(1:nspecies), A_Phiy(i,hi(2)+1,1:nspecies))
          Epot(i,hi(2)+1) = Epot(i,hi(2)) + gradphi*dx(2)/2.d0
       end do
    end if

  end subroutine fill_phi_bc_eln_reservoir_2d

  subroutine fill_phi_bc_eln_reservoir_3d(Epot,ng_e, &
                                          diffstoch_mass_fluxx, &
                                          diffstoch_mass_fluxy, &
                                          diffstoch_mass_fluxz,ng_d, &
                                          A_Phix,A_Phiy,A_Phiz,ng_a, &
                                          z_dot_Ax,z_dot_Ay,z_dot_Az,ng_z, &
                                          lo,hi,bc,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng_e,ng_d,ng_a,ng_z
    real(kind=dp_t), intent(inout) ::                 Epot(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(inout) :: diffstoch_mass_fluxx(lo(1)-ng_d:,lo(2)-ng_d:,lo(3)-ng_d:,:)
    real(kind=dp_t), intent(inout) :: diffstoch_mass_fluxy(lo(1)-ng_d:,lo(2)-ng_d:,lo(3)-ng_d:,:)
    real(kind=dp_t), intent(inout) :: diffstoch_mass_fluxz(lo(1)-ng_d:,lo(2)-ng_d:,lo(3)-ng_d:,:)
    real(kind=dp_t), intent(inout) ::               A_Phix(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:,:)
    real(kind=dp_t), intent(inout) ::               A_Phiy(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:,:)
    real(kind=dp_t), intent(inout) ::               A_Phiz(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:,:)
    real(kind=dp_t), intent(inout) ::             z_dot_Ax(lo(1)-ng_z:,lo(2)-ng_z:,lo(3)-ng_z:)
    real(kind=dp_t), intent(inout) ::             z_dot_Ay(lo(1)-ng_z:,lo(2)-ng_z:,lo(3)-ng_z:)
    real(kind=dp_t), intent(inout) ::             z_dot_Az(lo(1)-ng_z:,lo(2)-ng_z:,lo(3)-ng_z:)
    integer        , intent(in   ) :: bc(:,:)
    real(kind=dp_t)                :: dx(:)    
    
    integer :: i,j,k
    real(kind=dp_t) :: zTF, gradphi

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. NO_SLIP_RESERVOIR .or. bc(1,1) .eq. SLIP_RESERVOIR) then
       do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          zTF=dot_product(charge_per_mass(1:nspecies), diffstoch_mass_fluxx(lo(1),j,k,1:nspecies))
          gradphi= -zTF / dot_product(charge_per_mass(1:nspecies), A_Phix(lo(1),j,k,1:nspecies))
          Epot(lo(1)-1,j,k) = Epot(lo(1),j,k) - gradphi*dx(1)/2.d0
       end do
       end do
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. NO_SLIP_RESERVOIR .or. bc(1,2) .eq. SLIP_RESERVOIR) then
       do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          zTF=dot_product(charge_per_mass(1:nspecies), diffstoch_mass_fluxx(hi(1)+1,j,k,1:nspecies))
          gradphi= -zTF / dot_product(charge_per_mass(1:nspecies), A_Phix(hi(1)+1,j,k,1:nspecies))
          Epot(hi(1)+1,j,k) = Epot(hi(1),j,k) + gradphi*dx(1)/2.d0
       end do
       end do
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. NO_SLIP_RESERVOIR .or. bc(2,1) .eq. SLIP_RESERVOIR) then
       do k=lo(3),hi(3)
       do i=lo(1),hi(1)
          zTF=dot_product(charge_per_mass(1:nspecies), diffstoch_mass_fluxy(i,lo(2),k,1:nspecies))
          gradphi= -zTF / dot_product(charge_per_mass(1:nspecies), A_Phiy(i,lo(2),k,1:nspecies))
          Epot(i,lo(2)-1,k) = Epot(i,lo(2),k) - gradphi*dx(2)/2.d0
       end do
       end do
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. NO_SLIP_RESERVOIR .or. bc(2,2) .eq. SLIP_RESERVOIR) then
       do k=lo(3),hi(3)
       do i=lo(1),hi(1)
          zTF=dot_product(charge_per_mass(1:nspecies), diffstoch_mass_fluxy(i,hi(2)+1,k,1:nspecies))
          gradphi= -zTF / dot_product(charge_per_mass(1:nspecies), A_Phiy(i,hi(2)+1,k,1:nspecies))
          Epot(i,hi(2)+1,k) = Epot(i,hi(2),k) + gradphi*dx(2)/2.d0
       end do
       end do
    end if

!!!!!!!!!!!!!!!!!!
! lo-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,1) .eq. NO_SLIP_RESERVOIR .or. bc(3,1) .eq. SLIP_RESERVOIR) then
       do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          zTF=dot_product(charge_per_mass(1:nspecies), diffstoch_mass_fluxz(i,j,lo(3),1:nspecies))
          gradphi= -zTF / dot_product(charge_per_mass(1:nspecies), A_Phiz(i,j,lo(3),1:nspecies))
          Epot(i,j,lo(3)-1) = Epot(i,j,lo(3)) - gradphi*dx(3)/2.d0
       end do
       end do
    end if

!!!!!!!!!!!!!!!!!!
! hi-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,2) .eq. NO_SLIP_RESERVOIR .or. bc(3,2) .eq. SLIP_RESERVOIR) then
       do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          zTF=dot_product(charge_per_mass(1:nspecies), diffstoch_mass_fluxz(i,j,hi(3)+1,1:nspecies))
          gradphi= -zTF / dot_product(charge_per_mass(1:nspecies), A_Phiz(i,j,hi(3)+1,1:nspecies))
          Epot(i,j,hi(3)+1) = Epot(i,j,hi(3)) + gradphi*dx(3)/2.d0
       end do
       end do
    end if

  end subroutine fill_phi_bc_eln_reservoir_3d

!=================================
! This is old code only used here to print out the values of the electrodiffusive fluxes
! and compare to the new code
  subroutine print_flux_reservoir(electro_mass_flux,grad_Epot,A_phi,z_dot_A, &
                                    start_comp,num_comp,the_bc_level)

    ! sets the electro_mass_flux at reservoirs equal to  F_e = - A_Phi (z^T F_diffstoch ) / (z^T A_Phi)
    ! This ensures that z^T*(F_diffstoch+F_e) = z^T*(F_diffstoch + A_phi*grad(phi)) = 0 on the reservoir walls

    type(multifab) , intent(inout) :: electro_mass_flux(:)
    type(multifab) , intent(inout) :: grad_Epot(:)
    type(multifab) , intent(inout) :: A_phi(:)
    type(multifab) , intent(inout) :: z_dot_A(:)
    integer        , intent(in   ) :: start_comp,num_comp
    type(bc_level) , intent(in   ) :: the_bc_level

    ! Local
    integer                  :: lo(get_dim(electro_mass_flux(1))),hi(get_dim(electro_mass_flux(1)))
    integer                  :: ng_e,ng_d,ng_a,ng_z,i,dm,comp
    real(kind=dp_t), pointer :: epx(:,:,:,:), epy(:,:,:,:), epz(:,:,:,:)
    real(kind=dp_t), pointer :: dpx(:,:,:,:), dpy(:,:,:,:), dpz(:,:,:,:)
    real(kind=dp_t), pointer :: apx(:,:,:,:), apy(:,:,:,:), apz(:,:,:,:)
    real(kind=dp_t), pointer :: zpx(:,:,:,:), zpy(:,:,:,:), zpz(:,:,:,:)

    dm = get_dim(electro_mass_flux(1))
    ng_e = nghost(electro_mass_flux(1))
    ng_d = nghost(grad_Epot(1))
    ng_a = nghost(A_phi(1))
    ng_z = nghost(z_dot_A(1))
    
    do i=1,nfabs(electro_mass_flux(1))
       epx => dataptr(electro_mass_flux(1),i)
       epy => dataptr(electro_mass_flux(2),i)
       dpx => dataptr(grad_Epot(1),i)
       dpy => dataptr(grad_Epot(2),i)
       apx => dataptr(A_Phi(1),i)
       apy => dataptr(A_Phi(2),i)
       zpx => dataptr(z_dot_A(1),i)
       zpy => dataptr(z_dot_A(2),i)

       lo = lwb(get_box(electro_mass_flux(1),i))
       hi = upb(get_box(electro_mass_flux(1),i))
       write(*,*) "Starting to print reservoir fluxes"
       do comp=start_comp,start_comp+num_comp-1
          write(*,*) "---------------------------"
          write(*,*) "species=", comp-start_comp+1  
          select case (dm)
          case (2)
             call print_flux_reservoir_2d(epx(:,:,1,comp), epy(:,:,1,comp), ng_e, &
                                          dpx(:,:,1,1), dpy(:,:,1,1), ng_d, &
                                          apx(:,:,1,comp), apy(:,:,1,comp), ng_a, &
                                          zpx(:,:,1,1), zpy(:,:,1,1), ng_z, &
                                          lo, hi, &
                                          the_bc_level%phys_bc_level_array(i,:,:))
          case (3)
             epz => dataptr(electro_mass_flux(3),i)
             dpz => dataptr(grad_Epot(3),i)
             apz => dataptr(A_Phi(3),i)
             zpz => dataptr(z_dot_A(3),i)
             call print_flux_reservoir_3d(epx(:,:,:,comp), epy(:,:,:,comp), epz(:,:,:,comp), ng_e, &
                                          dpx(:,:,:,1), dpy(:,:,:,1), dpz(:,:,:,1), ng_d, &
                                          apx(:,:,:,comp), apy(:,:,:,comp), apz(:,:,:,comp), ng_a, &
                                          zpx(:,:,:,1), zpy(:,:,:,1), zpz(:,:,:,1), ng_z, &
                                          lo, hi, &
                                          the_bc_level%phys_bc_level_array(i,:,:))
          end select
       end do
       write(*,*) "---------------------------"
    end do
 
  end subroutine print_flux_reservoir

  subroutine print_flux_reservoir_2d(electro_mass_fluxx,electro_mass_fluxy,ng_e, &
                                     grad_Epotx,grad_Epoty,ng_d, &
                                     A_Phix,A_Phiy,ng_a, &
                                     z_dot_Ax,z_dot_Ay,ng_z, &
                                     lo,hi,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng_e,ng_d,ng_a,ng_z
    real(kind=dp_t), intent(inout) :: electro_mass_fluxx(lo(1)-ng_e:,lo(2)-ng_e:)
    real(kind=dp_t), intent(inout) :: electro_mass_fluxy(lo(1)-ng_e:,lo(2)-ng_e:)
    real(kind=dp_t), intent(inout) ::         grad_Epotx(lo(1)-ng_d:,lo(2)-ng_d:)
    real(kind=dp_t), intent(inout) ::         grad_Epoty(lo(1)-ng_d:,lo(2)-ng_d:)
    real(kind=dp_t), intent(inout) ::             A_Phix(lo(1)-ng_a:,lo(2)-ng_a:)
    real(kind=dp_t), intent(inout) ::             A_Phiy(lo(1)-ng_a:,lo(2)-ng_a:)
    real(kind=dp_t), intent(inout) ::           z_dot_Ax(lo(1)-ng_z:,lo(2)-ng_z:)
    real(kind=dp_t), intent(inout) ::           z_dot_Ay(lo(1)-ng_z:,lo(2)-ng_z:)
    integer        , intent(in   ) :: bc(:,:)
    
    integer :: i,j
    real(kind=dp_t) :: zTF

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. NO_SLIP_RESERVOIR .or. bc(1,1) .eq. SLIP_RESERVOIR) then
       do j=lo(2),hi(2)
          write(*,*) "LO-X grad_phi=", grad_Epotx(lo(1),j),"F_el=", electro_mass_fluxx(lo(1),j)
       end do
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. NO_SLIP_RESERVOIR .or. bc(1,2) .eq. SLIP_RESERVOIR) then
       do j=lo(2),hi(2)
          write(*,*) "HI-X grad_phi=", grad_Epotx(hi(1)+1,j),"F_el=", electro_mass_fluxx(hi(1)+1,j)
       end do
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. NO_SLIP_RESERVOIR .or. bc(2,1) .eq. SLIP_RESERVOIR) then
       do i=lo(1),hi(1) 
          write(*,*) "LO-Y grad_phi=", grad_Epoty(i,lo(2)),"F_el=", electro_mass_fluxy(i,lo(2))
       end do
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. NO_SLIP_RESERVOIR .or. bc(2,2) .eq. SLIP_RESERVOIR) then
       do i=lo(1),hi(1)
          write(*,*) "HI-Y grad_phi=", grad_Epoty(i,hi(2)+1),"F_el=", electro_mass_fluxy(i,hi(2)+1)
       end do
    end if

  end subroutine print_flux_reservoir_2d

  subroutine print_flux_reservoir_3d(electro_mass_fluxx,electro_mass_fluxy,electro_mass_fluxz,ng_e, &
                                     grad_Epotx,grad_Epoty,grad_Epotz,ng_d, &
                                     A_Phix,A_Phiy,A_Phiz,ng_a, &
                                     z_dot_Ax,z_dot_Ay,z_dot_Az,ng_z, &
                                     lo,hi,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng_e,ng_d,ng_a,ng_z
    real(kind=dp_t), intent(inout) :: electro_mass_fluxx(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(inout) :: electro_mass_fluxy(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(inout) :: electro_mass_fluxz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(inout) ::         grad_Epotx(lo(1)-ng_d:,lo(2)-ng_d:,lo(3)-ng_d:)
    real(kind=dp_t), intent(inout) ::         grad_Epoty(lo(1)-ng_d:,lo(2)-ng_d:,lo(3)-ng_d:)
    real(kind=dp_t), intent(inout) ::         grad_Epotz(lo(1)-ng_d:,lo(2)-ng_d:,lo(3)-ng_d:)
    real(kind=dp_t), intent(inout) ::             A_Phix(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
    real(kind=dp_t), intent(inout) ::             A_Phiy(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
    real(kind=dp_t), intent(inout) ::             A_Phiz(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
    real(kind=dp_t), intent(inout) ::           z_dot_Ax(lo(1)-ng_z:,lo(2)-ng_z:,lo(3)-ng_z:)
    real(kind=dp_t), intent(inout) ::           z_dot_Ay(lo(1)-ng_z:,lo(2)-ng_z:,lo(3)-ng_z:)
    real(kind=dp_t), intent(inout) ::           z_dot_Az(lo(1)-ng_z:,lo(2)-ng_z:,lo(3)-ng_z:)
    integer        , intent(in   ) :: bc(:,:)
    
    integer :: i,j,k
    real(kind=dp_t) :: zTF

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. NO_SLIP_RESERVOIR .or. bc(1,1) .eq. SLIP_RESERVOIR) then
       do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          write(*,*) "LO-X grad_phi=", grad_Epotx(lo(1),j,k),"F_el=", electro_mass_fluxx(lo(1),j,k)
       end do
       end do
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. NO_SLIP_RESERVOIR .or. bc(1,2) .eq. SLIP_RESERVOIR) then
       do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          write(*,*) "HI-X grad_phi=", grad_Epotx(hi(1)+1,j,k),"F_el=", electro_mass_fluxx(hi(1)+1,j,k)
       end do
       end do
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. NO_SLIP_RESERVOIR .or. bc(2,1) .eq. SLIP_RESERVOIR) then
       do k=lo(3),hi(3)
       do i=lo(1),hi(1)
          write(*,*) "LO-Y grad_phi=", grad_Epoty(i,lo(2),k),"F_el=", electro_mass_fluxy(i,lo(2),k)
       end do
       end do
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. NO_SLIP_RESERVOIR .or. bc(2,2) .eq. SLIP_RESERVOIR) then
       do k=lo(3),hi(3)
       do i=lo(1),hi(1)
          write(*,*) "HI-Y grad_phi=", grad_Epoty(i,hi(2)+1,k),"F_el=", electro_mass_fluxy(i,hi(2)+1,k)
       end do
       end do
    end if

!!!!!!!!!!!!!!!!!!
! lo-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,1) .eq. NO_SLIP_RESERVOIR .or. bc(3,1) .eq. SLIP_RESERVOIR) then
       do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          write(*,*) "LO-Z grad_phi=", grad_Epotz(i,j,lo(3)),"F_el=", electro_mass_fluxz(i,j,lo(3))
       end do
       end do
    end if

!!!!!!!!!!!!!!!!!!
! hi-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,2) .eq. NO_SLIP_RESERVOIR .or. bc(3,2) .eq. SLIP_RESERVOIR) then
       do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          write(*,*) "HI-Z grad_phi=", grad_Epotz(i,j,hi(3)+1),"F_el=", electro_mass_fluxz(i,j,hi(3)+1)
       end do
       end do
    end if

  end subroutine print_flux_reservoir_3d

  ! 0 emf for negative densities
  subroutine limit_emf(rho, electro_mass_flux, grad_Epot)

    type(multifab) , intent(inout) :: electro_mass_flux(:,:)
    type(multifab) , intent(in) :: rho(:)
    type(multifab) , intent(in) :: grad_Epot(:,:)

    ! Local
    integer                  :: lo(get_dim(rho(1))),hi(get_dim(rho(1)))
    integer                  :: ng_r,ng_e,ng_g,i,dm
    real(kind=dp_t), pointer :: rp(:,:,:,:)
    real(kind=dp_t), pointer :: epx(:,:,:,:), epy(:,:,:,:), epz(:,:,:,:)
    real(kind=dp_t), pointer :: gpx(:,:,:,:), gpy(:,:,:,:), gpz(:,:,:,:)

    dm = get_dim(rho(1))
    ng_r = nghost(rho(1))
    ng_e = nghost(electro_mass_flux(1,1))
    ng_g = nghost(grad_Epot(1,1))
    
    do i=1,nfabs(rho(1))
        rp => dataptr(rho(1),i)
       epx => dataptr(electro_mass_flux(1,1),i)
       epy => dataptr(electro_mass_flux(1,2),i)
       gpx => dataptr(grad_Epot(1,1),i)
       gpy => dataptr(grad_Epot(1,2),i)

       lo = lwb(get_box(rho(1),i))
       hi = upb(get_box(rho(1),i))
          select case (dm)
          case (2)
             call limit_emf_2d(rp(:,:,1,:), ng_r, &
                                          epx(:,:,1,:), epy(:,:,1,:), ng_e, &
                                          gpx(:,:,1,1), gpy(:,:,1,1), ng_g, &
                                          lo, hi)
          case (3)
             epz => dataptr(electro_mass_flux(1,3),i)
             gpz => dataptr(grad_Epot(1,3),i)
             call limit_emf_3d(rp(:,:,:,:), ng_r, &
                                          epx(:,:,:,:), epy(:,:,:,:), epz(:,:,:,:), ng_e, &
                                          gpx(:,:,:,1), gpy(:,:,:,1), gpz(:,:,:,1), ng_g, &
                                          lo, hi)
          end select
    end do
 
  end subroutine limit_emf


  subroutine limit_emf_2d(rho, ng_r,electro_mass_fluxx,electro_mass_fluxy,ng_e, &
                                     grad_Epotx,grad_Epoty,ng_g, &
                                     lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:),ng_e,ng_g,ng_r
    real(kind=dp_t), intent(in) :: rho(lo(1)-ng_r:,lo(2)-ng_r:,:)
    real(kind=dp_t), intent(inout) :: electro_mass_fluxx(lo(1)-ng_e:,lo(2)-ng_e:,:)
    real(kind=dp_t), intent(inout) :: electro_mass_fluxy(lo(1)-ng_e:,lo(2)-ng_e:,:)
    real(kind=dp_t), intent(in) ::         grad_Epotx(lo(1)-ng_g:,lo(2)-ng_g:)
    real(kind=dp_t), intent(in) ::         grad_Epoty(lo(1)-ng_g:,lo(2)-ng_g:)

    integer i,j,n

    
    do j = lo(2),hi(2)
        do i = lo(1),hi(1)+1
            do n=1,nspecies

                if(rho(i,j,n) .le. 0) then            

                    if(electro_mass_fluxx(i,j,n)>0) then
                        electro_mass_fluxx(i,j,1:nspecies) = 0
                    end if

		end if

                if(rho(i-1,j,n) .le. 0) then
                    if(electro_mass_fluxx(i,j,n)<0) then
                        electro_mass_fluxx(i,j,1:nspecies) = 0
                    end if
		end if
            
                end do
            end do
        end do       

    do j = lo(2),hi(2)+1
        do i = lo(1),hi(1)
            do n=1,nspecies


                if(rho(i,j,n) .le. 0) then            
                    if(electro_mass_fluxy(i,j,n)>0) then
                        electro_mass_fluxy(i,j,1:nspecies) = 0
                    end if
		end if

                if(rho(i,j-1,n) .le. 0) then
                    if(electro_mass_fluxy(i,j,n)<0) then
                        electro_mass_fluxy(i,j,1:nspecies) = 0
                    end if

                end if
            
                end do
            end do
        end do       
    


    end subroutine limit_emf_2d

  subroutine limit_emf_3d(rho, ng_r,electro_mass_fluxx,electro_mass_fluxy,electro_mass_fluxz,ng_e, &
                                     grad_Epotx,grad_Epoty,grad_Epotz,ng_g, &
                                     lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:),ng_e,ng_g,ng_r
    real(kind=dp_t), intent(in) :: rho(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:,:)
    real(kind=dp_t), intent(inout) :: electro_mass_fluxx(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
    real(kind=dp_t), intent(inout) :: electro_mass_fluxy(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
    real(kind=dp_t), intent(inout) :: electro_mass_fluxz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
    real(kind=dp_t), intent(in) ::         grad_Epotx(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:)
    real(kind=dp_t), intent(in) ::         grad_Epoty(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:)
    real(kind=dp_t), intent(in) ::         grad_Epotz(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:)

    integer i,j,k,n



    do k = lo(3),hi(3)
    do j = lo(2),hi(2)

        do i = lo(1),hi(1)

            do n=1,nspecies

                if(rho(i,j,k,n) .le. 0) then            

                    if(electro_mass_fluxx(i,j,k,n)>0) then

                        electro_mass_fluxx(i,j,k,1:nspecies) = 0

                    end if

                    if(electro_mass_fluxx(i+1,j,k,n)<0) then

                        electro_mass_fluxx(i+1,j,k,1:nspecies) = 0

                    end if

                    if(electro_mass_fluxy(i,j,k,n)>0) then

                        electro_mass_fluxy(i,j,k,1:nspecies) = 0

                    end if

                    if(electro_mass_fluxy(i,j+1,k,n)<0) then

                        electro_mass_fluxy(i,j+1,k,1:nspecies) = 0

                    end if

                    if(electro_mass_fluxz(i,j,k,n)>0) then

                        electro_mass_fluxz(i,j,k,1:nspecies) = 0

                    end if

                    if(electro_mass_fluxz(i,j,k+1,n)<0) then

                        electro_mass_fluxz(i,j,k+1,1:nspecies) = 0

                    end if

                end if
            
                end do
            end do
        end do       
        end do       

    end subroutine limit_emf_3d

end module electrodiffusive_mass_fluxdiv_module
