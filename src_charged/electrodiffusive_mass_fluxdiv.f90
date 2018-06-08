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
  use probin_common_module, only: nspecies, variance_coef_mass, shift_cc_to_boundary
  use probin_charged_module, only: Epot_wall_bc_type, Epot_wall, E_ext_type, electroneutral, &
                                   zero_eps_on_wall_type, epot_mg_verbose, epot_mg_abs_tol, &
                                   epot_mg_rel_tol, charge_per_mass
  
  use fabio_module

  implicit none

  private

  public :: electrodiffusive_mass_flux, electrodiffusive_mass_fluxdiv, &
       inhomogeneous_neumann_fix

contains

  subroutine electrodiffusive_mass_fluxdiv(mla,rho,Temp,rhoWchi, &
                                           diff_mass_flux,diff_mass_fluxdiv, &
                                           stoch_mass_fluxdiv, &
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
    type(multifab) , intent(in   )  :: stoch_mass_fluxdiv(:)
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
                                    diff_mass_fluxdiv,stoch_mass_fluxdiv,dx,the_bc_tower, &
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
                                        diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                        dx,the_bc_tower,charge,grad_Epot,Epot, &
                                        permittivity,dt,zero_initial_Epot)

    ! this computes "-F = A_Phi grad Phi"

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: rho(:) 
    type(multifab) , intent(in   ) :: Temp(:)  
    type(multifab) , intent(in   ) :: rhoWchi(:)
    type(multifab) , intent(inout) :: electro_mass_flux(:,:)
    type(multifab) , intent(in   ) :: diff_mass_fluxdiv(:)
    type(multifab) , intent(in   ) :: stoch_mass_fluxdiv(:)
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
    type(multifab)  :: charge_coef_face(mla%nlevel,mla%dim)
    type(multifab)  ::  permittivity_fc(mla%nlevel,mla%dim)
    type(multifab)  ::            E_ext(mla%nlevel,mla%dim)
    type(multifab) ::             A_Phi(mla%nlevel,mla%dim)

    type(bndry_reg) :: fine_flx(mla%nlevel)
  
    real(kind=dp_t) :: sum

    real(kind=dp_t) :: Epot_wall_save(2,mla%dim)

    real(kind=dp_t) :: z_temp(nspecies), norm, epot_mg_abs_tol_temp

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
          call multifab_build_edge(charge_coef_face(n,i),mla%la(n),   nspecies,0,i)
          call multifab_build_edge(    rhoWchi_face(n,i),mla%la(n),nspecies**2,0,i)
          call multifab_build_edge( permittivity_fc(n,i),mla%la(n),          1,0,i)
          call multifab_build_edge(           E_Ext(n,i),mla%la(n),          1,0,i)
          call multifab_build_edge(           A_Phi(n,i),mla%la(n),   nspecies,0,i)
       end do
    end do

    ! build the boundary flux register
    do n=1,nlevs
       call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
    end do

    ! if periodic, ensure charge sums to zero by subtracting off the averaeg
    if (all(mla%pmask(1:dm))) then
       sum = multifab_sum_c(charge(1),1,1) / multifab_volume(charge(1))
       call multifab_sub_sub_s_c(charge(1),1,sum,1,0)
       if (abs(sum) .gt. 1.d-12) then
          if (parallel_IOProcessor()) then
             print*,'average charge =',sum
          end if
          call bl_warn("Warning: average charge is not zero")
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

       ! compute A_Phi for Poisson solve (does not have z^T)
       call implicit_potential_coef(mla,rho,Temp,A_Phi,the_bc_tower)
              
       ! compute z^T A_Phi^n, store in solver_beta
       call dot_with_z_face(mla,A_Phi,beta)

       ! compute RHS = div (z^T (F_d + F_s))

       ! first set rhsvec = F_d + F_s
       do n=1,nlevs
          call multifab_copy_c(rhsvec(n),1,diff_mass_fluxdiv(n),1,nspecies,0)
          if (variance_coef_mass .ne. 0.d0) then
             call multifab_plus_plus_c(rhsvec(n),1,stoch_mass_fluxdiv(n),1,nspecies,0)
          end if
       end do

       ! save original charge_per_mass
       z_temp(1:nspecies) = charge_per_mass(1:nspecies)

       ! take absolute value of charge per mass
       charge_per_mass(1:nspecies) = abs(charge_per_mass(1:nspecies))

       ! dot abs(z) with F
       call dot_with_z(mla,rhsvec,rhs)
       ! compute norm
       norm = multifab_norm_inf(rhsvec(1))

       ! set absolute tolerance to be the norm*epot_mg_rel_tol
       epot_mg_abs_tol_temp = epot_mg_abs_tol
       epot_mg_abs_tol = norm*epot_mg_rel_tol

       ! restore charge per mass
       charge_per_mass(1:nspecies) = z_temp(1:nspecies)

       ! compute rhs for Poisson zolve, z dot F
       call dot_with_z(mla,rhsvec,rhs)

    else

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
       if (any(Epot_wall_bc_type(1:2,1:dm).eq. 2)) then

          ! save the numerical values for the Dirichlet and Neumann conditions
          Epot_wall_save(1:2,1:dm) = Epot_wall(1:2,1:dm)

          ! for Dirichlet conditions, temporarily set the numerical values to zero
          ! so we can put the Neumann boundaries into homogeneous form
          do comp=1,dm
             do i=1,dm
                if (Epot_wall_bc_type(comp,1) .eq. 1) then
                   Epot_wall(comp,1) = 0.d0
                end if
             end do
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

       ! compute div (epsilon*E_ext) and add it to solver rhs
       call compute_div(mla,permittivity_fc,rhs,dx,1,1,1,increment_in=.true.)

    end if

    ! solve (alpha - del dot beta grad) Epot = charge
    call ml_cc_solve(mla,rhs,Epot,fine_flx,alpha,beta,dx(:,1:dm),the_bc_tower,Epot_bc_comp, &
                     eps=epot_mg_rel_tol, &
                     abs_eps=epot_mg_abs_tol, &
                     verbose=epot_mg_verbose, &
                     ok_to_fix_singular=.false.)

    if (electroneutral) then
       epot_mg_abs_tol = epot_mg_abs_tol_temp
    end if
    

    ! for periodic problems subtract off the average of Epot
    ! we can generalize this later for walls
    if ( all(mla%pmask(1:dm)) ) then
       sum = multifab_sum(Epot(1)) / boxarray_dvolume(get_boxarray(Epot(1)))
       call multifab_sub_sub_s(Epot(1),sum)
       do n=1,nlevs
          call multifab_physbc(Epot(n),1,Epot_bc_comp,1,the_bc_tower%bc_tower_array(n), &
                               dx_in=dx(n,:))
          call multifab_fill_boundary(Epot(n))
       end do
    else
       ! we need to fill the ghost cells so the inhomogeneous Neumann phi case
       ! has properly filled ghost cells (since the solver assumed homogeneous BC's)
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
       do n=1,nlevs
          do i=1,dm
             call multifab_plus_plus_c(grad_Epot(n,i),1,E_ext(n,i),1,1,0)
          end do
       end do
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

    ! compute -rhoWchi * (... ) on faces
    do n=1,nlevs
       do i=1,dm
          call matvec_mul(mla, electro_mass_flux(n,i), rhoWchi_face(n,i), nspecies)
       end do
    end do    

    ! zero the total mass flux on walls to make sure 
    ! that the potential gradient matches the species gradient
    do n=1,nlevs
       call zero_edgeval_walls(electro_mass_flux(n,:),1,nspecies, &
                               the_bc_tower%bc_tower_array(n))
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
          call multifab_destroy(charge_coef_face(n,i))
          call multifab_destroy(permittivity_fc(n,i))
          call multifab_destroy(E_ext(n,i))
          call multifab_destroy(A_Phi(n,i))
       end do
    end do

    call destroy(bpt)

  end subroutine electrodiffusive_mass_flux
  
  ! We would like to solve A x = b with inhomogeneous bc's
  ! Here, "A" is -epsilon * Lap
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

end module electrodiffusive_mass_fluxdiv_module
