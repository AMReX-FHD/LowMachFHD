module Epot_mass_fluxdiv_module

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
  use probin_multispecies_module, only: nspecies
  use probin_gmres_module, only: mg_verbose
  use probin_charged_module, only: dielectric_const, Epot_wall_bc_type
  
  implicit none

  private

  public :: Epot_mass_flux, Epot_mass_fluxdiv, inhomogeneous_neumann_fix

contains

  subroutine Epot_mass_fluxdiv(mla,rho,Epot_fluxdiv,Temp,rhoWchi, &
                               flux_total,dx,the_bc_tower,charge,grad_Epot,Epot)

    ! this computes "Epot_fluxdiv = -div(F) = div(A_Phi grad Phi)"
    !               "grad_Epot    = grad Phi"

    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: rho(:)
    type(multifab) , intent(inout)  :: Epot_fluxdiv(:)
    type(multifab) , intent(in   )  :: Temp(:)
    type(multifab) , intent(in   )  :: rhoWchi(:)
    type(multifab) , intent(inout)  :: flux_total(:,:)
    real(kind=dp_t), intent(in   )  :: dx(:,:)
    type(bc_tower) , intent(in   )  :: the_bc_tower
    type(multifab) , intent(inout)  :: charge(:)
    type(multifab) , intent(inout)  :: grad_Epot(:,:)
    type(multifab) , intent(inout)  :: Epot(:)

    ! local variables
    integer i,dm,n,nlevs

    ! local array of multifabs for grad and div; one for each direction
    type(multifab) :: flux(mla%nlevel,mla%dim)
    
    type(bl_prof_timer), save :: bpt

    call build(bpt, "Epot_mass_fluxdiv")
    
    nlevs = mla%nlevel  ! number of levels 
    dm    = mla%dim     ! dimensionality
 
    ! build the local multifabs
    do n=1,nlevs
       do i=1,dm
          ! flux(i) is face-centered, has nspecies component, zero ghost 
          ! cells & nodal in direction i
          call multifab_build_edge(flux(n,i),mla%la(n),nspecies,0,i)
       end do
    end do   
    
    ! compute the face-centered flux (each direction: cells+1 faces while 
    ! cells contain interior+2 ghost cells) 
    call Epot_mass_flux(mla,rho,Temp,rhoWchi,flux,dx,the_bc_tower,charge,grad_Epot,Epot)
    
    ! add fluxes to flux_total
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(flux_total(n,i),1,flux(n,i),1,nspecies,0)
       end do
    end do

    ! compute divergence of determinstic flux 
    call compute_div(mla,flux,Epot_fluxdiv,dx,1,1,nspecies)
    
    ! destroy the multifab to free the memory
    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(flux(n,i))
       end do
    end do

    call destroy(bpt)

  end subroutine Epot_mass_fluxdiv
 
  subroutine Epot_mass_flux(mla,rho,Temp,rhoWchi,flux,dx, &
                            the_bc_tower,charge,grad_Epot,Epot)

    ! this computes "-F = A_Phi grad Phi"

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: rho(:) 
    type(multifab) , intent(in   ) :: Temp(:)  
    type(multifab) , intent(in   ) :: rhoWchi(:)
    type(multifab) , intent(inout) :: flux(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(inout) :: charge(:)
    type(multifab) , intent(inout) :: grad_Epot(:,:)
    type(multifab) , intent(inout) :: Epot(:)

    ! local variables
    integer :: n,i,s,dm,nlevs
 
    ! local face-centered multifabs 
    type(multifab)  :: rhoWchi_face(mla%nlevel,mla%dim)

    ! for electric potential Poisson solve
    type(multifab)  :: alpha(mla%nlevel)
    type(multifab)  :: beta(mla%nlevel,mla%dim)
    type(multifab)  :: charge_coef(mla%nlevel)
    type(multifab)  :: charge_coef_face(mla%nlevel,mla%dim)
    type(bndry_reg) :: fine_flx(mla%nlevel)
  
    real(kind=dp_t) :: avg_charge

    type(bl_prof_timer), save :: bpt
    
    call build(bpt,"Epot_mass_flux")

    dm    = mla%dim     ! dimensionality
    nlevs = mla%nlevel  ! number of levels 
    
    do n=1,nlevs
       call multifab_build(alpha(n),mla%la(n),1,0)
       call multifab_build(charge_coef(n),mla%la(n),nspecies,1)
       do i=1,dm
          call multifab_build_edge(beta(n,i),mla%la(n),1,0,i)
          call multifab_build_edge(charge_coef_face(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge(rhoWchi_face(n,i),mla%la(n),nspecies**2,0,i)
       end do
    end do

    ! compute face-centered rhoWchi from cell-centered values 
    call average_cc_to_face(nlevs, rhoWchi, rhoWchi_face, 1, tran_bc_comp, &
                            nspecies**2, the_bc_tower%bc_tower_array, .false.) 

    ! solve poisson equation for phi (the electric potential)
    ! -del dot epsilon grad Phi = charge
    do n=1,nlevs

       call setval(Epot(n),0.d0,all=.true.)

       ! fill ghost cells for Epot at walls using Dirichlet value
       call multifab_physbc(Epot(n),1,Epot_bc_comp,1,the_bc_tower%bc_tower_array(n), &
                            dx_in=dx(n,:))

       ! set alpha=0
       call setval(alpha(n),0.d0,all=.true.)

       ! set beta=dielectric_const
       do i=1,dm
          call setval(beta(n,i),dielectric_const,all=.true.)
       end do

    end do

    ! build the boundary flux register
    do n=1,nlevs
       call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
    end do

    if (all(mla%pmask(1:dm))) then
       avg_charge = multifab_sum_c(charge(1),1,1) / multifab_volume(charge(1))
       call multifab_sub_sub_s_c(charge(1),1,avg_charge,1,0)
       if (abs(avg_charge) .gt. 1.d-12) then
          if (parallel_IOProcessor()) then
             print*,'average charge =',avg_charge
          end if
          call bl_warn("Warning: average charge is not zero")
       end if
    end if

    ! for inhomogeneous Neumann bc's for electric potential, put in homogeneous form
    if (Epot_wall_bc_type .eq. 2) then
       call inhomogeneous_neumann_fix(mla,charge,dx,the_bc_tower)
    end if

    ! solve (alpha - del dot beta grad) Epot = charge
    call ml_cc_solve(mla,charge,Epot,fine_flx,alpha,beta,dx,the_bc_tower,Epot_bc_comp, &
                     verbose=mg_verbose)

    do n=1,nlevs
       call bndry_reg_destroy(fine_flx(n))
    end do

    ! compute the gradient of the electric potential
    call compute_grad(mla,Epot,grad_Epot,dx,1,Epot_bc_comp,1,1,the_bc_tower%bc_tower_array)

    ! hack - add potential in periodic problem
!    do n=1,nlevs
!       call multifab_plus_plus_s(grad_Epot(n,1),-1.d10,0)
!    end do

    ! compute the charge flux coefficient
    call compute_charge_coef(mla,rho,Temp,charge_coef)

    ! average charge flux coefficient to faces, store in flux
    call average_cc_to_face(nlevs,charge_coef,flux,1,c_bc_comp,nspecies, &
                            the_bc_tower%bc_tower_array,.true.)

    ! multiply flux coefficient by gradient of electric potential
    do n=1,nlevs
       do i=1,dm
          do s=1,nspecies
             call multifab_mult_mult_c(flux(n,i), s, grad_Epot(n,i), 1, 1)
          end do
       end do
    end do

    ! compute -rhoWchi * (... ) on faces
    do n=1,nlevs
       do i=1,dm
          call matvec_mul(mla, flux(n,i), rhoWchi_face(n,i), nspecies)
          call multifab_mult_mult_s(flux(n,i),-1.d0)
       end do
    end do    

    ! zero the total mass flux on walls to make sure 
    ! that the potential gradient matches the species gradient
    do n=1,nlevs
       call zero_edgeval_walls(flux(n,:),1,nspecies, &
                               the_bc_tower%bc_tower_array(n))
    end do

    do n=1,nlevs
       call multifab_destroy(alpha(n))
       call multifab_destroy(charge_coef(n))
       do i=1,dm
          call multifab_destroy(rhoWchi_face(n,i))
          call multifab_destroy(beta(n,i))
          call multifab_destroy(charge_coef_face(n,i))
       end do
    end do

    call destroy(bpt)

  end subroutine Epot_mass_flux
  
  ! We would like to solve A x = b with inhomogeneous bc's
  ! Here, "A" is -dielectric_const * Lap
  ! This is equivalent to A_H x = b - A x_H, where
  !   A   is the inhomogeneous operator
  !   A_H is the homogeneous operator
  !   x_H is a multifab filled with zeros, but ghost cells filled to respect bc's
  ! We use this for walls with inhomogeneous Neumann conditions on the electric potential
  subroutine inhomogeneous_neumann_fix(mla,charge,dx,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: charge(:)
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

       ! multiply zerofab everywhere (including ghost cells) by dielectric_const since
       ! we are incrementing the RHS (charge) by -A x_H
       call multifab_mult_mult_s(zerofab(n), dielectric_const, 1)

    end do

    ! compute gradient of zerofab
    call compute_grad(mla,zerofab,gradphi,dx,1,Epot_bc_comp,1,1,the_bc_tower%bc_tower_array)

    ! increment charge with negative divergence
    call compute_div(mla,gradphi,charge,dx,1,1,1,increment_in=.true.)

    do n=1,nlevs
       call multifab_fill_boundary(charge(n))
    end do

    do n=1,nlevs
       call multifab_destroy(zerofab(n))
       do i=1,dm
          call multifab_destroy(gradphi(n,i))
       end do
    end do

  end subroutine inhomogeneous_neumann_fix

end module Epot_mass_fluxdiv_module
