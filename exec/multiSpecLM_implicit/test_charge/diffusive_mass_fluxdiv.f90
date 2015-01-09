module diffusive_mass_fluxdiv_module

  use multifab_module
  use define_bc_module
  use bc_module
  use div_and_grad_module
  use probin_multispecies_module, only: nspecies, is_nonisothermal, mol_frac_bc_comp, &
                                        nspecies, correct_flux, use_charged_fluid, dielectric_const, &
                                        rho_part_bc_comp
  use probin_common_module, only: temp_bc_comp, barodiffusion_type
  use mass_flux_utilities_module
  use ml_layout_module
  use convert_stag_module
  use matvec_mul_module
  use matmat_mul_module
  use correction_flux_module
  use zero_edgeval_module
  use F95_LAPACK

  ! for charged fluid
  use fluid_charge_module
  use bndry_reg_module
  use ml_solve_module
  
  implicit none

  private

  public :: diffusive_mass_flux, diffusive_mass_fluxdiv

contains

  subroutine diffusive_mass_fluxdiv(mla,rho,rhotot,molarconc,rhoWchi,Gama,&
                                    diff_fluxdiv,Temp,zeta_by_Temp,gradp_baro, &
                                    flux_total,dx,the_bc_tower, &
                                    charge,grad_Epot)

    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: rho(:)
    type(multifab) , intent(in   )  :: rhotot(:)
    type(multifab) , intent(in   )  :: molarconc(:)
    type(multifab) , intent(in   )  :: rhoWchi(:)
    type(multifab) , intent(in   )  :: Gama(:)
    type(multifab) , intent(inout)  :: diff_fluxdiv(:)
    type(multifab) , intent(in   )  :: Temp(:)
    type(multifab) , intent(in   )  :: zeta_by_Temp(:)
    type(multifab) , intent(in   )  :: gradp_baro(:,:)
    type(multifab) , intent(inout)  :: flux_total(:,:)
    real(kind=dp_t), intent(in   )  :: dx(:,:)
    type(bc_tower) , intent(in   )  :: the_bc_tower
    type(multifab) , intent(inout)  :: charge(:)
    type(multifab) , intent(inout)  :: grad_Epot(:,:)

    ! local variables
    integer i,dm,n,nlevs

    ! local array of multifabs for grad and div; one for each direction
    type(multifab) :: flux(mla%nlevel,mla%dim)
    
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
    call diffusive_mass_flux(mla,rho,rhotot,molarconc,rhoWchi,Gama,Temp,&
                             zeta_by_Temp,gradp_baro,flux,dx,the_bc_tower, &
                             charge,grad_Epot)
    
    ! add fluxes to flux_total
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(flux_total(n,i),1,flux(n,i),1,nspecies,0)
       end do
    end do

    ! compute divergence of determinstic flux 
    call compute_div(mla,flux,diff_fluxdiv,dx,1,1,nspecies)
    
    ! destroy the multifab to free the memory
    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(flux(n,i))
       end do
    end do

  end subroutine diffusive_mass_fluxdiv
 
  subroutine diffusive_mass_flux(mla,rho,rhotot,molarconc,rhoWchi,Gama, &
                                 Temp,zeta_by_Temp,gradp_baro,flux,dx,the_bc_tower, &
                                 charge,grad_Epot)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: rho(:) 
    type(multifab) , intent(in   ) :: rhotot(:) 
    type(multifab) , intent(in   ) :: molarconc(:) 
    type(multifab) , intent(in   ) :: rhoWchi(:)  
    type(multifab) , intent(in   ) :: Gama(:)  
    type(multifab) , intent(in   ) :: Temp(:)  
    type(multifab) , intent(in   ) :: zeta_by_Temp(:)  
    type(multifab) , intent(in   ) :: gradp_baro(:,:)
    type(multifab) , intent(inout) :: flux(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(inout) :: charge(:)
    type(multifab) , intent(inout) :: grad_Epot(:,:)

    ! local variables
    integer :: n,i,s,dm,nlevs
 
    ! local face-centered multifabs 
    type(multifab)  :: rhoWchi_face(mla%nlevel,mla%dim)
    type(multifab)  :: Gama_face(mla%nlevel,mla%dim)
    type(multifab)  :: zeta_by_Temp_face(mla%nlevel,mla%dim)
    type(multifab)  :: flux_Temp(mla%nlevel,mla%dim)
    type(multifab)  :: baro_coef(mla%nlevel)
    type(multifab)  :: baro_coef_face(mla%nlevel,mla%dim)
  
    ! for electric potential Poisson solve
    type(multifab)  :: alpha(mla%nlevel)
    type(multifab)  :: beta(mla%nlevel,mla%dim)
    type(multifab)  :: Epot(mla%nlevel)
    type(multifab)  :: charge_coef(mla%nlevel)
    type(multifab)  :: charge_coef_face(mla%nlevel,mla%dim)
    type(bndry_reg) :: fine_flx(2:mla%nlevel)

    integer :: Epot_bc_comp

    ! NOTE THE BOUNDARY CONDITION COMPONENT ISN'T RIGHT
    Epot_bc_comp = scal_bc_comp
 
    dm    = mla%dim     ! dimensionality
    nlevs = mla%nlevel  ! number of levels 

    ! build local face-centered multifab with nspecies^2 component, zero ghost cells 
    ! and nodal in direction i
    do n=1,nlevs
       call multifab_build(baro_coef(n),mla%la(n),nspecies,rho(n)%ng)
       do i=1,dm
          call multifab_build_edge(rhoWchi_face(n,i),     mla%la(n),nspecies**2,0,i)
          call multifab_build_edge(Gama_face(n,i),        mla%la(n),nspecies**2,0,i)
          call multifab_build_edge(zeta_by_Temp_face(n,i),mla%la(n),nspecies,   0,i)
          call multifab_build_edge(baro_coef_face(n,i),mla%la(n),nspecies,   0,i)
          call multifab_build_edge(flux_Temp(n,i),        mla%la(n),1,          0,i)
       end do
    end do 

    ! compute face-centered rhoWchi from cell-centered values 
    call average_cc_to_face(nlevs, rhoWchi, rhoWchi_face, 1, tran_bc_comp, &
                            nspecies**2, the_bc_tower%bc_tower_array, .false.) 

    !==================================!
    ! compute flux-piece from molarconc
    !==================================! 

    ! calculate face-centrered grad(molarconc) 
    call compute_grad(mla, molarconc, flux, dx, 1, mol_frac_bc_comp, 1, nspecies, & 
                      the_bc_tower%bc_tower_array)

    ! compute face-centered Gama from cell-centered values 
    call average_cc_to_face(nlevs, Gama, Gama_face, 1, tran_bc_comp, &
                            nspecies**2, the_bc_tower%bc_tower_array, .false.)

    ! compute Gama*grad(molarconc): Gama is nspecies^2 matrix; grad(x) is nspecies component vector 
    do n=1,nlevs
       do i=1,dm
          call matvec_mul(mla, flux(n,i), Gama_face(n,i), nspecies)
       end do
    end do    

    if(is_nonisothermal) then
    
       !====================================!
       ! compute flux-piece from Temperature 
       !====================================! 
 
       ! calculate face-centrered grad(T) 
       call compute_grad(mla, Temp, flux_Temp, dx, 1, temp_bc_comp, 1, 1, the_bc_tower%bc_tower_array)
    
       ! compute face-centered zeta_by_T from cell-centered values 
       call average_cc_to_face(nlevs, zeta_by_Temp, zeta_by_Temp_face, 1, tran_bc_comp, &
                               nspecies, the_bc_tower%bc_tower_array, .false.) 

       ! compute zeta_by_T*grad(T): zeta_by_T is nspecies component vector; grad(T) is scalar
       do n=1,nlevs
          do i=1,dm
             do s=1,nspecies
                call multifab_mult_mult_c(zeta_by_Temp_face(n,i), s, flux_Temp(n,i), 1, 1)
             end do
          end do
       end do  
    
       !===============================!
       ! assemble different flux-pieces 
       !===============================! 
       do n=1,nlevs
          do i=1,dm
             call multifab_plus_plus(flux(n,i), zeta_by_Temp_face(n,i), 0)
          end do
       end do  
   
    end if

    if (barodiffusion_type .gt. 0) then
    
       !====================================!
       ! compute flux-piece from barodiffusion
       !====================================! 

       ! compute cell-centered barodiffusion coefficient, (phi-w) / (n kB T)
       call compute_baro_coef(mla,baro_coef,rho,rhotot,Temp)

       ! average baro_coef to faces
       call average_cc_to_face(nlevs, baro_coef, baro_coef_face, 1, scal_bc_comp, &
                               nspecies, the_bc_tower%bc_tower_array, .false.)

       ! store the fluxes, baro_coef(1:nspecies) * gradp_baro, in baro_coef_face
       do n=1,nlevs
          do i=1,dm
             do s=1,nspecies
                call multifab_mult_mult_c(baro_coef_face(n,i), s, gradp_baro(n,i), 1, 1)
             end do
          end do
       end do
    
       !===============================!
       ! assemble different flux-pieces 
       !===============================! 
       do n=1,nlevs
          do i=1,dm
             call multifab_plus_plus(flux(n,i), baro_coef_face(n,i), 0)
          end do
       end do  

    end if

    if (use_charged_fluid) then

       ! compute total charge
       call compute_total_charge(mla,rho,charge)

       ! solve poisson equation for phi (the electric potential)
       ! -del dot epsilon grad Phi = charge
       do n=1,nlevs

          call multifab_build(Epot(n),mla%la(n),1,1)

          ! set alpha=0
          call multifab_build(alpha(n),mla%la(n),1,0)
          call setval(alpha(n),0.d0,all=.true.)

          ! set beta=dielectric_const
          do i=1,dm
             call multifab_build_edge(beta(n,i),mla%la(n),1,0,i)
             call setval(beta(n,i),dielectric_const,all=.true.)
          end do

       end do

       ! build the boundary flux register
       do n=2,nlevs
          call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
       end do

       ! solve (alpha - del dot beta grad) Epot = charge
       call ml_cc_solve(mla,charge,Epot,fine_flx,alpha,beta,dx,the_bc_tower,Epot_bc_comp)
          
       do n=1,nlevs
          call multifab_destroy(alpha(n))
          do i=1,dm
             call multifab_destroy(beta(n,i))
          end do
       end do
          
       ! compute the gradient of the electric potential
       call compute_grad(mla,Epot,grad_Epot,dx,1,Epot_bc_comp,1,1,the_bc_tower%bc_tower_array)
       do n=1,nlevs
          call multifab_destroy(Epot(n))
       end do

       ! compute the charge flux coefficient
       do n=1,nlevs
          call multifab_build(charge_coef(n),mla%la(n),nspecies,1)
       end do
       call compute_charge_coef(mla,rho,rhotot,Temp,charge,charge_coef)

       ! average charge flux coefficient to faces
       do n=1,nlevs
          do i=1,dm
             call multifab_build_edge(charge_coef_face(n,i),mla%la(n),nspecies,0,i)
          end do
       end do
       call average_cc_to_face(nlevs,charge_coef,charge_coef_face,1,rho_part_bc_comp,nspecies, &
                               the_bc_tower%bc_tower_array,.true.)

       ! multiply flux coefficient by gradient of electric potential
       do n=1,nlevs
          do i=1,dm
             do s=1,nspecies
                call multifab_mult_mult_c(charge_coef_face(n,i), s, grad_Epot(n,i), 1, 1)
             end do
          end do
       end do

       ! add charge flux to the running total
       do n=1,nlevs
          do i=1,dm
             call multifab_plus_plus(flux(n,i), charge_coef_face(n,i), 0)
          end do
       end do  

       do n=1,nlevs
          call multifab_destroy(charge_coef(n))
          do i=1,dm
             call multifab_destroy(charge_coef_face(n,i))
          end do
       end do       

    end if

    ! compute rhoWchi * totalflux (on faces) 
    do n=1,nlevs
       do i=1,dm
          call matvec_mul(mla, flux(n,i), rhoWchi_face(n,i), nspecies)
       end do
    end do    

    ! If there are walls with zero-flux boundary conditions
    if(is_nonisothermal) then
       do n=1,nlevs
          call zero_edgeval_walls(flux(n,:),1,nspecies,the_bc_tower%bc_tower_array(n))
       end do   
    end if

    !correct fluxes to ensure mass conservation to roundoff
    if (correct_flux .and. (nspecies .gt. 1)) then
       !write(*,*) "Checking conservation of deterministic fluxes"
       call correction_flux(mla, rho, rhotot, flux, the_bc_tower%bc_tower_array)
    end if
    
    ! destroy B^(-1)*Gama multifab to prevent leakage in memory
    do n=1,nlevs
       call multifab_destroy(baro_coef(n))
       do i=1,dm
          call multifab_destroy(rhoWchi_face(n,i))
          call multifab_destroy(Gama_face(n,i))
          call multifab_destroy(zeta_by_Temp_face(n,i))
          call multifab_destroy(flux_Temp(n,i))
          call multifab_destroy(baro_coef_face(n,i))
       end do
    end do

  end subroutine diffusive_mass_flux

end module diffusive_mass_fluxdiv_module
