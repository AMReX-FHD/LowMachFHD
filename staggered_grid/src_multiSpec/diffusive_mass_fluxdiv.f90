module diffusive_mass_fluxdiv_module

  use multifab_module
  use define_bc_module
  use bc_module
  use debug_module
  use div_and_grad_module
  use probin_multispecies_module, only: is_nonisothermal, correct_flux, use_multiphase
  use probin_common_module, only: barodiffusion_type, nspecies, shift_cc_to_boundary
  use mass_flux_utilities_module
  use ml_layout_module
  use convert_stag_module
  use matvec_mul_module
  use matmat_mul_module
  use correction_flux_module
  use zero_edgeval_module
  
  implicit none

  private

  public :: diffusive_mass_flux, diffusive_mass_fluxdiv

contains

  subroutine diffusive_mass_fluxdiv(mla,rho,rhotot,molarconc,rhoWchi,Gama,&
                                    diff_mass_fluxdiv,Temp,zeta_by_Temp,gradp_baro, &
                                    diff_mass_flux,dx,the_bc_tower)

    ! this computes divergence of "-F = rho*W*chi*Gamma*grad(x) - ..."

    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: rho(:)
    type(multifab) , intent(in   )  :: rhotot(:)
    type(multifab) , intent(in   )  :: molarconc(:)
    type(multifab) , intent(in   )  :: rhoWchi(:)
    type(multifab) , intent(in   )  :: Gama(:)
    type(multifab) , intent(inout)  :: diff_mass_fluxdiv(:)
    type(multifab) , intent(in   )  :: Temp(:)
    type(multifab) , intent(in   )  :: zeta_by_Temp(:)
    type(multifab) , intent(in   )  :: gradp_baro(:,:)
    type(multifab) , intent(inout)  :: diff_mass_flux(:,:)
    real(kind=dp_t), intent(in   )  :: dx(:,:)
    type(bc_tower) , intent(in   )  :: the_bc_tower

    ! local variables
    type(bl_prof_timer), save :: bpt

    call build(bpt, "diffusive_mass_fluxdiv")

    ! compute the face-centered flux (each direction: cells+1 faces while 
    ! cells contain interior+2 ghost cells) 
    call diffusive_mass_flux(mla,rho,rhotot,molarconc,rhoWchi,Gama,Temp,&
                             zeta_by_Temp,gradp_baro,diff_mass_flux,dx,the_bc_tower)
    
    ! compute divergence of determinstic flux 
    call compute_div(mla,diff_mass_flux,diff_mass_fluxdiv,dx,1,1,nspecies)

    call destroy(bpt)

  end subroutine diffusive_mass_fluxdiv
 
  subroutine diffusive_mass_flux(mla,rho,rhotot,molarconc,rhoWchi,Gama, &
                                 Temp,zeta_by_Temp,gradp_baro,diff_mass_flux,dx, &
                                 the_bc_tower)

    ! this computes "-F = rho*W*chi*Gamma*grad(x) - ..."

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: rho(:) 
    type(multifab) , intent(in   ) :: rhotot(:) 
    type(multifab) , intent(in   ) :: molarconc(:) 
    type(multifab) , intent(in   ) :: rhoWchi(:)  
    type(multifab) , intent(in   ) :: Gama(:)  
    type(multifab) , intent(in   ) :: Temp(:)  
    type(multifab) , intent(in   ) :: zeta_by_Temp(:)  
    type(multifab) , intent(in   ) :: gradp_baro(:,:)
    type(multifab) , intent(inout) :: diff_mass_flux(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: n,i,s,dm,nlevs
 
    ! local face-centered multifabs 
    type(multifab)  :: rhoWchi_face(mla%nlevel,mla%dim)
    type(multifab)  :: Gama_face(mla%nlevel,mla%dim)
    type(multifab)  :: zeta_by_Temp_face(mla%nlevel,mla%dim)
    type(multifab)  :: thermodiff_mass_flux(mla%nlevel,mla%dim)
    type(multifab)  :: baro_coef(mla%nlevel)
    type(multifab)  :: baro_coef_face(mla%nlevel,mla%dim)
  
    type(bl_prof_timer), save :: bpt
    
    call build(bpt,"diffusive_mass_flux")

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
          call multifab_build_edge(thermodiff_mass_flux(n,i),        mla%la(n),1,          0,i)
       end do
    end do 

    ! compute face-centered rhoWchi from cell-centered values 
    if (any(shift_cc_to_boundary(:,:) .eq. 1)) then
       call shift_cc_to_boundary_face(nlevs, rhoWchi, rhoWchi_face, 1, tran_bc_comp, &
                                      nspecies**2, the_bc_tower%bc_tower_array, .false.) 
    else
       call average_cc_to_face(nlevs, rhoWchi, rhoWchi_face, 1, tran_bc_comp, &
                               nspecies**2, the_bc_tower%bc_tower_array, .false.) 
    end if

!    print*, 'rhoWchi'
!    call print_edge(rhoWchi_face,2,8,1,1)
!    call print_edge(rhoWchi_face,2,8,1,2)
!    call print_edge(rhoWchi_face,2,8,1,3)
!    call print_edge(rhoWchi_face,2,8,1,4)
  !!  print*, 'rhoWchi'
  !!  print*, rhoWchi

    !==================================!
    ! compute flux-piece from molarconc
    !==================================! 

    ! calculate face-centrered grad(molarconc) 
    call compute_grad(mla, molarconc, diff_mass_flux, dx, 1, mol_frac_bc_comp, 1, nspecies, & 
                      the_bc_tower%bc_tower_array)

!    print*, 'diff mass flux at grad x'
!    call print_edge(diff_mass_flux,2,8,1,1)

    ! compute face-centered Gama from cell-centered values 
    if (any(shift_cc_to_boundary(:,:) .eq. 1)) then
       call shift_cc_to_boundary_face(nlevs, Gama, Gama_face, 1, tran_bc_comp, &
                                      nspecies**2, the_bc_tower%bc_tower_array, .false.)
    else
       call average_cc_to_face(nlevs, Gama, Gama_face, 1, tran_bc_comp, &
                               nspecies**2, the_bc_tower%bc_tower_array, .false.)
    end if

    !KK I think we need to insert Gama limiting here

    call limit_Gama(molarconc,Gama_face)

!    print*, 'Gama face'
!    call print_edge(Gama_face,2,8,1,1)
!    call print_edge(Gama_face,2,8,1,2)
!    call print_edge(Gama_face,2,8,1,3)
!    call print_edge(Gama_face,2,8,1,4)

    ! compute Gama*grad(molarconc): Gama is nspecies^2 matrix; grad(x) is nspecies component vector 
    do n=1,nlevs
       do i=1,dm
          call matvec_mul(mla, diff_mass_flux(n,i), Gama_face(n,i), nspecies)
       end do
    end do   
 
!    print*, 'diff mass flux at GAMA grad x'
!    call print_edge(diff_mass_flux,2,8,1,1)

    if(use_multiphase)then

        call compute_higher_order_term(mla, molarconc, diff_mass_flux, dx, 1, mol_frac_bc_comp, 1, nspecies, & 
                      the_bc_tower%bc_tower_array)

    endif
     

    if(is_nonisothermal) then
    
       !====================================!
       ! compute flux-piece from Temperature 
       !====================================! 
 
       ! calculate face-centrered grad(T) 
       call compute_grad(mla, Temp, thermodiff_mass_flux, dx, 1, temp_bc_comp, 1, 1, the_bc_tower%bc_tower_array)
    
       ! compute face-centered zeta_by_T from cell-centered values 
       if (any(shift_cc_to_boundary(:,:) .eq. 1)) then
          call shift_cc_to_boundary_face(nlevs, zeta_by_Temp, zeta_by_Temp_face, 1, tran_bc_comp, &
                                         nspecies, the_bc_tower%bc_tower_array, .false.) 
       else
          call average_cc_to_face(nlevs, zeta_by_Temp, zeta_by_Temp_face, 1, tran_bc_comp, &
                                  nspecies, the_bc_tower%bc_tower_array, .false.) 
       end if

       ! compute zeta_by_T*grad(T): zeta_by_T is nspecies component vector; grad(T) is scalar
       do n=1,nlevs
          do i=1,dm
             do s=1,nspecies
                call multifab_mult_mult_c(zeta_by_Temp_face(n,i), s, thermodiff_mass_flux(n,i), 1, 1)
             end do
          end do
       end do  
    
       !===============================!
       ! assemble different flux-pieces 
       !===============================! 
       do n=1,nlevs
          do i=1,dm
             call multifab_plus_plus(diff_mass_flux(n,i), zeta_by_Temp_face(n,i), 0)
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
       if (any(shift_cc_to_boundary(:,:) .eq. 1)) then
          call shift_cc_to_boundary_face(nlevs, baro_coef, baro_coef_face, 1, scal_bc_comp, &
                                         nspecies, the_bc_tower%bc_tower_array, .false.)
       else
          call average_cc_to_face(nlevs, baro_coef, baro_coef_face, 1, scal_bc_comp, &
                                  nspecies, the_bc_tower%bc_tower_array, .false.)
       end if

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
             call multifab_plus_plus(diff_mass_flux(n,i), baro_coef_face(n,i), 0)
          end do
       end do  

    end if

    ! compute -rhoWchi * (Gamma*grad(x) + ... ) on faces
    do n=1,nlevs
       do i=1,dm
          call matvec_mul(mla, diff_mass_flux(n,i), rhoWchi_face(n,i), nspecies)
       end do
    end do    
 
!    print*, 'diff mass flux at rhowchi GAMA grad x'
!!    call print_edge(diff_mass_flux,2,8,1,1)

    ! If there are walls with zero-flux boundary conditions
    if(is_nonisothermal) then
       do n=1,nlevs
          call zero_edgeval_walls(diff_mass_flux(n,:),1,nspecies,the_bc_tower%bc_tower_array(n))
       end do   
    end if

    !correct fluxes to ensure mass conservation to roundoff
    if (correct_flux .and. (nspecies .gt. 1)) then
       !write(*,*) "Checking conservation of deterministic fluxes"
       call correction_flux(mla, rho, rhotot, diff_mass_flux, the_bc_tower%bc_tower_array)
    end if
    
    ! destroy B^(-1)*Gama multifab to prevent leakage in memory
    do n=1,nlevs
       call multifab_destroy(baro_coef(n))
       do i=1,dm
          call multifab_destroy(rhoWchi_face(n,i))
          call multifab_destroy(Gama_face(n,i))
          call multifab_destroy(zeta_by_Temp_face(n,i))
          call multifab_destroy(thermodiff_mass_flux(n,i))
          call multifab_destroy(baro_coef_face(n,i))
       end do
    end do

    call destroy(bpt)

  end subroutine diffusive_mass_flux

!make Gama identity for small concentrations
  subroutine limit_Gama(molarconc, Gama)

    type(multifab) , intent(inout) :: Gama(:,:)
    type(multifab) , intent(in) :: molarconc(:)

    ! Local
    integer                  :: lo(get_dim(molarconc(1))),hi(get_dim(molarconc(1)))
    integer                  :: ng_m,ng_G,i,dm
    real(kind=dp_t), pointer :: mp(:,:,:,:)
    real(kind=dp_t), pointer :: Gpx(:,:,:,:), Gpy(:,:,:,:), Gpz(:,:,:,:)

    dm = get_dim(molarconc(1))
    ng_m = nghost(molarconc(1))
    ng_G = nghost(Gama(1,1))
    
    do i=1,nfabs(molarconc(1))
        mp => dataptr(molarconc(1),i)
       Gpx => dataptr(Gama(1,1),i)
       Gpy => dataptr(Gama(1,2),i)

       lo = lwb(get_box(molarconc(1),i))
       hi = upb(get_box(molarconc(1),i))
          select case (dm)
          case (2)
             call limit_Gama_2d(mp(:,:,1,:), ng_m, &
                                          Gpx(:,:,1,:), Gpy(:,:,1,:), ng_G, &
                                          lo, hi)
          case (3)
             Gpz => dataptr(Gama(1,3),i)
             call limit_Gama_3d(mp(:,:,:,:), ng_m, &
                                          Gpx(:,:,:,:), Gpy(:,:,:,:), Gpz(:,:,:,:), ng_G, &
                                          lo, hi)
          end select
    end do
 
  end subroutine limit_Gama

  subroutine limit_Gama_2d(molarconc, ng_m,Gpx,Gpy,ng_G, &
                                     lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:),ng_m,ng_G
    real(kind=dp_t), intent(in) :: molarconc(lo(1)-ng_m:,lo(2)-ng_m:,:)
    real(kind=dp_t), intent(inout) :: Gpx(lo(1)-ng_G:,lo(2)-ng_G:,:)
    real(kind=dp_t), intent(inout) :: Gpy(lo(1)-ng_G:,lo(2)-ng_G:,:)

    integer i,j,n
    real(kind=dp_t) eepsilon

    eepsilon = 0.0001

    do j = lo(2),hi(2)
        do i = lo(1),hi(1)+1
            do n=1,nspecies

                if(min(molarconc(i,j,n),molarconc(i-1,j,n)) .le. eepsilon) then !need to make this an input parameter?            

                    Gpx(i,j,1)=1.d0
                    Gpx(i,j,2)=0.d0
                    Gpx(i,j,3)=0.d0
                    Gpx(i,j,4)=1.d0

                endif

                end do
            end do
        end do       

    do j = lo(2),hi(2)+1
        do i = lo(1),hi(1)
            do n=1,nspecies

                if(min(molarconc(i,j,n),molarconc(i,j-1,n)) .le. eepsilon) then !need to make this an input parameter?            
                    Gpy(i,j,1)=1.d0
                    Gpy(i,j,2)=0.d0
                    Gpy(i,j,3)=0.d0
                    Gpy(i,j,4)=1.d0

                end if

                end do
            end do
        end do       
    

    end subroutine limit_Gama_2d

  subroutine limit_Gama_3d(molarconc, ng_m,Gpx,Gpy,Gpz,ng_G, &
                                     lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:),ng_m,ng_G
    real(kind=dp_t), intent(in) :: molarconc(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:,:)
    real(kind=dp_t), intent(inout) :: Gpx(lo(1)-ng_G:,lo(2)-ng_G:,lo(3)-ng_G:,:)
    real(kind=dp_t), intent(inout) :: Gpy(lo(1)-ng_G:,lo(2)-ng_G:,lo(3)-ng_G:,:)
    real(kind=dp_t), intent(inout) :: Gpz(lo(1)-ng_G:,lo(2)-ng_G:,lo(3)-ng_G:,:)

    integer i,j,k,n
    real(kind=dp_t) eepsilon

    eepsilon = 0.0001

    do k = lo(3),hi(3)
    do j = lo(2),hi(2)
        do i = lo(1),hi(1)+1
            do n=1,nspecies

                if(min(molarconc(i,j,k,n),molarconc(i-1,j,k,n)) .le. eepsilon) then !need to make this an input parameter?            

                    Gpx(i,j,k,1)=1.d0
                    Gpx(i,j,k,2)=0.d0
                    Gpx(i,j,k,3)=0.d0
                    Gpx(i,j,k,4)=1.d0

                endif

                end do
            end do
        end do       
        end do       

    do k = lo(3),hi(3)
    do j = lo(2),hi(2)+1
        do i = lo(1),hi(1)
            do n=1,nspecies

                if(min(molarconc(i,j,k,n),molarconc(i,j-1,k,n)) .le. eepsilon) then !need to make this an input parameter?            
                    Gpy(i,j,k,1)=1.d0
                    Gpy(i,j,k,2)=0.d0
                    Gpy(i,j,k,3)=0.d0
                    Gpy(i,j,k,4)=1.d0

                end if

                end do
            end do
        end do       
        end do       

    do k = lo(3),hi(3)+1
    do j = lo(2),hi(2)
        do i = lo(1),hi(1)
            do n=1,nspecies

                if(min(molarconc(i,j,k,n),molarconc(i,j,k-1,n)) .le. eepsilon) then !need to make this an input parameter?            
                    Gpz(i,j,k,1)=1.d0
                    Gpz(i,j,k,2)=0.d0
                    Gpz(i,j,k,3)=0.d0
                    Gpz(i,j,k,4)=1.d0

                end if

                end do
            end do
        end do       
        end do       
    

    end subroutine limit_Gama_3d



end module diffusive_mass_fluxdiv_module
