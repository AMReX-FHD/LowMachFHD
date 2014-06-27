module compute_mixture_properties_module
  ! The purpose of the fluid model is to provide concentration-dependent transport coefficients
  ! which should change with different problem types and number of species
  ! depending on the exact physical system being simulated. 

  use multifab_module
  use define_bc_module
  use ml_layout_module
  use bc_module
  use convert_stag_module
  use probin_common_module, only: molmass, prob_type, visc_coef
  use probin_multispecies_module, only: nspecies, is_ideal_mixture, Dbar, &
                                        Dtherm, H_offdiag, H_diag
 
  implicit none

  private

  public :: compute_mixture_properties, compute_eta, compute_kappa
    
contains
  
  subroutine compute_mixture_properties(mla,rho,rhotot,molarconc,molmtot,D_bar,D_therm,Hessian,Temp,the_bc_level)

    type(ml_layout), intent(in   )  :: mla
    type(multifab),  intent(in   )  :: rho(:) 
    type(multifab),  intent(in   )  :: rhotot(:) 
    type(multifab),  intent(in   )  :: molarconc(:) 
    type(multifab),  intent(in   )  :: molmtot(:) 
    type(multifab),  intent(inout)  :: D_bar(:)      ! MS diffusion constants 
    type(multifab),  intent(inout)  :: D_therm(:)    ! thermo diffusion constants 
    type(multifab),  intent(inout)  :: Hessian(:)    ! Non-ideality coefficient 
    type(multifab) , intent(in   )  :: Temp(:) 
    type(bc_level),  intent(in   )  :: the_bc_level(:)
 
    ! local variables
    integer :: lo(rho(1)%dim), hi(rho(1)%dim)
    integer :: n,i,ng,dm,nlevs

    ! assign pointers for multifabs to be passed
    real(kind=dp_t), pointer        :: dp(:,:,:,:)   ! for rho    
    real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for rhotot
    real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for molarconc
    real(kind=dp_t), pointer        :: dp3(:,:,:,:)  ! for molmtot
    real(kind=dp_t), pointer        :: dp4(:,:,:,:)  ! for D_bar
    real(kind=dp_t), pointer        :: dp5(:,:,:,:)  ! for D_therm 
    real(kind=dp_t), pointer        :: dp6(:,:,:,:)  ! for Hessian
    real(kind=dp_t), pointer        :: dp7(:,:,:,:)  ! for Temp

    dm    = mla%dim     ! dimensionality
    ng    = rho(1)%ng   ! number of ghost cells 
    nlevs = mla%nlevel  ! number of levels 
 
    ! loop over all boxes 
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp => dataptr(rho(n),i)
          dp1 => dataptr(rhotot(n),i)
          dp2 => dataptr(molarconc(n),i)
          dp3 => dataptr(molmtot(n),i)
          dp4 => dataptr(D_bar(n),i)
          dp5 => dataptr(D_therm(n),i)
          dp6 => dataptr(Hessian(n),i)
          dp7 => dataptr(Temp(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          
          select case(dm)
          case (2)
             call mixture_properties_mass_2d(dp(:,:,1,:),dp1(:,:,1,1),dp2(:,:,1,:),&
                  dp3(:,:,1,1),dp4(:,:,1,:),dp5(:,:,1,:),dp6(:,:,1,:),dp7(:,:,1,1),ng,lo,hi) 
          case (3)
             call mixture_properties_mass_3d(dp(:,:,:,:),dp1(:,:,:,1),dp2(:,:,:,:),&
                  dp3(:,:,:,1),dp4(:,:,:,:),dp5(:,:,:,:),dp6(:,:,:,:),dp7(:,:,:,1),ng,lo,hi) 
          end select
       end do
    end do
  
  end subroutine compute_mixture_properties
  
  subroutine mixture_properties_mass_2d(rho,rhotot,molarconc,molmtot,D_bar,D_therm,Hessian,Temp,ng,lo,hi)

    integer          :: lo(2), hi(2), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)        ! density; last dimension for species
    real(kind=dp_t)  :: rhotot(lo(1)-ng:,lo(2)-ng:)       ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,:)  ! molar concentration 
    real(kind=dp_t)  :: molmtot(lo(1)-ng:,lo(2)-ng:)      ! total molar mass 
    real(kind=dp_t)  :: D_bar(lo(1)-ng:,lo(2)-ng:,:)      ! last dimension for nspecies^2
    real(kind=dp_t)  :: D_therm(lo(1)-ng:,lo(2)-ng:,:)    ! last dimension for nspecies
    real(kind=dp_t)  :: Hessian(lo(1)-ng:,lo(2)-ng:,:)    ! last dimension for nspecies^2
    real(kind=dp_t)  :: Temp(lo(1)-ng:,lo(2)-ng:)         ! Temperature 

    ! local varialbes
    integer          :: i,j

    ! for specific box, now start loops over alloted cells 
    do j=lo(2)-ng,hi(2)+ng
       do i=lo(1)-ng,hi(1)+ng
       
          call mixture_properties_mass_local(rho(i,j,:),rhotot(i,j),molarconc(i,j,:),&
                                             molmtot(i,j),D_bar(i,j,:),D_therm(i,j,:),&
                                             Hessian(i,j,:),Temp(i,j))
       end do
    end do
   
  end subroutine mixture_properties_mass_2d

  subroutine mixture_properties_mass_3d(rho,rhotot,molarconc,molmtot,D_bar,D_therm,Hessian,Temp,ng,lo,hi)
 
    integer          :: lo(3), hi(3), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)       ! density; last dimension for species
    real(kind=dp_t)  :: rhotot(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)      ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:) ! molar concentration; 
    real(kind=dp_t)  :: molmtot(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)     ! total molar mass 
    real(kind=dp_t)  :: D_bar(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)     ! last dimension for nspecies^2
    real(kind=dp_t)  :: D_therm(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)   ! last dimension for nspecies
    real(kind=dp_t)  :: Hessian(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)   ! last dimension for nspecies^2
    real(kind=dp_t)  :: Temp(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)        ! Temperature 
 
    ! local varialbes
    integer          :: i,j,k

    ! for specific box, now start loops over alloted cells 
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng

             call mixture_properties_mass_local(rho(i,j,k,:),rhotot(i,j,k),molarconc(i,j,k,:),&
                                                molmtot(i,j,k),D_bar(i,j,k,:),D_therm(i,j,k,:),&
                                                Hessian(i,j,k,:),Temp(i,j,k))
          end do
       end do
    end do
   
  end subroutine mixture_properties_mass_3d

  ! The default case should be to simply set D_bar, D_therm and H to constants, read from the input file
  subroutine mixture_properties_mass_local(rho,rhotot,molarconc,molmtot,D_bar,D_therm,Hessian,Temp)
   
    real(kind=dp_t), intent(in)   :: rho(nspecies)        
    real(kind=dp_t), intent(in)   :: rhotot
    real(kind=dp_t), intent(in)   :: molarconc(nspecies)
    real(kind=dp_t), intent(in)   :: molmtot
    real(kind=dp_t), intent(out)  :: D_bar(nspecies,nspecies) 
    real(kind=dp_t), intent(out)  :: D_therm(nspecies) 
    real(kind=dp_t), intent(out)  :: Hessian(nspecies,nspecies)
    real(kind=dp_t), intent(in)   :: Temp
 
    ! local variables
    integer :: n,row,column
    real(kind=dp_t) :: massfrac(nspecies), c_loc

    ! Local values of transport and thermodynamic coefficients (functions of composition!):
    real(kind=dp_t), dimension(nspecies*(nspecies-1)/2) :: D_bar_local, H_offdiag_local ! off-diagonal components of symmetric matrices
    real(kind=dp_t), dimension(nspecies)    :: H_diag_local ! Diagonal component
    
    massfrac = rho/rhotot;
    
    select case (abs(prob_type))
    case (9)

       ! This is where formula for chi as a function of concentration goes
       ! We assume nspecies=2
       ! Dbar(1) = chi0 in the binary notation
       if (nspecies .ne. 2) then
          call bl_error("mixture_properties_mass_local assumes nspecies=2 if prob_type=9 (water-glycerol)")
       end if
       
       call chi_water_glycerol(D_bar_local(1), rho, rhotot, Temp)
       
    case default

       D_bar_local(1:nspecies*(nspecies-1)/2) = Dbar(1:nspecies*(nspecies-1)/2) ! Keep it constant

    end select
    
    ! For now we only encode constant Hessian matrices since we do not have any thermodynamic models coded up (Wilson, NTLR, UNIQUAC, etc.)
    H_diag_local(1:nspecies) = H_diag(1:nspecies)
    H_offdiag_local(1:nspecies*(nspecies-1)/2) = H_offdiag(1:nspecies*(nspecies-1)/2)
    
    ! Complete the process by filling the matrices using generic formulae -- this part should not change
    ! populate D_bar and Hessian matrix 
    n=0; 
    do row=1, nspecies  
       do column=1, row-1
          n=n+1
          D_bar(row, column) = D_bar_local(n)                ! SM-diffcoeff's read from input
          D_bar(column, row) = D_bar(row, column)     ! symmetric
          
          if(.not. is_ideal_mixture) then
             Hessian(row, column) = H_offdiag_local(n)         ! positive semidefinite matrix read from input
             Hessian(column, row) = Hessian(row,column)  ! Hessian is symmetric
          end if
       end do
       
       ! populate diagonals 
       D_bar(row, row) = 0.d0           ! as self-diffusion is zero
       D_therm(row)    = Dtherm(row)    ! thermal diffcoeff's read from input
       if(.not. is_ideal_mixture) then 
          Hessian(row, row) = H_diag_local(row)     
       else 
          Hessian = 0.d0     ! set matrix to zero for ideal-mixture
       end if
    
    end do
    
  end subroutine mixture_properties_mass_local

  subroutine compute_eta(mla,eta,eta_ed,rho,rhotot,Temp,pres,dx,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: eta(:)
    type(multifab) , intent(inout) :: eta_ed(:,:) ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(in   ) :: rho(:)
    type(multifab) , intent(in   ) :: rhotot(:)
    type(multifab) , intent(in   ) :: Temp(:)
    type(multifab) , intent(in   ) :: pres(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    integer :: nlevs,dm,i,n,ng_e,ng_r,ng_m,ng_t,ng_p
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: ep(:,:,:,:)
    real(kind=dp_t), pointer :: rp(:,:,:,:)
    real(kind=dp_t), pointer :: mp(:,:,:,:)
    real(kind=dp_t), pointer :: tp(:,:,:,:)
    real(kind=dp_t), pointer :: pp(:,:,:,:)

    ng_e  = eta(1)%ng
    ng_r  = rho(1)%ng
    ng_m = rhotot(1)%ng
    ng_t  = Temp(1)%ng
    ng_p  = pres(1)%ng

    nlevs = mla%nlevel
    dm = mla%dim

    do n=1,nlevs
       do i=1,nfabs(eta(n))
          ep  => dataptr(eta(n), i)
          rp  => dataptr(rho(n), i)
          mp => dataptr(rhotot(n), i)
          tp  => dataptr(Temp(n), i)
          pp  => dataptr(pres(n), i)
          lo = lwb(get_box(eta(n), i))
          hi = upb(get_box(eta(n), i))
          select case (dm)
          case (2)
             call compute_eta_2d(ep(:,:,1,1),ng_e,rp(:,:,1,:),ng_r,mp(:,:,1,1),ng_m, &
                                 tp(:,:,1,1),ng_t,pp(:,:,1,1),ng_p,lo,hi,dx(n,:))
          case (3)
             call compute_eta_3d(ep(:,:,:,1),ng_e,rp(:,:,:,:),ng_r,mp(:,:,:,1),ng_m, &
                                 tp(:,:,:,1),ng_t,pp(:,:,:,1),ng_p,lo,hi,dx(n,:))
          end select
       end do
    end do

    if (dm .eq. 2) then
       call average_cc_to_node(nlevs,eta,eta_ed(:,1),1,tran_bc_comp,1,the_bc_level)
    else if (dm .eq. 3) then
       call average_cc_to_edge(nlevs,eta,eta_ed,1,tran_bc_comp,1,the_bc_level)
    end if

  end subroutine compute_eta

  subroutine compute_eta_2d(eta,ng_e,rho,ng_r,rhotot,ng_m,Temp,ng_t,pres,ng_p,lo,hi,dx)

    ! compute eta in valid AND ghost regions
    ! the ghost cells for rho, Temp, etc., have already been filled properly

    integer        , intent(in   ) :: lo(:), hi(:), ng_e, ng_r, ng_m, ng_t, ng_p
    real(kind=dp_t), intent(inout) ::    eta(lo(1)-ng_e:,lo(2)-ng_e:)
    real(kind=dp_t), intent(inout) ::    rho(lo(1)-ng_r:,lo(2)-ng_r:,:)
    real(kind=dp_t), intent(inout) :: rhotot(lo(1)-ng_m:,lo(2)-ng_m:)
    real(kind=dp_t), intent(inout) ::   Temp(lo(1)-ng_t:,lo(2)-ng_t:)
    real(kind=dp_t), intent(inout) ::   pres(lo(1)-ng_p:,lo(2)-ng_p:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! local
    integer :: i,j

    select case (abs(prob_type))
    case (9)
       
       do j=lo(2)-ng_e,hi(2)+ng_e
          do i=lo(1)-ng_e,hi(1)+ng_e

             call eta_water_glycerol(eta(i,j),rho(i,j,:),rhotot(i,j),Temp(i,j))

          end do
       end do

    case default

       eta = visc_coef

    end select

  end subroutine compute_eta_2d

  subroutine compute_eta_3d(eta,ng_e,rho,ng_r,rhotot,ng_m,Temp,ng_t,pres,ng_p,lo,hi,dx)

    ! compute eta in valid AND ghost regions
    ! the ghost cells for rho, Temp, etc., have already been filled properly

    integer        , intent(in   ) :: lo(:), hi(:), ng_e, ng_r, ng_m, ng_t, ng_p
    real(kind=dp_t), intent(inout) ::    eta(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(inout) ::    rho(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:,:)
    real(kind=dp_t), intent(inout) :: rhotot(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real(kind=dp_t), intent(inout) ::   Temp(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
    real(kind=dp_t), intent(inout) ::   pres(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! local
    integer :: i,j,k

    select case (abs(prob_type))
    case (9)

       do k=lo(3)-ng_e,hi(3)+ng_e
          do j=lo(2)-ng_e,hi(2)+ng_e
             do i=lo(1)-ng_e,hi(1)+ng_e

                call eta_water_glycerol(eta(i,j,k),rho(i,j,k,:),rhotot(i,j,k),Temp(i,j,k))

             end do
          end do
       end do

    case default

       eta = visc_coef

    end select

  end subroutine compute_eta_3d

  !=================================================
  ! Water-glycerol mixtures near room temperature
  !=================================================
  
  subroutine chi_water_glycerol(chi,rho,rhotot,Temp) ! This only works for room temperature for now

    real(kind=dp_t), intent(inout) :: chi
    real(kind=dp_t), intent(in   ) :: rho(:)
    real(kind=dp_t), intent(in   ) :: rhotot
    real(kind=dp_t), intent(in   ) :: Temp

    ! local
    real(kind=dp_t) :: c_loc
    
    ! mass fraction of glycerol
    c_loc = rho(1)/rhotot

    ! chi = chi0 * rational function
    chi = Dbar(1)*(1.024d0-1.001692692d0*c_loc)/(1.d0+0.6632641981d0*c_loc)
        
  end subroutine chi_water_glycerol

  subroutine eta_water_glycerol(eta,rho,rhotot,Temp)

    real(kind=dp_t), intent(inout) :: eta
    real(kind=dp_t), intent(in   ) :: rho(:)
    real(kind=dp_t), intent(in   ) :: rhotot
    real(kind=dp_t), intent(in   ) :: Temp
    
    ! local
    real(kind=dp_t) :: c_loc, nu_g, nu_w, x, T

    ! mass fraction of glycerol
    c_loc = rho(1)/rhotot

    ! convert temperature to Celsius
    T = Temp - 273.d0 

    ! viscosities of pure glycerol and water
    nu_g = exp(9.09309 - 0.11397*T + 0.00054*T**2)
    nu_w = exp(0.55908 - 0.03051*T + 0.00015*T**2)

    x = c_loc*(1.d0 + (1.d0-c_loc)*(-0.727770 - 0.04943*c_loc - 1.2038*c_loc**2))

    ! visc_coef is the scaling prefactor (for unit conversion or non-unity scaling)
    eta = visc_coef*exp(-x*(-log(nu_g)+log(nu_w)))*nu_w*rhotot

  end subroutine eta_water_glycerol

  subroutine compute_kappa(mla,kappa)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: kappa(:)

    integer :: n,nlevs

    nlevs = mla%nlevel

    do n=1,nlevs
       call multifab_setval(kappa(n), 1.d0, all=.true.)
    end do

  end subroutine compute_kappa

end module compute_mixture_properties_module
