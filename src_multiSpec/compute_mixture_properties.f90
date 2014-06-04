module compute_mixture_properties_module

  use multifab_module
  use define_bc_module
  use ml_layout_module
  use probin_common_module, only: molmass, prob_type
  use probin_multispecies_module, only: nspecies, is_ideal_mixture, Dbar, &
                                        Dtherm, H_offdiag, H_diag
 
  implicit none

  private

  public :: compute_mixture_properties ! DONEV: Why would mixture_properties_mass_local be public?
  
  ! Donev:
  ! The purpose of the fluid model is to provide concentration-dependent transport coefficients
  ! which should change with different problem types and number of species
  ! depending on the exact physical system being simulated. 
  ! A lot of the stuff below is general and thus belongs with the other general routines such as compute_chi
  ! In particular the main piece below computes Gamma (the rest is fluff code) from H
  ! call H "Hessian" to make it clear it is symmetric
  ! This code should be separated from here and moved to convert_mass_variables.f90
  ! (which should be renamed to something more meaningful, like mass_flux_utilities.f90)
  
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
    integer :: n,i,ng,dm,nlevs,row,column

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

  ! Donev: This is the key routine here
  ! It should have a case statement, in which different things are done depending on prob_type (Andy can organize that part)
  ! For now the default case should be to simply set D_bar, D_therm and H to constants, read from the input file, as done below
  ! Rename this routine mixture_properties_mass (for now) to make it clear what this does
  
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
    integer                       :: n,row,column
    real(kind=dp_t) :: massfrac(nspecies)

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
          call bl_error("mixture_properties_mass_local assumes nspecies=2")
       end if
       
       D_bar_local(1) = Dbar(1) ! * function(massfrac(1))

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

end module compute_mixture_properties_module
