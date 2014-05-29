module fluid_model_module

  use multifab_module
  use define_bc_module
  use ml_layout_module
  use probin_common_module, only: molmass
  use probin_multispecies_module, only: nspecies, is_ideal_mixture, Dbar, &
      Dtherm, H_offdiag, H_diag
 
  implicit none

  private

  public :: fluid_model, compute_D_MSGama_local
  
contains
  
  subroutine fluid_model(mla,rho,rhotot,molarconc,molmtot,D_MS,D_therm,Gama,the_bc_level)

    type(ml_layout), intent(in   )  :: mla
    type(multifab),  intent(in   )  :: rho(:) 
    type(multifab),  intent(in   )  :: rhotot(:) 
    type(multifab),  intent(in   )  :: molarconc(:) 
    type(multifab),  intent(in   )  :: molmtot(:) 
    type(multifab),  intent(inout)  :: D_MS(:)      ! MS diffusion constants 
    type(multifab),  intent(inout)  :: D_therm(:)   ! thermo diffusion constants 
    type(multifab),  intent(inout)  :: Gama(:)      ! Non-ideality coefficient 
    type(bc_level),  intent(in   )  :: the_bc_level(:)
 
    ! local variables
    integer :: lo(rho(1)%dim), hi(rho(1)%dim)
    integer :: n,i,ng,dm,nlevs,row,column

    ! assign pointers for multifabs to be passed
    real(kind=dp_t), pointer        :: dp(:,:,:,:)   ! for rho    
    real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for rhotot
    real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for molarconc
    real(kind=dp_t), pointer        :: dp3(:,:,:,:)  ! for molmtot
    real(kind=dp_t), pointer        :: dp4(:,:,:,:)  ! for D_MS
    real(kind=dp_t), pointer        :: dp5(:,:,:,:)  ! for D_therm 
    real(kind=dp_t), pointer        :: dp6(:,:,:,:)  ! for Gama 

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
          dp4 => dataptr(D_MS(n),i)
          dp5 => dataptr(D_therm(n),i)
          dp6 => dataptr(Gama(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          
          select case(dm)
          case (2)
             call compute_D_MSGama_2d(dp(:,:,1,:),dp1(:,:,1,1),dp2(:,:,1,:),&
                  dp3(:,:,1,1),dp4(:,:,1,:),dp5(:,:,1,:),dp6(:,:,1,:),ng,lo,hi) 
          case (3)
             call compute_D_MSGama_3d(dp(:,:,:,:),dp1(:,:,:,1),dp2(:,:,:,:),&
                  dp3(:,:,:,1),dp4(:,:,:,:),dp5(:,:,:,:),dp6(:,:,:,:),ng,lo,hi) 
          end select
       end do
    end do
  
  end subroutine fluid_model
  
  subroutine compute_D_MSGama_2d(rho,rhotot,molarconc,molmtot,D_MS,D_therm,Gama,ng,lo,hi)

    integer          :: lo(2), hi(2), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)        ! density; last dimension for species
    real(kind=dp_t)  :: rhotot(lo(1)-ng:,lo(2)-ng:)      ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,:)  ! molar concentration 
    real(kind=dp_t)  :: molmtot(lo(1)-ng:,lo(2)-ng:)      ! total molar mass 
    real(kind=dp_t)  :: D_MS(lo(1)-ng:,lo(2)-ng:,:)       ! last dimension for nspecies^2
    real(kind=dp_t)  :: D_therm(lo(1)-ng:,lo(2)-ng:,:)    ! last dimension for nspecies
    real(kind=dp_t)  :: Gama(lo(1)-ng:,lo(2)-ng:,:)       ! last dimension for nspecies^2

    ! local varialbes
    integer          :: i,j

    ! for specific box, now start loops over alloted cells 
    do j=lo(2)-ng,hi(2)+ng
       do i=lo(1)-ng,hi(1)+ng
       
          call compute_D_MSGama_local(rho(i,j,:),rhotot(i,j),molarconc(i,j,:),&
                                      molmtot(i,j),D_MS(i,j,:),D_therm(i,j,:),Gama(i,j,:))
       end do
    end do
   
  end subroutine compute_D_MSGama_2d

  subroutine compute_D_MSGama_3d(rho,rhotot,molarconc,molmtot,D_MS,D_therm,Gama,ng,lo,hi)
 
    integer          :: lo(3), hi(3), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)        ! density; last dimension for species
    real(kind=dp_t)  :: rhotot(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)      ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  ! molar concentration; 
    real(kind=dp_t)  :: molmtot(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)      ! total molar mass 
    real(kind=dp_t)  :: D_MS(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)       ! last dimension for nspecies^2
    real(kind=dp_t)  :: D_therm(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)    ! last dimension for nspecies
    real(kind=dp_t)  :: Gama(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)       ! last dimension for nspecies^2

    ! local varialbes
    integer          :: i,j,k

    ! for specific box, now start loops over alloted cells 
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng

             call compute_D_MSGama_local(rho(i,j,k,:),rhotot(i,j,k),molarconc(i,j,k,:),&
                                         molmtot(i,j,k),D_MS(i,j,k,:),D_therm(i,j,k,:),Gama(i,j,k,:))
          end do
       end do
    end do
   
  end subroutine compute_D_MSGama_3d

  subroutine compute_D_MSGama_local(rho,rhotot,molarconc,molmtot,D_MS,D_therm,Gama)
   
    real(kind=dp_t), intent(in)   :: rho(nspecies)        
    real(kind=dp_t), intent(in)   :: rhotot
    real(kind=dp_t), intent(in)   :: molarconc(nspecies)
    real(kind=dp_t), intent(in)   :: molmtot
    real(kind=dp_t), intent(out)  :: D_MS(nspecies,nspecies) 
    real(kind=dp_t), intent(out)  :: D_therm(nspecies) 
    real(kind=dp_t), intent(out)  :: Gama(nspecies,nspecies)
 
    ! local variables
    integer                                       :: n,row,column
    real(kind=dp_t), dimension(nspecies,nspecies) :: H, I, X_xxT   

    ! free the memory
    H=0; I=0; X_xxT=0;

    ! populate D_MS, H, I and X_xxT where  X = molmtot*W*M^(-1) 
    n=0; 
    do row=1, nspecies  
       do column=1, row-1
          n=n+1
          D_MS(row, column) = Dbar(n)                              ! SM-diffcoeff's read from input
          D_MS(column, row) = D_MS(row, column)                    ! symmetric
          I(row, column)    = 0.d0      
          I(column, row)    = I(row, column) ! symmetric
          
          if(is_ideal_mixture .eqv. .false.) then
             X_xxT(row,column)   = -molarconc(row)*molarconc(column)  ! form x*transpose(x) off diagonals 
             X_xxT(column, row)  = X_xxT(row, column)                 ! symmetric
             H(row, column) = H_offdiag(n)           ! positive semidefinite matrix read from input
             H(column, row) = H(row,column)          ! H is symmetric 
          end if
       end do
       
       ! populate diagonals 
       D_MS(row, row) = 0.d0           ! as self-diffusion is zero
       D_therm(row)   = Dtherm(row)    ! thermal diffcoeff's read from input
       I(row, row) = 1.d0        ! unit matrix for ideal mixture
       if(is_ideal_mixture .eqv. .false.) then
          X_xxT(row,row) = molmtot*rho(row)/(rhotot*molmass(row)) - molarconc(row)**2 
          H(row, row)    = H_diag(n)      
       end if
    end do

    if(is_ideal_mixture) then
       Gama = I
    else 
       Gama = I + matmul(X_xxT, H)     ! non-ideal mixture
    end if

  end subroutine compute_D_MSGama_local

end module fluid_model_module
