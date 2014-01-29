module fluid_model_module

  use multifab_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use ml_layout_module
  use probin_multispecies_module
 
  implicit none

  private

  public :: fluid_model
  
contains
  
  subroutine fluid_model(mla,rho,rho_tot,molarconc,molmtot,D_MS,Gama,the_bc_level)

    type(ml_layout), intent(in   )  :: mla
    type(multifab),  intent(in   )  :: rho(:) 
    type(multifab),  intent(in   )  :: rho_tot(:) 
    type(multifab),  intent(in   )  :: molarconc(:) 
    type(multifab),  intent(in   )  :: molmtot(:) 
    type(multifab),  intent(inout)  :: D_MS(:)      ! MS diffusion constants 
    type(multifab),  intent(inout)  :: Gama(:)      ! Non-ideality coefficient 
    type(bc_level),  intent(in   )  :: the_bc_level(:)
 
    ! local variables
    integer :: lo(rho(1)%dim), hi(rho(1)%dim)
    integer :: n,i,ng,dm,nlevs,row,column

    ! assign pointers for multifabs to be passed
    real(kind=dp_t), pointer        :: dp(:,:,:,:)   ! for rho    
    real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for rho_tot
    real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for molarconc
    real(kind=dp_t), pointer        :: dp3(:,:,:,:)  ! for molmtot
    real(kind=dp_t), pointer        :: dp4(:,:,:,:)  ! for D_MS
    real(kind=dp_t), pointer        :: dp5(:,:,:,:)  ! for Gama

    dm    = mla%dim     ! dimensionality
    ng    = rho(1)%ng   ! number of ghost cells 
    nlevs = mla%nlevel  ! number of levels 
 
    ! loop over all boxes 
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp => dataptr(rho(n),i)
          dp1 => dataptr(rho_tot(n),i)
          dp2 => dataptr(molarconc(n),i)
          dp3 => dataptr(molmtot(n),i)
          dp4 => dataptr(D_MS(n),i)
          dp5 => dataptr(Gama(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          
          select case(dm)
          case (2)
             call compute_D_MSGama_2d(dp(:,:,1,:),dp1(:,:,1,1),dp2(:,:,1,:),&
                  dp3(:,:,1,1),dp4(:,:,1,:),dp5(:,:,1,:),ng,lo,hi) 
          case (3)
             call compute_D_MSGama_3d(dp(:,:,:,:),dp1(:,:,:,1),dp2(:,:,:,:),&
                  dp3(:,:,:,1),dp4(:,:,:,:),dp5(:,:,:,:),ng,lo,hi) 
          end select
       end do
    end do
  
  end subroutine fluid_model

  subroutine compute_D_MSGama_2d(rho,rho_tot,molarconc,molmtot,D_MS,Gama,ng,lo,hi)

    integer          :: lo(2), hi(2), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)        ! density; last dimension for species
    real(kind=dp_t)  :: rho_tot(lo(1)-ng:,lo(2)-ng:)      ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,:)  ! molar concentration 
    real(kind=dp_t)  :: molmtot(lo(1)-ng:,lo(2)-ng:)      ! total molar mass 
    real(kind=dp_t)  :: D_MS(lo(1)-ng:,lo(2)-ng:,:)       ! last dimension for nspecies^2
    real(kind=dp_t)  :: Gama(lo(1)-ng:,lo(2)-ng:,:)       ! last dimension for nspecies^2

    ! local varialbes; vectors and matrices to be used by D_MS, Gama 
    real(kind=dp_t), dimension(nspecies,nspecies) :: D_MS_local, Gama_local
    integer                                       :: n,i,j,row,column

    ! for specific box, now start loops over alloted cells 
    do j=lo(2)-ng,hi(2)+ng
       do i=lo(1)-ng,hi(1)+ng
        
          ! populate D_MS,Gama; for initial case doesn't change in each cell. 
          n=0; 
          do row=1, nspecies  
             do column=1, row-1
                n=n+1
                D_MS_local(row, column) = Dbar_in(n)
                D_MS_local(column, row) = D_MS_local(row, column) ! symmetric
                Gama_local(row, column) = 0.d0       
                Gama_local(column, row) = Gama_local(row, column) ! symmetric
             enddo
             D_MS_local(row, row) = 0.d0             ! self-diffusion is zero
             Gama_local(row, row) = 1.d0             ! set to unit matrix for time being
          enddo

          ! do the rank conversion 
          call set_Xij(D_MS(i,j,:), D_MS_local)
          call set_Xij(Gama(i,j,:), Gama_local)
 
       end do
    end do
   
    ! use contained (internal) subroutine to do the copy without doing index algebra:
    contains 
     subroutine set_Xij(Xout_ij, Xin_ij)
        real(kind=dp_t), dimension(nspecies,nspecies), intent(in)  :: Xin_ij
        real(kind=dp_t), dimension(nspecies,nspecies), intent(out) :: Xout_ij  
        Xout_ij = Xin_ij
     end subroutine 

  end subroutine compute_D_MSGama_2d
 
  subroutine compute_D_MSGama_3d(rho,rho_tot,molarconc,molmtot,D_MS,Gama,ng,lo,hi)
 
    integer          :: lo(3), hi(3), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)        ! density; last dimension for species
    real(kind=dp_t)  :: rho_tot(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)      ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  ! molar concentration; 
    real(kind=dp_t)  :: molmtot(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)      ! total molar mass 
    real(kind=dp_t)  :: D_MS(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)       ! last dimension for nspecies^2
    real(kind=dp_t)  :: Gama(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)       ! last dimension for nspecies^2

    ! local varialbes; vectors and matrices to be used by D_MS, Gama 
    real(kind=dp_t), dimension(nspecies,nspecies) :: D_MS_local, Gama_local
    integer                                       :: n,i,j,k,row,column

    ! for specific box, now start loops over alloted cells 
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng

             ! populate D_MS,Gama; for initial case doesn't change in each cell. 
             n=0; 
             do row=1, nspecies  
                do column=1, row-1
                   n=n+1
                   D_MS_local(row, column) = Dbar_in(n)
                   D_MS_local(column, row) = D_MS_local(row, column) ! symmetric
                   Gama_local(row, column) = 0.d0       
                   Gama_local(column, row) = Gama_local(row, column) ! symmetric
                enddo
                D_MS_local(row, row) = 0.d0             ! self-diffusion is zero
                Gama_local(row, row) = 1.d0             ! set to unit matrix for time being
             enddo

             ! do the rank conversion 
             call set_Xij(D_MS(i,j,k,:), D_MS_local)
             call set_Xij(Gama(i,j,k,:), Gama_local)
 
          end do
       end do
    end do
   
    ! use contained (internal) subroutine to do the copy without doing index algebra:
    contains 
     subroutine set_Xij(Xout_ij, Xin_ij)
        real(kind=dp_t), dimension(nspecies,nspecies), intent(in)  :: Xin_ij
        real(kind=dp_t), dimension(nspecies,nspecies), intent(out) :: Xout_ij  
        Xout_ij = Xin_ij
     end subroutine 

  end subroutine compute_D_MSGama_3d

end module fluid_model_module
