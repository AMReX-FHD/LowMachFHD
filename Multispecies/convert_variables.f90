module convert_variables_module

  use multifab_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use ml_layout_module
  use probin_multispecies_module
  use matrix_utilities 
  use F95_LAPACK

  implicit none

  private

  public :: convert_cons_to_prim, compute_BinvGamma, compute_coefficient

contains

  subroutine convert_cons_to_prim(mla,rho,rho_tot,molarconc,mass,molmtot,the_bc_level)
   
   type(ml_layout), intent(in   )  :: mla
   type(multifab) , intent(in)     :: rho(:) 
   type(multifab) , intent(inout)  :: rho_tot(:) 
   type(multifab) , intent(inout)  :: molarconc(:) 
   real(kind=dp_t), intent(in   )  :: mass(:)
   type(multifab) , intent(inout)  :: molmtot(:) 
   type(bc_level) , intent(in   )  :: the_bc_level(:)

   ! local variables
   integer :: lo(rho(1)%dim), hi(rho(1)%dim)
   integer :: n,i,ng,dm,nlevs
 
   ! pointer for rho(nspecies), rho_tot(1), molarconc(nspecies) 
   real(kind=dp_t), pointer        :: dp(:,:,:,:)   ! for rho    
   real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for rho_tot
   real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for molar concentrations
   real(kind=dp_t), pointer        :: dp3(:,:,:,:)  ! for total molar concentrations

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
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          
          select case(dm)
          case (2)
             call compute_molconc_rhotot_2d(dp(:,:,1,:),dp1(:,:,1,1),&
                           dp2(:,:,1,:),mass(:),dp3(:,:,1,1),ng,lo,hi) 
          case (3)
             call compute_molconc_rhotot_3d(dp(:,:,:,:),dp1(:,:,:,1),&
                           dp2(:,:,:,:),mass(:),dp3(:,:,:,1),ng,lo,hi) 
          end select
       end do
    end do

  end subroutine convert_cons_to_prim

  subroutine compute_molconc_rhotot_2d(rho,rho_tot,molarconc,mass,molmtot,ng,lo,hi)
 
    integer          :: lo(2), hi(2), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)       ! density- last dim for #species
    real(kind=dp_t)  :: rho_tot(lo(1)-ng:,lo(2)-ng:)     ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,:) ! molar concentration
    real(kind=dp_t)  :: mass(:)                          ! species molar mass 
    real(kind=dp_t)  :: molmtot(lo(1)-ng:,lo(2)-ng:)     ! total molar mass 
    real(kind=dp_t), dimension(nspecies) :: W            ! mass fraction w_i = rho_i/rho 
    
    ! local variables
    integer          :: i,j,n
    real(kind=dp_t)  :: Sum_woverm, rho_tot_local
    
    ! for specific box, now start loops over alloted cells    
    do j=lo(2)-ng, hi(2)+ng
       do i=lo(1)-ng, hi(1)+ng

          ! calculate total density inside each cell
          rho_tot_local=0.d0 
          do n=1, nspecies  
             rho_tot_local = rho_tot_local + rho(i,j,n)
          enddo         
          rho_tot(i,j) = rho_tot_local
 
          ! calculate mass fraction and total molar mass (1/m=Sum(w_i/m_i))
          Sum_woverm=0.d0
          do n=1, nspecies  
             W(n) = rho(i,j,n)/rho_tot(i,j)
             Sum_woverm = Sum_woverm + W(n)/mass(n)
          enddo
          molmtot(i,j) = 1.0d0/Sum_woverm 
  
          ! calculate molar concentrations in each cell (x_i=m*w_i/m_i) 
          do n=1, nspecies 
             molarconc(i,j,n) = molmtot(i,j)*W(n)/mass(n)
          enddo

       enddo
    enddo
 
  end subroutine compute_molconc_rhotot_2d

  subroutine compute_molconc_rhotot_3d(rho,rho_tot,molarconc,mass,molmtot,ng,lo,hi)
 
    integer          :: lo(3), hi(3), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)       ! density- last dim for #species
    real(kind=dp_t)  :: rho_tot(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)     ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:) ! molar concentration
    real(kind=dp_t)  :: mass(:)                                    ! species molar mass 
    real(kind=dp_t)  :: molmtot(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)     ! total molar mass 
    real(kind=dp_t), dimension(nspecies) :: W                      ! mass fraction w_i = rho_i/rho 
    
    ! local variables
    integer          :: i,j,k,n
    real(kind=dp_t)  :: Sum_woverm, rho_tot_local
    
    ! for specific box, now start loops over alloted cells    
    do k=lo(3)-ng, hi(3)+ng
       do j=lo(2)-ng, hi(2)+ng
          do i=lo(1)-ng, hi(1)+ng

             ! calculate total density inside each cell
             rho_tot_local=0.d0 
             do n=1, nspecies  
                rho_tot_local = rho_tot_local + rho(i,j,k,n)
             enddo         
             rho_tot(i,j,k) = rho_tot_local

             ! calculate mass fraction and total molar mass (1/m=Sum(w_i/m_i))
             Sum_woverm=0.d0
             do n=1, nspecies  
                W(n) = rho(i,j,k,n)/rho_tot(i,j,k)
                Sum_woverm = Sum_woverm + W(n)/mass(n)
             enddo
             molmtot(i,j,k) = 1.0d0/Sum_woverm 
  
             ! calculate molar concentrations in each cell (x_i=m*w_i/m_i) 
             do n=1, nspecies 
                molarconc(i,j,k,n) = molmtot(i,j,k)*W(n)/mass(n)
             enddo
        
          enddo
       enddo
    enddo
 
  end subroutine compute_molconc_rhotot_3d

  subroutine compute_BinvGamma(mla,rho,rho_tot,molarconc,BinvGamma,Dbar,Gama,mass,molmtot,the_bc_level)
   
    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: rho(:)
    type(multifab) , intent(in   )  :: rho_tot(:) 
    type(multifab) , intent(in   )  :: molarconc(:) 
    type(multifab) , intent(inout)  :: BinvGamma(:) 
    real(kind=dp_t), intent(in   )  :: Dbar(:,:) 
    real(kind=dp_t), intent(in   )  :: Gama(:,:) 
    real(kind=dp_t), intent(in   )  :: mass(:)
    type(multifab) , intent(in   )  :: molmtot(:) 
    type(bc_level) , intent(in   )  :: the_bc_level(:)

    ! local variables
    integer :: lo(rho(1)%dim), hi(rho(1)%dim)
    integer :: n,i,ng,dm,nlevs
 
    ! pointer for rho(nspecies), molarconc(nspecies) 
    real(kind=dp_t), pointer        :: dp(:,:,:,:)   ! for rho    
    real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for rho_tot
    real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for molarconc
    real(kind=dp_t), pointer        :: dp3(:,:,:,:)  ! for BinvGamma 
    real(kind=dp_t), pointer        :: dp4(:,:,:,:)  ! for molmtot

    dm = mla%dim        ! dimensionality
    ng = rho(1)%ng      ! number of ghost cells 
    nlevs = mla%nlevel  ! number of levels 
 
    ! loop over all boxes 
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp => dataptr(rho(n),i)
          dp1 => dataptr(rho_tot(n),i)
          dp2 => dataptr(molarconc(n),i)
          dp3 => dataptr(BinvGamma(n),i)
          dp4 => dataptr(molmtot(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          
          select case(dm)
          case (2)
             call compute_BinvGamma_2d(dp(:,:,1,:),dp1(:,:,1,1),dp2(:,:,1,:),&
                  dp3(:,:,1,:),Dbar(:,:),Gama(:,:),mass(:),dp4(:,:,1,1),ng,lo,hi) 
          case (3)
             call compute_BinvGamma_3d(dp(:,:,:,:),dp1(:,:,:,1),dp2(:,:,:,:),&
                  dp3(:,:,:,:),Dbar(:,:),Gama(:,:),mass(:),dp4(:,:,:,1),ng,lo,hi) 
          end select
       end do
    end do

   end subroutine compute_BinvGamma

   subroutine compute_BinvGamma_2d(rho,rho_tot,molarconc,BinvGamma,Dbar,Gama,mass,molmtot,ng,lo,hi)

    integer          :: lo(2), hi(2), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)        ! density; last dimension for species
    real(kind=dp_t)  :: rho_tot(lo(1)-ng:,lo(2)-ng:)      ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,:)  ! molar concentration; 
    real(kind=dp_t)  :: BinvGamma(lo(1)-ng:,lo(2)-ng:,:)  ! B^(-1)*Gamma; last dimension for nspecies^2
    real(kind=dp_t)  :: Dbar(:,:)                         ! SM diffusion constants 
    real(kind=dp_t)  :: Gama(:,:)                         ! non-ideality coefficient 
    real(kind=dp_t)  :: mass(:)                           ! species molar mass 
    real(kind=dp_t)  :: molmtot(lo(1)-ng:,lo(2)-ng:)      ! total molar mass 

    ! local variables
    integer          :: i,j,row,column
    real(kind=dp_t)  :: tolerance, Sum_knoti              ! tolerance set for pinverse 

    ! vectors and matrices to be used by LAPACK 
    real(kind=dp_t), dimension(nspecies,nspecies) :: Bijprime, Lonsager
    real(kind=dp_t), dimension(nspecies,nspecies) :: BdagGamma,BdagCapW 
    real(kind=dp_t), dimension(nspecies)          :: W 

    tolerance = 1e-13

    ! for specific box, now start loops over alloted cells 
    do j=lo(2)-ng,hi(2)+ng
       do i=lo(1)-ng,hi(1)+ng
        
          ! free up the memory  
          Bijprime=0.d0
          BdagGamma = 0.d0
          W=0.d0
          Lonsager=0.d0         
 
          ! change 0 with tolerance to prevent division by zero in case species
          ! density, molar concentration or total density = 0. 
          do row=1, nspecies
             if(molarconc(i,j,row) .lt. tolerance) then
                molarconc(i,j,row) = tolerance
                rho(i,j,row)       = tolerance*rho_tot(i,j)                
             endif
          enddo

          ! compute Bijprime matrix and massfraction W_i = rho_i/rho; molarconc is 
          ! expressed in terms of molmtot,mi,rhotot etc. 
          do row=1, nspecies  
             do column=1, row-1
                Bijprime(row, column) = rho(i,j,row)*molmtot(i,j)**2/(mass(row)* &
                                        mass(column)*Dbar(row, column)*rho_tot(i,j)**2) 
                Bijprime(column, row) = rho(i,j,column)*molmtot(i,j)**2/(mass(row)* &
                                        mass(column)*Dbar(column, row)*rho_tot(i,j)**2) 
             enddo
            
             Sum_knoti=0.d0 
             do column=1, nspecies
                if (column.ne.row) then
                   Sum_knoti = Sum_knoti - rho(i,j,column)/(mass(column)*Dbar(row,column))
                endif
                Bijprime(row, row) = Sum_knoti*molmtot(i,j)**2/(mass(row)*rho_tot(i,j)**2)
             enddo
             
             W(row) = rho(i,j,row)/rho_tot(i,j)
          enddo

          ! compute Binverse*Gamma and other relevant matrices
          call compute_BdagGamma(Bijprime(:,:),BdagGamma(:,:),BdagCapW(:,:),Gama(:,:),W(:),tolerance)

          ! compute Onsager matrix L
          Lonsager = -rho_tot(i,j)*Temp*BdagCapW/Press

          !print*, i,j,"mol1=",molarconc(i,j,1),"mol2=",molarconc(i,j,2)
          if(.false.) then
            if(i.eq.7 .and. j.eq.14) then
               print*, "rho1=",rho(i,j,1),"rho2=",rho(i,j,2),"rho1+rho2=",rho_tot(i,j)
               print*, "rho1/rho=",W(1),"rho2/rho=",W(2),"w1+w2=",W(1)+W(2)
               print*, "mol1=",molarconc(i,j,1),"mol2=",molarconc(i,j,2),"mol1+mol2=",molmtot(i,j)
              do row=1, nspecies
                do column=1, nspecies
                   print*, BdagGamma(row, column)/rho_tot(i,j) 
                   !print*, BdagGamma(row, column)
                   !print*, Lonsager(row, column)
                   !print*, Gama(row, column)
                enddo
                !print*, Bdagw(row) 
                print*, '' 
              enddo
            endif
          endif

          ! do the rank conversion 
          call set_Bij(BinvGamma(i,j,:), BdagGamma)
 
      end do
    end do
   
    ! use contained (internal) subroutine to do the copy without doing index algebra:
    contains 
     subroutine set_Bij(BinvGamma_ij, BdagGamma_ij)
        real(kind=dp_t), dimension(nspecies,nspecies), intent(in)  :: BdagGamma_ij
        real(kind=dp_t), dimension(nspecies,nspecies), intent(out) :: BinvGamma_ij  
        BinvGamma_ij = BdagGamma_ij
     end subroutine 

  end subroutine compute_BinvGamma_2d

  subroutine compute_BinvGamma_3d(rho,rho_tot,molarconc,BinvGamma,Dbar,Gama,mass,molmtot,ng,lo,hi)
   
    integer          :: lo(3), hi(3), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)        ! density; last dimension for species
    real(kind=dp_t)  :: rho_tot(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)      ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  ! molar concentration; 
    real(kind=dp_t)  :: BinvGamma(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  ! B^(-1)*Gamma; last dimension for nspecies^2
    real(kind=dp_t)  :: Dbar(:,:)                                   ! SM diffusion constants 
    real(kind=dp_t)  :: Gama(:,:)                                   ! non-ideality coefficient 
    real(kind=dp_t)  :: mass(:)                                     ! species molar mass 
    real(kind=dp_t)  :: molmtot(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)      ! total molar mass 
    
    ! local variables
    integer          :: i,j,k,row,column
    real(kind=dp_t)  :: tolerance, Sum_knoti   

    ! vectors and matrices to be used by LAPACK 
    real(kind=dp_t), dimension(nspecies,nspecies) :: Bijprime, Lonsager
    real(kind=dp_t), dimension(nspecies,nspecies) :: BdagGamma,BdagCapW 
    real(kind=dp_t), dimension(nspecies)          :: W 

    tolerance = 1e-13

    ! for specific box, now start loops over alloted cells 
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
        
             ! free up the memory  
             Bijprime=0.d0
             BdagGamma = 0.d0
             W=0.d0
             Lonsager=0.d0         
 
             ! change 0 with tolerance to prevent division by zero in case species
             ! density, molar concentration or total density = 0. 
             do row=1, nspecies
                if(molarconc(i,j,k,row) .lt. tolerance) then
                   molarconc(i,j,k,row) = tolerance
                   rho(i,j,k,row)       = tolerance*rho_tot(i,j,k)                
                endif
             enddo

             ! compute Bijprime matrix and massfraction W_i = rho_i/rho; molarconc is 
             ! expressed in terms of molmtot,mi,rhotot etc. 
             do row=1, nspecies  
                do column=1, row-1
                   Bijprime(row, column) = rho(i,j,k,row)*molmtot(i,j,k)**2/(mass(row)* &
                                           mass(column)*Dbar(row, column)*rho_tot(i,j,k)**2) 
                   Bijprime(column, row) = rho(i,j,k,column)*molmtot(i,j,k)**2/(mass(row)* &
                                           mass(column)*Dbar(column, row)*rho_tot(i,j,k)**2) 
                enddo
            
                Sum_knoti=0.d0 
                do column=1, nspecies
                   if (column.ne.row) then
                      Sum_knoti = Sum_knoti - rho(i,j,k,column)/(mass(column)*Dbar(row,column))
                   endif
                   Bijprime(row, row) = Sum_knoti*molmtot(i,j,k)**2/(mass(row)*rho_tot(i,j,k)**2)
                enddo
             
                W(row) = rho(i,j,k,row)/rho_tot(i,j,k)
             enddo

             ! compute Binverse*Gamma and other relevant matrices
             call compute_BdagGamma(Bijprime(:,:),BdagGamma(:,:),BdagCapW(:,:),Gama(:,:),W(:),tolerance)

             ! compute Onsager matrix L
             Lonsager = -rho_tot(i,j,k)*Temp*BdagCapW/Press

             ! do the rank conversion 
             call set_Bij(BinvGamma(i,j,k,:), BdagGamma)
 
          end do
       end do
    end do
   
    ! use contained (internal) subroutine to do the copy without doing index algebra:
    contains 
     subroutine set_Bij(BinvGamma_ij, BdagGamma_ij)
        real(kind=dp_t), dimension(nspecies,nspecies), intent(in)  :: BdagGamma_ij
        real(kind=dp_t), dimension(nspecies,nspecies), intent(out) :: BinvGamma_ij  
        BinvGamma_ij = BdagGamma_ij
     end subroutine 

  end subroutine compute_BinvGamma_3d

  subroutine compute_BdagGamma(Bijprime,BdagGamma,BdagCapW,Gama,W,tolerance)
         
    real(kind=dp_t)  :: Bijprime(:,:)
    real(kind=dp_t)  :: BdagGamma(:,:)
    real(kind=dp_t)  :: BdagCapW(:,:)
    real(kind=dp_t)  :: Gama(:,:)
    real(kind=dp_t)  :: W(:)
    real(kind=dp_t)  :: tolerance 
 
    ! local variables
    integer          :: row,column,info
    real(kind=dp_t)  :: alpha1    

    ! vectors and matrices to be used by LAPACK
    real(kind=dp_t), dimension(nspecies,nspecies) :: Bdag, Sdag
    real(kind=dp_t), dimension(nspecies,nspecies) :: U, UT, V, VT
    real(kind=dp_t), dimension(nspecies)          :: S, alpha, work 
    real(kind=dp_t), dimension(nspecies)          :: Bdagw 
    integer,         dimension(nspecies)          :: ipiv

    ! free up the memory  
    Bdag   = 0.d0
    Sdag   = 0.d0
    U      = 0.d0
    UT     = 0.d0
    V      = 0.d0
    VT     = 0.d0
    S      = 0.d0
    alpha  = 0.d0
    work   = 0.d0
    Bdagw  = 0.d0
    alpha1 = 0.d0
        
    ! adjust parameter alpha1 to max val of Bijprime
    alpha1 = maxval(abs(Bijprime))

    ! transform Bijprime to Bij by adding alpha1*matrix of ones(denoted here as Bijprime)
    do row=1, nspecies
       do column=1, nspecies
          Bijprime(row, column) = Bijprime(row, column) + alpha1
       enddo
    enddo
        
    !===============================================================          
    ! select LAPACK inversion type, 1=inverse, 2=pseudo inverse 
    !===============================================================          
    select case(inverse_type) 
           
    case(1)
    !==========================================================
    ! Using Inverse 
    !==========================================================
 
    ! compute A^(-1)*B = c;  
    !call la_gesvx(A=Bijprime, B=Gama, X=BdagGamma)
            
    ! compute Bijprime inverse through LU factorization. 
    ! Bijprime^(-1) is stored in Bijprime.
    call dgetrf(nspecies, nspecies, Bijprime, nspecies, ipiv, info) 
    call dgetri(nspecies, Bijprime, nspecies, ipiv, work, nspecies, info) 
 
    ! populate Bdagger with B^(-1)
    Bdag = Bijprime
     
    case(2) 
    !==========================================================
    ! Using pseudoinverse 
    !==========================================================

    ! SVD decomposition of Bijprime = U * S * VTranspose; note that Bijprime 
    ! is changed. also V=(VT)T, UT = (U)T are needed for pseudoinverse of Bprime.
    call la_gesvd(Bijprime, S, U, VT)
    V = transpose(VT)
    UT = transpose(U)
   
    ! populate diagonal matrix Sdag = 1/S with diagonal=0 below a chosen tolerance
    do row=1, nspecies
       do column=1,nspecies
          Sdag(row,column) = 0.0d0
       enddo
       
       if(S(row).gt.tolerance) then 
          Sdag(row,row) = 1.0d0/S(row)
       else
          Sdag(row,row) = 0.0d0
       endif 
    enddo

    ! compute Bdag = V*Sdag*UT, the pseudoinverse of Bprime & alpha.
    Bdag = matmul(V, matmul(Sdag, UT))

    end select
    !===============================================================          
 
    ! compute alpha 
    alpha = matmul(Bdag, W)

    ! substract alpha from every row element of Bdag to get Bdag*W=0
    do row=1, nspecies
       do column=1, nspecies
          Bdag(row, column) = Bdag(row, column) - alpha(row)
       enddo
    enddo
          
    ! compute B^(-1)*Gamma = Bdag*Gamma. is_ideal_mixture = .true. for ideal mixture
    if(is_ideal_mixture) then
       BdagGamma = Bdag
    else
       BdagGamma = matmul(Bdag, Gama)
    endif 

    ! compute BdagCapW for Onsager matrix L
    do column=1, nspecies
       do row=1, nspecies
          BdagCapW(row, column) = Bdag(row,column)*W(column)
       enddo
    enddo

    ! check Bdag*w = 0 
    if(.false.) then
    Bdagw = matmul(Bdag, W)
       do row=1, nspecies
          do column=1, nspecies
             print*, BdagGamma(row, column)
             print*, Gama(row, column)
          enddo
          print*, Bdagw(row) 
          print*, '' 
       enddo
    endif

  end subroutine compute_BdagGamma

  subroutine compute_coefficient(mla,rho,rho_tot,molarconc,chi,rhoWchiGama,Dbar,Gama,mass,molmtot,the_bc_level)
   
    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: rho(:)
    type(multifab) , intent(in   )  :: rho_tot(:) 
    type(multifab) , intent(in   )  :: molarconc(:) 
    type(multifab) , intent(inout)  :: chi(:) 
    type(multifab) , intent(inout)  :: rhoWchiGama(:) 
    real(kind=dp_t), intent(in   )  :: Dbar(:,:) 
    real(kind=dp_t), intent(in   )  :: Gama(:,:) 
    real(kind=dp_t), intent(in   )  :: mass(:)
    type(multifab) , intent(in   )  :: molmtot(:) 
    type(bc_level) , intent(in   )  :: the_bc_level(:)

    ! local variables
    integer :: lo(rho(1)%dim), hi(rho(1)%dim)
    integer :: n,i,ng,dm,nlevs
 
    ! pointer for rho(nspecies), molarconc(nspecies) 
    real(kind=dp_t), pointer        :: dp(:,:,:,:)   ! for rho    
    real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for rho_tot
    real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for molarconc
    real(kind=dp_t), pointer        :: dp3(:,:,:,:)  ! for chi
    real(kind=dp_t), pointer        :: dp4(:,:,:,:)  ! for molmtot
    real(kind=dp_t), pointer        :: dp5(:,:,:,:)  ! for rhoWchiGama

    dm = mla%dim        ! dimensionality
    ng = rho(1)%ng      ! number of ghost cells 
    nlevs = mla%nlevel  ! number of levels 
 
    ! loop over all boxes 
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp => dataptr(rho(n),       i)
          dp1 => dataptr(rho_tot(n),  i)
          dp2 => dataptr(molarconc(n),i)
          dp3 => dataptr(chi(n),i)
          dp4 => dataptr(molmtot(n),  i)
          dp5 => dataptr(rhoWchiGama(n),  i)
          lo  =  lwb(get_box(rho(n), i))
          hi  =  upb(get_box(rho(n), i))
          
          select case(dm)
          case (2)
             call compute_coefficient_2d(dp(:,:,1,:),dp1(:,:,1,1),dp2(:,:,1,:),&
                  dp3(:,:,1,:),Dbar(:,:),Gama(:,:),mass(:),dp4(:,:,1,1),dp5(:,:,1,:),ng,lo,hi) 
          case (3)
             call compute_coefficient_3d(dp(:,:,:,:),dp1(:,:,:,1),dp2(:,:,:,:),&
                  dp3(:,:,:,:),Dbar(:,:),Gama(:,:),mass(:),dp4(:,:,:,1),dp5(:,:,:,:),ng,lo,hi) 
          end select
       end do
    end do

   end subroutine compute_coefficient
 
subroutine compute_coefficient_2d(rho,rho_tot,molarconc,chi,Dbar,Gama,mass,molmtot,rhoWchiGama,ng,lo,hi)

    integer          :: lo(2), hi(2), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)        ! density; last dimension for species
    real(kind=dp_t)  :: rho_tot(lo(1)-ng:,lo(2)-ng:)      ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,:)  ! molar concentration; 
    real(kind=dp_t)  :: chi(lo(1)-ng:,lo(2)-ng:,:)        ! last dimension for nspecies^2
    real(kind=dp_t)  :: Dbar(:,:)                         ! SM diffusion constants 
    real(kind=dp_t)  :: Gama(:,:)                         ! non-ideality coefficient 
    real(kind=dp_t)  :: mass(:)                           ! species molar mass 
    real(kind=dp_t)  :: molmtot(lo(1)-ng:,lo(2)-ng:)      ! total molar mass 
    real(kind=dp_t)  :: rhoWchiGama(lo(1)-ng:,lo(2)-ng:,:)! last dimension for nspecies^2

    ! local variables
    integer          :: i,j,row,column
    real(kind=dp_t)  :: tolerance,Sum_knoti   ! tolerance set for pinverse 

    ! vectors and matrices to be used by LAPACK 
    real(kind=dp_t), dimension(nspecies,nspecies) :: Lonsager, Lambda
    real(kind=dp_t), dimension(nspecies,nspecies) :: chidag,rhoWchiGamaloc
    real(kind=dp_t), dimension(nspecies)          :: W 

    tolerance = 1e-13

    ! for specific box, now start loops over alloted cells 
    do j=lo(2)-ng,hi(2)+ng
       do i=lo(1)-ng,hi(1)+ng
        
          ! free up memory  
          Lonsager       = 0.d0         
          Lambda         = 0.d0         
          chidag         = 0.d0         
          W              = 0.d0
          rhoWchiGamaloc =0.d0
  
          ! change 0 with tolerance to prevent division by zero in case species
          ! density, molar concentration or total density = 0. 
          do row=1, nspecies
             if(molarconc(i,j,row) .lt. tolerance) then
                molarconc(i,j,row) = tolerance
                rho(i,j,row)       = tolerance*rho_tot(i,j)                
             endif
          enddo

          ! compute Lambda_ij matrix and massfraction W_i = rho_i/rho; molarconc is 
          ! expressed in terms of molmtot,mi,rhotot etc. 
          do row=1, nspecies  
             do column=1, row-1
                Lambda(row, column) = -molarconc(i,j,row)*molarconc(i,j,column)/Dbar(row,column)
                Lambda(column, row) = Lambda(row, column) 
             enddo
             W(row) = rho(i,j,row)/rho_tot(i,j)
          enddo

          ! compute Lambda_ii
          do row=1,nspecies
             Sum_knoti = 0.d0
             do column=1,nspecies
                if(column.ne.row) then
                   Sum_knoti = Sum_knoti - Lambda(row,column)
                endif
                Lambda(row,row) = Sum_knoti
             enddo
          enddo

          ! check Lambda is numerically correct
          if(.false.) then 
          if(i.eq.7 .and. j.eq.14) then
            do row=1, nspecies
               Sum_knoti=0.d0
               do column=1, nspecies 
                  Sum_knoti = Sum_knoti + Lambda(row, column)
                  print*, Lambda(row, column)
               enddo
                  !print*, Sum_knoti 
               print*, ''
            enddo
          endif
          endif

          ! compute chi either selecting inverse/pseudoinverse or iterative methods 
          if(use_lapack) then
             call populate_coefficient(Lambda(:,:),chidag(:,:),rhoWchiGamaloc(:,:),Gama(:,:),W(:),tolerance)
          else
             call Dbar2chi_iterative(nspecies,3,Dbar(:,:),W(:),molarconc(i,j,:),chidag(:,:))
          endif

          ! compute CapW*chi*CapW and Onsager matrix L
          do column=1, nspecies
             do row=1, nspecies
                Lonsager(row, column) = rho_tot(i,j)*rho_tot(i,j)*Temp*W(row)*chidag(row,column)*W(column)/Press
             enddo
          enddo
          !Lonsager = rho_tot(i,j)*rho_tot(i,j)*Temp*CapWchiCapW/Press

          ! compute chi*Gamma here 
          if(is_ideal_mixture) then
            rhoWchiGamaloc = chidag
          else
            rhoWchiGamaloc = matmul(chidag, Gama)
          endif 
          
          ! populate -rho*W*chi*Gama = -rho_i*chi*Gamma
          do row=1, nspecies
             do column=1, nspecies
                rhoWchiGamaloc(row,column) = -rho(i,j,row)*rhoWchiGamaloc(row,column)  
             enddo
          enddo

          ! print to match with previous code
          !if(.false.) then 
          if(i.eq.7 .and. j.eq.14) then
             if(use_lapack) then 
                print*, 'printing rho*W*chi*Gama via inverse/p-inverse'
             else 
                print*, 'printing rho*W*chi*Gama via iterative methods'
             endif 
             do row=1, nspecies
                do column=1, nspecies
                   print*, rhoWchiGamaloc(row, column)
                enddo
                   print*, ''
             enddo
          endif
          !endif

          ! do the rank conversion 
          call set_Bij(chi(i,j,:),         chidag)
          call set_Bij(rhoWchiGama(i,j,:), rhoWchiGamaloc)
 
      end do
    end do
   
    ! use contained (internal) subroutine to do the copy without doing index algebra:
    contains 
     subroutine set_Bij(BinvGamma_ij, BdagGamma_ij)
        real(kind=dp_t), dimension(nspecies,nspecies), intent(in)  :: BdagGamma_ij
        real(kind=dp_t), dimension(nspecies,nspecies), intent(out) :: BinvGamma_ij  
        BinvGamma_ij = BdagGamma_ij
     end subroutine 

  end subroutine compute_coefficient_2d

  subroutine compute_coefficient_3d(rho,rho_tot,molarconc,chi,Dbar,Gama,mass,molmtot,rhoWchiGama,ng,lo,hi)
   
    integer          :: lo(3), hi(3), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)        ! density; last dimension for species
    real(kind=dp_t)  :: rho_tot(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)      ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  ! molar concentration; 
    real(kind=dp_t)  :: chi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)        ! last dimension for nspecies^2
    real(kind=dp_t)  :: Dbar(:,:)                                   ! SM diffusion constants 
    real(kind=dp_t)  :: Gama(:,:)                                   ! non-ideality coefficient 
    real(kind=dp_t)  :: mass(:)                                     ! species molar mass 
    real(kind=dp_t)  :: molmtot(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)      ! total molar mass 
    real(kind=dp_t)  :: rhoWchiGama(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)! last dimension for nspecies^2
    
    ! local variables
    integer          :: i,j,k,row,column
    real(kind=dp_t)  :: tolerance, Sum_knoti   

    ! vectors and matrices to be used by LAPACK 
    real(kind=dp_t), dimension(nspecies,nspecies) :: Lonsager, Lambda
    real(kind=dp_t), dimension(nspecies,nspecies) :: chidag, rhoWchiGamaloc 
    real(kind=dp_t), dimension(nspecies)          :: W 

    tolerance = 1e-13

    ! for specific box, now start loops over alloted cells 
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
        
             ! free up the memory  
             Lonsager       = 0.d0         
             Lambda         = 0.d0         
             chidag         = 0.d0         
             W              = 0.d0
             rhoWchiGamaloc = 0.d0
 
             ! change 0 with tolerance to prevent division by zero in case species
             ! density, molar concentration or total density = 0. 
             do row=1, nspecies
                if(molarconc(i,j,k,row) .lt. tolerance) then
                   molarconc(i,j,k,row) = tolerance
                   rho(i,j,k,row)       = tolerance*rho_tot(i,j,k)                
                endif
             enddo

             ! compute Lambda_ij matrix and massfraction W_i = rho_i/rho; molarconc is 
             ! expressed in terms of molmtot,mi,rhotot etc. 
             do row=1, nspecies  
                do column=1, row-1
                   Lambda(row, column) = -molarconc(i,j,k,row)*molarconc(i,j,k,column)/Dbar(row,column)
                   Lambda(column, row) = Lambda(row, column) 
                enddo
                W(row) = rho(i,j,k,row)/rho_tot(i,j,k)
             enddo
             
             ! compute Lambda_ii
             do row=1,nspecies
                Sum_knoti = 0.d0
                do column=1,nspecies
                   if(column.ne.row) then
                      Sum_knoti = Sum_knoti - Lambda(row,column)
                   endif
                   Lambda(row,row) = Sum_knoti
                enddo
             enddo
             
             ! compute chi  
             call populate_coefficient(Lambda(:,:),chidag(:,:),rhoWchiGamaloc(:,:),Gama(:,:),W(:),tolerance)

             ! compute CapW*chi*CapW and Onsager matrix L
             do column=1, nspecies
                do row=1, nspecies
                   Lonsager(row, column) = rho_tot(i,j,k)*rho_tot(i,j,k)*Temp*W(row)*chidag(row,column)*W(column)/Press
                enddo
             enddo

             ! populate rho*W*chi*Gama (chiGama matrix rows * rho_i)
             do row=1, nspecies
                do column=1, nspecies
                   rhoWchiGamaloc(row,column) = rho(i,j,k,row)*rhoWchiGamaloc(row,column)  
                enddo
             enddo

             ! do the rank conversion 
             call set_Bij(chi(i,j,k,:),         chidag)
             call set_Bij(rhoWchiGama(i,j,k,:), rhoWchiGamaloc)
 
          end do
       end do
    end do
   
    ! use contained (internal) subroutine to do the copy without doing index algebra:
    contains 
     subroutine set_Bij(BinvGamma_ij, BdagGamma_ij)
        real(kind=dp_t), dimension(nspecies,nspecies), intent(in)  :: BdagGamma_ij
        real(kind=dp_t), dimension(nspecies,nspecies), intent(out) :: BinvGamma_ij  
        BinvGamma_ij = BdagGamma_ij
     end subroutine 

  end subroutine compute_coefficient_3d

subroutine populate_coefficient(Lambda,chidag,rhoWchiGamaloc,Gama,W,tolerance)
         
    real(kind=dp_t)  :: Lambda(:,:)
    real(kind=dp_t)  :: chidag(:,:)
    real(kind=dp_t)  :: rhoWchiGamaloc(:,:)
    real(kind=dp_t)  :: Gama(:,:)
    real(kind=dp_t)  :: W(:)
    real(kind=dp_t)  :: tolerance 
 
    ! local variables
    integer          :: row,column,info
    real(kind=dp_t)  :: alpha    

    ! vectors and matrices to be used by LAPACK
    real(kind=dp_t), dimension(nspecies,nspecies) :: Sdag,chilocal
    real(kind=dp_t), dimension(nspecies,nspecies) :: U, UT, V, VT
    real(kind=dp_t), dimension(nspecies)          :: S, work 
    integer,         dimension(nspecies)          :: ipiv

    ! free up the memory  
    Sdag     = 0.d0
    U        = 0.d0
    UT       = 0.d0
    V        = 0.d0
    VT       = 0.d0
    S        = 0.d0
    work     = 0.d0
    alpha    = 0.d0
    chilocal = 0.d0
 
    ! calculate trace(Lambda)
    alpha = 0.d0
    do row=1, nspecies
       alpha = alpha + Lambda(row,row)
    enddo
 
    ! calculate Lambda + alpha*W*WT (Equation 6) 
    do row=1, nspecies
       do column=1, nspecies
          chilocal(row,column) = alpha*W(row)*W(column) + Lambda(row,column)
       enddo
    enddo

    !===============================================================          
    ! select LAPACK inversion type, 1=inverse, 2=pseudo inverse 
    !===============================================================          
    select case(inverse_type) 
           
    case(1)
    !==========================================================
    ! Using Inverse 
    !==========================================================
 
    ! compute chilocal inverse
    call dgetrf(nspecies, nspecies, chilocal, nspecies, ipiv, info) 
    call dgetri(nspecies, chilocal, nspecies, ipiv, work, nspecies, info) 

    ! populate chidag with B^(-1)
    chidag = chilocal   
 
    case(2) 
    !==========================================================
    ! Using pseudoinverse 
    !==========================================================

    ! SVD decomposition of chilocal = U * S * VTranspose; note that chilocal 
    ! is changed. also V=(VT)T, UT = (U)T are needed for pseudoinverse of chilocal.
    call la_gesvd(chilocal, S, U, VT)
    V = transpose(VT)
    UT = transpose(U)
   
    ! populate diagonal matrix Sdag = 1/S with diagonal=0 below a chosen tolerance
    do row=1, nspecies
       do column=1,nspecies
          Sdag(row,column) = 0.0d0
       enddo
       
       if(S(row).gt.tolerance) then 
          Sdag(row,row) = 1.0d0/S(row)
       else
          Sdag(row,row) = 0.0d0
       endif 
    enddo

    ! compute chidag = V*Sdag*UT, the pseudoinverse of chilocal 
    chidag = matmul(V, matmul(Sdag, UT))

    end select
    !===============================================================          
 
    ! compute chi as equation (6) 
    do row=1, nspecies
       do column=1, nspecies
          chidag(row,column) = chidag(row,column) - 1.0d0/alpha
       enddo
    enddo
          
  end subroutine populate_coefficient

end module convert_variables_module
