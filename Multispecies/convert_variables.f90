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

  public :: convert_cons_to_prim, compute_chi, populate_chi, compute_rhoWchiGama

contains

  subroutine convert_cons_to_prim(mla,rho,rho_tot,molarconc,molmtot,molmass,the_bc_level)
   
   type(ml_layout), intent(in   )  :: mla
   type(multifab) , intent(in)     :: rho(:) 
   type(multifab) , intent(inout)  :: rho_tot(:) 
   type(multifab) , intent(inout)  :: molarconc(:) 
   type(multifab) , intent(inout)  :: molmtot(:) 
   real(kind=dp_t), intent(in   )  :: molmass(:)
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
                           dp2(:,:,1,:),molmass(:),dp3(:,:,1,1),ng,lo,hi) 
          case (3)
             call compute_molconc_rhotot_3d(dp(:,:,:,:),dp1(:,:,:,1),&
                           dp2(:,:,:,:),molmass(:),dp3(:,:,:,1),ng,lo,hi) 
          end select
       end do
    end do

  end subroutine convert_cons_to_prim

! Donev: See my comment below regarding routine compute_chi_2d about code duplication between 2d and 3d.
! Same applies here

  subroutine compute_molconc_rhotot_2d(rho,rho_tot,molarconc,molmass,molmtot,ng,lo,hi)
 
    integer          :: lo(2), hi(2), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)       ! density- last dim for #species
    real(kind=dp_t)  :: rho_tot(lo(1)-ng:,lo(2)-ng:)     ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,:) ! molar concentration
    real(kind=dp_t)  :: molmass(:)                       ! species molar mass 
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
             Sum_woverm = Sum_woverm + W(n)/molmass(n)
          enddo
          molmtot(i,j) = 1.0d0/Sum_woverm 
  
          ! calculate molar concentrations in each cell (x_i=m*w_i/m_i) 
          do n=1, nspecies 
             molarconc(i,j,n) = molmtot(i,j)*W(n)/molmass(n)
          enddo

       enddo
    enddo
 
  end subroutine compute_molconc_rhotot_2d

  subroutine compute_molconc_rhotot_3d(rho,rho_tot,molarconc,molmass,molmtot,ng,lo,hi)
 
    integer          :: lo(3), hi(3), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)       ! density- last dim for #species
    real(kind=dp_t)  :: rho_tot(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)     ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:) ! molar concentration
    real(kind=dp_t)  :: molmass(:)                                 ! species molar mass 
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
                Sum_woverm = Sum_woverm + W(n)/molmass(n)
             enddo
             molmtot(i,j,k) = 1.0d0/Sum_woverm 
  
             ! calculate molar concentrations in each cell (x_i=m*w_i/m_i) 
             do n=1, nspecies 
                molarconc(i,j,k,n) = molmtot(i,j,k)*W(n)/molmass(n)
             enddo
        
          enddo
       enddo
    enddo
 
  end subroutine compute_molconc_rhotot_3d

  subroutine compute_chi(mla,rho,rho_tot,molarconc,chi,D_MS,the_bc_level)
   
    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: rho(:)
    type(multifab) , intent(in   )  :: rho_tot(:) 
    type(multifab) , intent(in   )  :: molarconc(:) 
    type(multifab) , intent(inout)  :: chi(:) 
    type(multifab) , intent(in   )  :: D_MS(:) 
    type(bc_level) , intent(in   )  :: the_bc_level(:)

    ! local variables
    integer :: lo(rho(1)%dim), hi(rho(1)%dim)
    integer :: n,i,ng,dm,nlevs
 
    ! pointer for rho(nspecies), molarconc(nspecies) 
    real(kind=dp_t), pointer        :: dp(:,:,:,:)   ! for rho    
    real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for rho_tot
    real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for molarconc
    real(kind=dp_t), pointer        :: dp3(:,:,:,:)  ! for chi
    real(kind=dp_t), pointer        :: dp5(:,:,:,:)  ! for D_MS

    dm = mla%dim        ! dimensionality
    ng = rho(1)%ng      ! number of ghost cells 
    nlevs = mla%nlevel  ! number of levels 
 
    ! loop over all boxes 
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp  => dataptr(rho(n), i)
          dp1 => dataptr(rho_tot(n), i)
          dp2 => dataptr(molarconc(n), i)
          dp3 => dataptr(chi(n), i)
          dp5 => dataptr(D_MS(n), i)
          lo  =  lwb(get_box(rho(n), i))
          hi  =  upb(get_box(rho(n), i))
          
          select case(dm)
          case (2)
             call compute_chi_2d(dp(:,:,1,:),dp1(:,:,1,1),dp2(:,:,1,:),dp3(:,:,1,:),&
                  dp5(:,:,1,:),ng,lo,hi) 
          case (3)
             call compute_chi_3d(dp(:,:,:,:),dp1(:,:,:,1),dp2(:,:,:,:),dp3(:,:,:,:),&
                  dp5(:,:,:,:),ng,lo,hi) 
          end select
       end do
    end do

   end subroutine compute_chi
 
! Donev:
! If you look below you see that a lot of code is duplicated in 2d and 3d.
! This code simply calculatest things in one cell, so it is dimension independent
! You should put the code into a helper routine:
! subroutine compute_chi_local(rho,rho_tot,molarconc,chi,D_MS)
! where now all of the arguments do not have (i,j,k) indices: They are just scalars, vectors of length nspecies, or nspecies^2 matrices
! Then you just call this one routine inside a do-loop in both 2d and 3d
! You should do this consistently in the code for any routine that has duplication in 2d and 3d.
! There should be NEVER code duplication unless it is trivial code, and this one is not trivial.
 
subroutine compute_chi_2d(rho,rho_tot,molarconc,chi,D_MS,ng,lo,hi)

    integer          :: lo(2), hi(2), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)         ! density; last dimension for species
    real(kind=dp_t)  :: rho_tot(lo(1)-ng:,lo(2)-ng:)       ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,:)   ! molar concentration; 
    real(kind=dp_t)  :: chi(lo(1)-ng:,lo(2)-ng:,:)         ! last dimension for nspecies^2
    real(kind=dp_t)  :: D_MS(lo(1)-ng:,lo(2)-ng:,:)        ! MS diff-coeffs 

    ! local variables
    integer          :: i,j,row,column
    real(kind=dp_t)  :: tolerance,Sum_knoti   ! tolerance set for pinverse 

    ! vectors and matrices to be used by LAPACK 
    real(kind=dp_t), dimension(nspecies,nspecies) :: Lonsager,Lambda,chidag,D_MS_loc
    real(kind=dp_t), dimension(nspecies)          :: W 

    tolerance = 1e-13 ! Donev: There is already a variable called fraction_tolerance in module matrix_utilities
    ! Donev: The best thing is to put this variable in probin and make it an input value, and use the same value everywhere

    ! for specific box, now start loops over alloted cells 
    do j=lo(2)-ng,hi(2)+ng
       do i=lo(1)-ng,hi(1)+ng
        
          ! free up memory  
          Lonsager = 0.d0         
          Lambda   = 0.d0         
          chidag   = 0.d0         
          D_MS_loc = 0.d0         
          W        = 0.d0
  
          ! do rank conversion to populate local matrix 
          call set_Xij(D_MS_loc, D_MS(i,j,:))
        
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
                Lambda(row, column) = -molarconc(i,j,row)*molarconc(i,j,column)/D_MS_loc(row,column)
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

          ! compute chi either selecting inverse/pseudoinverse or iterative methods 
          if(use_lapack) then
             call populate_chi(Lambda(:,:),chidag(:,:),W(:),tolerance)
          else
             call Dbar2chi_iterative(nspecies,3,D_MS_loc(:,:),W(:),molarconc(i,j,:),chidag(:,:))
          endif

          ! print chi 
          !if(.false.) then 
          if(i.eq.7 .and. j.eq.14) then
             if(use_lapack) then 
                print*, 'print chi via inverse/p-inverse'
             else 
                print*, 'print chi via iterative methods'
             endif 
             do row=1, nspecies
                do column=1, nspecies
                   print*, chidag(row, column)
                enddo
                   print*, ''
             enddo
          endif
          !endif

          ! do the rank conversion 
          call set_Xij(chi(i,j,:), chidag)
 
      end do
    end do
   
    ! use contained (internal) subroutine to do the copy without doing index algebra:
    contains 
     subroutine set_Xij(Xout_ij, Xin_ij)
        real(kind=dp_t), dimension(nspecies,nspecies), intent(in)  :: Xin_ij
        real(kind=dp_t), dimension(nspecies,nspecies), intent(out) :: Xout_ij 
        Xout_ij = Xin_ij
     end subroutine 

  end subroutine compute_chi_2d

  subroutine compute_chi_3d(rho,rho_tot,molarconc,chi,D_MS,ng,lo,hi)
   
    integer          :: lo(3), hi(3), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)         ! density; last dimension for species
    real(kind=dp_t)  :: rho_tot(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)       ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)   ! molar concentration; 
    real(kind=dp_t)  :: chi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)         ! last dimension for nspecies^2
    real(kind=dp_t)  :: D_MS(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)        ! SM diffusion constants 
    
    ! local variables
    integer          :: i,j,k,row,column
    real(kind=dp_t)  :: tolerance, Sum_knoti   

    ! vectors and matrices to be used by LAPACK 
    real(kind=dp_t), dimension(nspecies,nspecies) :: Lonsager,Lambda,chidag,D_MS_loc
    real(kind=dp_t), dimension(nspecies)          :: W 

    tolerance = 1e-13

    ! for specific box, now start loops over alloted cells 
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
        
             ! free up the memory  
             Lonsager = 0.d0         
             Lambda   = 0.d0         
             chidag   = 0.d0         
             D_MS_loc = 0.d0         
             W        = 0.d0
          
             ! do rank conversion to populate local matrix 
             call set_Xij(D_MS_loc, D_MS(i,j,k,:))
 
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
                   Lambda(row, column) = -molarconc(i,j,k,row)*molarconc(i,j,k,column)/D_MS_loc(row,column)
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
            
             ! compute chi either selecting inverse/pseudoinverse or iterative methods 
             if(use_lapack) then
                call populate_chi(Lambda(:,:),chidag(:,:),W(:),tolerance)
             else
                call Dbar2chi_iterative(nspecies,3,D_MS_loc(:,:),W(:),molarconc(i,j,k,:),chidag(:,:))
             endif
 
             ! do the rank conversion 
             call set_Xij(chi(i,j,k,:), chidag)
 
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

  end subroutine compute_chi_3d

subroutine populate_chi(Lambda,chidag,W,tolerance)
         
    real(kind=dp_t)  :: Lambda(:,:)
    real(kind=dp_t)  :: chidag(:,:)
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
          
  end subroutine populate_chi

  subroutine compute_rhoWchiGama(mla,rho,rho_tot,chi,Gama,rhoWchiGama,the_bc_level)
 
    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: rho(:)
    type(multifab) , intent(in   )  :: rho_tot(:) 
    type(multifab) , intent(in   )  :: chi(:) 
    type(multifab) , intent(in   )  :: Gama(:) 
    type(multifab) , intent(inout)  :: rhoWchiGama(:) 
    type(bc_level) , intent(in   )  :: the_bc_level(:)

    ! local variables
    integer :: lo(rho(1)%dim), hi(rho(1)%dim)
    integer :: n,i,ng,dm,nlevs
 
    ! pointer for rho(nspecies), molarconc(nspecies) 
    real(kind=dp_t), pointer        :: dp(:,:,:,:)   ! for rho    
    real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for rho_tot
    real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for chi
    real(kind=dp_t), pointer        :: dp3(:,:,:,:)  ! for Gama
    real(kind=dp_t), pointer        :: dp4(:,:,:,:)  ! for rhoWchiGama

    dm = mla%dim        ! dimensionality
    ng = rho(1)%ng      ! number of ghost cells 
    nlevs = mla%nlevel  ! number of levels 
 
    ! loop over all boxes 
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp => dataptr(rho(n), i)
          dp1 => dataptr(rho_tot(n), i)
          dp2 => dataptr(chi(n), i)
          dp3 => dataptr(Gama(n), i)
          dp4 => dataptr(rhoWchiGama(n), i)
          lo  =  lwb(get_box(rho(n), i))
          hi  =  upb(get_box(rho(n), i))
          
          select case(dm)
          case (2)
             call compute_rhoWchiGama_2d(dp(:,:,1,:),dp1(:,:,1,1),dp2(:,:,1,:),dp3(:,:,1,:),dp4(:,:,1,:),&
                                ng,lo,hi) 
          case (3)
             call compute_rhoWchiGama_3d(dp(:,:,:,:),dp1(:,:,:,1),dp2(:,:,:,:),dp3(:,:,:,:),dp4(:,:,:,:),&
                                ng,lo,hi) 
          end select
       end do
    end do
 
  end subroutine compute_rhoWchiGama
  
  ! Same code as before: Remove code duplication
 
  subroutine compute_rhoWchiGama_2d(rho,rho_tot,chi,Gama,rhoWchiGama,ng,lo,hi)
  
    integer          :: lo(2), hi(2), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)         ! density; last dimension for species
    real(kind=dp_t)  :: rho_tot(lo(1)-ng:,lo(2)-ng:)       ! total density in each cell 
    real(kind=dp_t)  :: chi(lo(1)-ng:,lo(2)-ng:,:)         ! last dimension for nspecies^2
    real(kind=dp_t)  :: Gama(lo(1)-ng:,lo(2)-ng:,:)        ! non-ideality coefficient 
    real(kind=dp_t)  :: rhoWchiGama(lo(1)-ng:,lo(2)-ng:,:) ! last dimension for nspecies^2

    ! local variables
    integer          :: i,j,row,column
    real(kind=dp_t)  :: tolerance,Sum_knoti   ! tolerance set for pinverse 

    ! vectors and matrices to be used by LAPACK 
    real(kind=dp_t), dimension(nspecies,nspecies) :: Lonsager,Gama_loc
    real(kind=dp_t), dimension(nspecies,nspecies) :: chidag,rhoWchiGamaloc
    real(kind=dp_t), dimension(nspecies)          :: W 

    tolerance = 1e-13

    ! for specific box, now start loops over alloted cells 
    do j=lo(2)-ng,hi(2)+ng
       do i=lo(1)-ng,hi(1)+ng
        
          ! free up memory  
          Lonsager       = 0.d0         
          chidag         = 0.d0         
          Gama_loc       = 0.d0
          rhoWchiGamaloc = 0.d0
          W              = 0.d0
  
          ! do rank conversion to populate local matrix 
          call set_Xij(chidag,   chi(i,j,:))
          call set_Xij(Gama_loc, Gama(i,j,:))
          
          ! compute massfraction W_i = rho_i/rho; 
          do row=1, nspecies  
             W(row) = rho(i,j,row)/rho_tot(i,j)
          enddo

          ! compute CapW*chi*CapW and Onsager matrix L
          ! Donev: Why are you computing Lonsager here?
          ! This should be a separate routine and there should be another multifab that stores Lonsager
          ! We don't need it yet but we will for fluctuations
          do column=1, nspecies
             do row=1, nspecies
                Lonsager(row, column) = rho_tot(i,j)*rho_tot(i,j)*Temp*W(row)*chidag(row,column)*W(column)/Press
             enddo
          enddo

          ! compute chi*Gamma 
          if(is_ideal_mixture) then
            rhoWchiGamaloc = chidag
          else
            rhoWchiGamaloc = matmul(chidag, Gama_loc)
          endif 
          
          ! populate -rho*W*chi*Gama = -rho_i*chi*Gamma
          do row=1, nspecies
             do column=1, nspecies
                rhoWchiGamaloc(row,column) = -rho(i,j,row)*rhoWchiGamaloc(row,column)  
             enddo
          enddo

          if(.false.) then
          if(i.eq.7 .and. j.eq.14) then
            do row=1, nspecies
               do column=1, nspecies
                  print*, rhoWchiGamaloc(row,column) 
               enddo
               print*, '' 
            enddo
          endif
          endif

          ! do the rank conversion 
          call set_Xij(rhoWchiGama(i,j,:), rhoWchiGamaloc)
 
      end do
    end do
   
    ! use contained (internal) subroutine to do the copy without doing index algebra:
    contains 
     subroutine set_Xij(Xout_ij, Xin_ij)
        real(kind=dp_t), dimension(nspecies,nspecies), intent(in)  :: Xin_ij
        real(kind=dp_t), dimension(nspecies,nspecies), intent(out) :: Xout_ij 
        Xout_ij = Xin_ij
     end subroutine
 
  end subroutine compute_rhoWchiGama_2d

  subroutine compute_rhoWchiGama_3d(rho,rho_tot,chi,Gama,rhoWchiGama,ng,lo,hi)

    integer          :: lo(3), hi(3), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)         ! density; last dimension for species
    real(kind=dp_t)  :: rho_tot(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)       ! total density in each cell 
    real(kind=dp_t)  :: chi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)         ! last dimension for nspecies^2
    real(kind=dp_t)  :: Gama(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)        ! non-ideality coefficient 
    real(kind=dp_t)  :: rhoWchiGama(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:) ! last dimension for nspecies^2
    
    ! local variables
    integer          :: i,j,k,row,column
    real(kind=dp_t)  :: tolerance,Sum_knoti   ! tolerance set for pinverse 

    ! vectors and matrices to be used by LAPACK 
    real(kind=dp_t), dimension(nspecies,nspecies) :: Lonsager,Gama_loc
    real(kind=dp_t), dimension(nspecies,nspecies) :: chidag,rhoWchiGamaloc
    real(kind=dp_t), dimension(nspecies)          :: W 

    tolerance = 1e-13

    ! for specific box, now start loops over alloted cells 
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
        
             ! free up memory  
             Lonsager       = 0.d0         
             chidag         = 0.d0         
             Gama_loc       = 0.d0
             rhoWchiGamaloc = 0.d0
             W              = 0.d0
  
             ! do rank conversion to populate local matrix 
             call set_Xij(chidag,   chi(i,j,k,:))
             call set_Xij(Gama_loc, Gama(i,j,k,:))
          
             ! compute massfraction W_i = rho_i/rho; 
             do row=1, nspecies  
                W(row) = rho(i,j,k,row)/rho_tot(i,j,k)
             enddo

             ! compute CapW*chi*CapW and Onsager matrix L
             do column=1, nspecies
                do row=1, nspecies
                   Lonsager(row, column) = rho_tot(i,j,k)*rho_tot(i,j,k)*Temp*W(row)*chidag(row,column)*W(column)/Press
                enddo
             enddo

             ! compute chi*Gamma 
             if(is_ideal_mixture) then
                rhoWchiGamaloc = chidag
             else
                rhoWchiGamaloc = matmul(chidag, Gama_loc)
             endif 
          
             ! populate -rho*W*chi*Gama = -rho_i*chi*Gamma
             do row=1, nspecies
                do column=1, nspecies
                   rhoWchiGamaloc(row,column) = -rho(i,j,k,row)*rhoWchiGamaloc(row,column) 
                enddo
             enddo

             ! do the rank conversion 
             call set_Xij(rhoWchiGama(i,j,k,:), rhoWchiGamaloc)
 
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

  end subroutine compute_rhoWchiGama_3d

end module convert_variables_module
