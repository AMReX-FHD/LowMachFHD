module convert_variables_module

  use multifab_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use ml_layout_module
  use probin_multispecies_module 
  use F95_LAPACK

  implicit none

  private

  public :: convert_cons_to_prim, convert_cons_to_BinvGamma

contains

  subroutine convert_cons_to_prim(mla, rho, rho_tot, molarconc, mass, molmtot, the_bc_level)
   
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

   dm = mla%dim        ! dimensionality
   ng = rho(1)%ng      ! number of ghost cells 
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
             call compute_molconc_rhotot_2d(dp(:,:,1,:), dp1(:,:,1,1), dp2(:,:,1,:), mass(:), & 
                                            dp3(:,:,1,1), ng, lo, hi) 
          case (3)
             stop
             !call init_rho_3d(dp(:,:,:,:), dp1(:,:,:,:), ng, lo, hi, prob_lo, prob_hi, dx(n,:))
          end select
       end do
    end do

  end subroutine convert_cons_to_prim

  subroutine compute_molconc_rhotot_2d(rho, rho_tot, molarconc, mass, molmtot, ng, lo, hi)
 
    integer          :: lo(2), hi(2), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)       ! density; last dim for numberof species
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
          if(i.eq.4 .and. j.eq.4)  write(*,*), "rho1=",rho(i,j,1),"rho2=",rho(i,j,2),"rho1+rho2=",rho(i,j,1)+rho(i,j,2)

          ! calculate mass fraction and total molar mass (1/m=Sum(w_i/m_i)
          Sum_woverm=0.d0
          do n=1, nspecies  
             W(n) = rho(i,j,n)/rho_tot(i,j)
             Sum_woverm = Sum_woverm + W(n)/mass(n)
          enddo
          molmtot(i,j) = 1.0d0/Sum_woverm 
 
          ! calculate molar concentrations in each cell 
          do n=1, nspecies 
             molarconc(i,j,n) = molmtot(i,j)*W(n)/mass(n)
          enddo
        
       enddo
    enddo
 
  end subroutine compute_molconc_rhotot_2d

  subroutine convert_cons_to_BinvGamma(mla, rho, rho_tot, molarconc, BinvGamma, Dbar, & 
             Gama, mass, molmtot, the_bc_level)
   
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
             call compute_BinvGamma_2d(dp(:,:,1,:), dp1(:,:,1,1), dp2(:,:,1,:), dp3(:,:,1,:), & 
                                       Dbar(:,:), Gama(:,:), mass(:), dp4(:,:,1,1), ng, lo, hi) 
          !case (3)
          !   call init_rho_3d(dp(:,:,:,:), dp1(:,:,:,:), ng, lo, hi, prob_lo, prob_hi, dx(n,:))
          end select
       end do
    end do

   end subroutine convert_cons_to_BinvGamma

   subroutine compute_BinvGamma_2d(rho, rho_tot, molarconc, BinvGamma, Dbar, Gama, & 
                                   mass, molmtot, ng, lo, hi)

    integer          :: lo(2), hi(2), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)        ! density; last dimension for species
    real(kind=dp_t)  :: rho_tot(lo(1)-ng:,lo(2)-ng:)      ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,:)  ! molar concentration; 
    real(kind=dp_t)  :: BinvGamma(lo(1)-ng:,lo(2)-ng:,:)  ! B^(-1)*Gamma; last dimension for nspecies^2
    real(kind=dp_t)  :: Dbar(:,:)                         ! SM diffusion constants 
    real(kind=dp_t)  :: Gama(:,:)                         ! non-ideality coefficient 
    real(kind=dp_t)  :: mass(:)                           ! species molar mass 
    real(kind=dp_t)  :: molmtot(lo(1)-ng:,lo(2)-ng:)      ! total molar mass 
    logical          :: is_ideal_mixture

    ! local variables
    integer          :: i,j,row,column,info
    real(kind=dp_t)  :: tolerance, Sum_knoti              ! tolerance set for pinverse 
    real(kind=dp_t)  :: Temp, Pres, alpha1    

    ! vectors and matrices to be used by LAPACK 
    real(kind=dp_t), dimension(nspecies,nspecies) :: Bijprime, Bdag, Sdag, Lonsager
    real(kind=dp_t), dimension(nspecies,nspecies) :: U, UT, V, VT, BdagGamma,BdagCapW 
    real(kind=dp_t), dimension(nspecies)          :: S, W, alpha, Checkmat, work 
    integer, dimension(nspecies)                  :: ipiv

    Temp      = 1.0d0
    Pres      = 1.0d0
    tolerance = 1e-13

    ! for specific box, now start loops over alloted cells    
    do j=lo(2)-ng,hi(2)+ng
       do i=lo(1)-ng,hi(1)+ng
        
          ! free up the memory  
          Bijprime=0.d0
          Bdag = 0.d0
          Sdag = 0.d0
          U = 0.d0
          UT = 0.d0
          V = 0.d0
          VT=0.d0
          BdagGamma = 0.d0
          S = 0.d0
          W=0.d0
          alpha=0.d0
          alpha1=0.d0
          Lonsager=0.d0         
          Checkmat=0.d0
 
          ! change 0 with tolerance to prevent division by zero in case species
          ! density, molar concentration or total density = 0. 
          do row=1, nspecies
             if(molarconc(i,j,row) .lt. tolerance) then
                molarconc(i,j,row) = tolerance
                rho(i,j,row)       = tolerance
             endif
             if(rho_tot(i,j) .lt. tolerance) then
                rho_tot(i,j) = tolerance
             endif
          enddo

          ! calculate Bijprime matrix and massfraction W_i = rho_i/rho, molarconc is 
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

          ! adjust parameter alpha1 to max val of Bijprime
          alpha1 = maxval(abs(Bijprime))

          ! transform Bijprime to Bij by adding alpha1*matrix of ones(denoted here as Bijprime)
          do row=1, nspecies
             do column=1, nspecies
                Bijprime(row, column) = Bijprime(row, column) + alpha1
             enddo
          enddo
        
          ! select LAPACK inversion type, 1=inverse, 2=pseudo inverse 
          select case(inverse_type) 
           
          !%%%%%%%%%%%%%%%%%%% Using Inverse %%%%%%%%%%%%%!
          case(1)
 
             ! compute A^(-1)*B = c;  
             !call la_gesvx(A=Bijprime, B=Gama, X=BdagGamma)
            
             ! compute Bijprime inverse through LU factorization. 
             ! Bijprime^(-1) is stored in Bijprime.
             call dgetrf(nspecies, nspecies, Bijprime, nspecies, ipiv, info) 
             call dgetri(nspecies, Bijprime, nspecies, ipiv, work, nspecies, info) 
  
             ! populate Bdagger with B^(-1)
             Bdag = Bijprime 
             
          !%%%%%%%%%%%%%%%%%%% Using pseudoinverse %%%%%%%%%%%%%!
          case(2) 

             ! SVD decomposition of Bijprime = U * S * VTranspose; note that Bijprime 
             ! is changed. also V=(VT)T, UT = (U)T are needed for pseudoinverse of
             ! Bprime.
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

             ! calculate Bdag = V*Sdag*UT, the pseudoinverse of Bprime & alpha.
             Bdag = matmul(V, matmul(Sdag, UT))

          end select
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
            
          alpha = matmul(Bdag, W)

          ! substract alpha from every row element of Bdag to get Bdag*W=0
          do row=1, nspecies
             do column=1, nspecies
                Bdag(row, column) = Bdag(row, column) - alpha(row)
             enddo
          enddo

          ! calculate Onsager matrix L
          do column=1, nspecies
             do row=1, nspecies
                BdagCapW(row, column) = Bdag(row,column)*W(column)
             enddo
          enddo
          Lonsager = -rho_tot(i,j) * Temp * BdagCapW/Pres
            
          ! compute B^(-1)*Gamma = Bdag*Gamma. is_ideal_mixture = .true. for ideal mixture
          is_ideal_mixture = .false.
          if(is_ideal_mixture) then
             BdagGamma = Bdag
          else
             BdagGamma = matmul(Bdag, Gama)
          endif 

          ! check Bdag*w = 0 
          if(.false.) then
          Checkmat = matmul(Bdag, W)
            if(i.eq.8 .and. j.eq.11) then
              do row=1, nspecies
                do column=1, nspecies
                   print*, BdagGamma(row, column)
                   print*, Lonsager(row, column)
                   print*, Gama(row, column)
                enddo
                print*, Checkmat(row) 
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

end module convert_variables_module
