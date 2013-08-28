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

  subroutine convert_cons_to_prim(mla,rho,rho_tot,molarconc,mass,mtot,the_bc_level)
   
   type(ml_layout), intent(in   )  :: mla
   type(multifab) , intent(inout)  :: rho(:) 
   type(multifab) , intent(inout)  :: rho_tot(:) 
   type(multifab) , intent(inout)  :: molarconc(:) 
   real(kind=dp_t), intent(in   )  :: mass(:)
   real(kind=dp_t), intent(inout)  :: mtot
   type(bc_level) , intent(in   )  :: the_bc_level(:)

   ! local variables
   integer :: lo(rho(1)%dim), hi(rho(1)%dim)
   integer :: n,i,ng,dm,nlevs
 
   ! pointer for rho(nspecies), rho_tot(1), molarconc(nspecies) 
   real(kind=dp_t), pointer        :: dp(:,:,:,:)   ! for rho    
   real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for rho_tot
   real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for molar concentrations

   dm = mla%dim        ! dimensionality
   ng = rho(1)%ng      ! number of ghost cells 
   nlevs = mla%nlevel  ! number of levels 

    ! loop over all boxes 
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp => dataptr(rho(n),i)
          dp1 => dataptr(rho_tot(n),i)
          dp2 => dataptr(molarconc(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          
          select case(dm)
          case (2)
             call compute_molconc_rhotot_2d(dp(:,:,1,:), dp1(:,:,1,1), dp2(:,:,1,:), mass(:), mtot, ng, lo, hi) 
          !case (3)
          !   call init_rho_3d(dp(:,:,:,:), dp1(:,:,:,:), ng, lo, hi, prob_lo, prob_hi, dx(n,:))
          end select
       end do

       ! filling up ghost cells for two adjacent grids at the same level
       ! including periodic domain boundary ghost cells
       call multifab_fill_boundary(rho(n))
       call multifab_fill_boundary(rho_tot(n))
       call multifab_fill_boundary(molarconc(n))

       ! fill non-periodic domain boundary ghost cells 
       !call multifab_physbc(rho(n),      1,1,nspecies,the_bc_level(n))
       !call multifab_physbc(rho_tot(n),  1,1,1       ,the_bc_level(n))
       !call multifab_physbc(molarconc(n),1,1,nspecies,the_bc_level(n))
    end do

  end subroutine convert_cons_to_prim

  subroutine compute_molconc_rhotot_2d(rho, rho_tot, molarconc, mass, mtot, ng, lo, hi)
 
    integer          :: lo(2), hi(2), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)       ! density; last dim for numberof species
    real(kind=dp_t)  :: rho_tot(lo(1)-ng:,lo(2)-ng:)     ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,:) ! molar concentration
    real(kind=dp_t)  :: mass(:)                          ! speciess mass 
    real(kind=dp_t)  :: mtot                             ! total mass 

    ! local variables
    integer          :: i,j,n
    real(kind=dp_t)  :: rho_tot_local

    ! calculate mtot 
    mtot = 0.d0
    do n=1, nspecies  
       mtot = mtot + mass(n)
    enddo
    
    ! for specific box, now start loops over alloted cells    
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          ! calculate total density inside each cell
          rho_tot_local=0.d0 
          do n=1, nspecies  
             rho_tot_local = rho_tot_local + rho(i,j,n)
          enddo         
          rho_tot(i,j) = rho_tot_local
          !print*, rho_tot_local
 
          ! calculate molar concentrations in each cell 
          do n=1, nspecies 
             molarconc(i,j,n) = mtot*rho(i,j,n)/(mass(n)*rho_tot_local)
             !print*, molarconc(i,j,1)
          enddo
        
       enddo
    enddo
 
  end subroutine compute_molconc_rhotot_2d

  subroutine convert_cons_to_BinvGamma(mla,rho,rho_tot,molarconc,BinvGamma,Dbar, & 
             Gama,mass,mtot,the_bc_level)
   
    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in)  :: rho(:)
    type(multifab) , intent(in)  :: rho_tot(:) 
    type(multifab) , intent(in)  :: molarconc(:) 
    type(multifab) , intent(inout)  :: BinvGamma(:) 
    real(kind=dp_t), intent(in   )  :: Dbar(:,:) 
    real(kind=dp_t), intent(in   )  :: Gama(:,:) 
    real(kind=dp_t), intent(in   )  :: mass(:)
    real(kind=dp_t), intent(in   )  :: mtot
    type(bc_level) , intent(in   )  :: the_bc_level(:)

    ! local variables
    integer :: lo(rho(1)%dim), hi(rho(1)%dim)
    integer :: n,i,ng,dm,nlevs
 
    ! pointer for rho(nspecies), molarconc(nspecies) 
    real(kind=dp_t), pointer        :: dp(:,:,:,:)   ! for rho    
    real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for rho_tot
    real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for molarconc
    real(kind=dp_t), pointer        :: dp3(:,:,:,:)  ! for BinvGamma 

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
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          
          select case(dm)
          case (2)
             call compute_BinvGamma_2d(dp(:,:,1,:), dp1(:,:,1,1), dp2(:,:,1,:), dp3(:,:,1,:), & 
                                            Dbar(:,:), Gama(:,:), mass(:), mtot, ng, lo, hi) 
          !case (3)
          !   call init_rho_3d(dp(:,:,:,:), dp1(:,:,:,:), ng, lo, hi, prob_lo, prob_hi, dx(n,:))
          end select
       end do

       ! filling up ghost cells for two adjacent grids at the same level
       ! including periodic domain boundary ghost cells
       !call multifab_fill_boundary(rho(n))
       !call multifab_fill_boundary(rho_tot(n))
       !call multifab_fill_boundary(molarconc(n))
       call multifab_fill_boundary(BinvGamma(n))

       ! fill non-periodic domain boundary ghost cells 
       !call multifab_physbc(rho(n),      1,1,nspecies,the_bc_level(n))
       !call multifab_physbc(rho_tot(n),  1,1,1       ,the_bc_level(n))
       !call multifab_physbc(molarconc(n),1,1,nspecies,the_bc_level(n))
    end do

   end subroutine convert_cons_to_BinvGamma

   subroutine compute_BinvGamma_2d(rho,rho_tot,molarconc,BinvGamma,Dbar,Gama, & 
          mass,mtot,ng,lo,hi)

    integer          :: lo(2), hi(2), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)        ! density; last dimension for species
    real(kind=dp_t)  :: rho_tot(lo(1)-ng:,lo(2)-ng:)      ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,:)  ! molar concentration; 
    real(kind=dp_t)  :: BinvGamma(lo(1)-ng:,lo(2)-ng:,:)  ! B^(-1)*Gamma; last dimension for nspecies^2
    real(kind=dp_t)  :: Dbar(:,:)                         ! SM diffusion constants 
    real(kind=dp_t)  :: Gama(:,:)                         ! non-ideality coefficient 
    real(kind=dp_t)  :: mass(:)                           ! speciess mass 
    real(kind=dp_t)  :: mtot                              ! total mass 

    ! local variables
    integer          :: i,j,k,n,row,column
    real(kind=dp_t)  :: tolerance                         ! tolerance set for pinverse 

    ! vectors and matrices to be used by LAPACK 
    real(kind=dp_t), dimension(nspecies,nspecies) :: Bijprime, Bdag, Sdag
    real(kind=dp_t), dimension(nspecies,nspecies) :: U, UT, V, VT, BdagGamma
    real(kind=dp_t), dimension(nspecies)          :: S, W, alpha 

    ! for specific box, now start loops over alloted cells    
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
         
          ! calculate Bijprime matrix and massfraction W_i = rho_i/rho.
          Bijprime=0.d0
          do row=1, nspecies  
             do column=1, row-1
                Bijprime(row, column) = rho(i,j,row)*mtot**2/(mass(row)* &
                       mass(column)*Dbar(row, column)*rho_tot(i,j)**2) 
                Bijprime(column, row) = rho(i,j,column)*mtot**2/(mass(row)* &
                       mass(column)*Dbar(column, row)*rho_tot(i,j)**2)   
             enddo
             
             do column=1, nspecies
                if (column.ne.row) then
                   Bijprime(row, row) = Bijprime(row, row) - mtot**2/(mass(row)*rho_tot(i,j)**2)* & 
                              (rho(i,j,column)/(mass(column)*Dbar(row,column)))
                endif
             enddo
             W(row) = rho(i,j,row)/rho_tot(i,j)
          enddo

          ! SVD decomposition of Bijprime = U * S * VTranspose; note that Bijprime 
          ! is changed. also V=(VT)T, UT = (U)T are needed for pseudoinverse of
          ! Bprime.
          call la_gesvd(Bijprime, S, U, VT)
          V = transpose(VT)
          UT = transpose(U)

          ! populate diagonal matrix Sdag = 1/S with diagonal=0 below a chosen tolerance
          tolerance = 1e-13
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
          alpha = matmul(Bdag, W)

          ! substract alpha from every row element of Bdag to get Bdag*W=0
          ! tested that psuedoinverse comes correct and Bdag*W=0. 
          do row=1, nspecies
             do column=1, nspecies
                Bdag(row, column) = Bdag(row, column) - alpha(row)
             enddo
          enddo

          ! compute B^(-1)*Gamma which is Bdag*Gamma, result is written to Gamma
          BdagGamma = matmul(Bdag, Gama); 

          ! do the rank conversion 
          call set_Bij(BinvGamma(i,j,:), BdagGamma)
         
       end do
    end do
   
    ! Use contained (internal) subroutine to do the copy without doing index algebra:
    contains 
     subroutine set_Bij(BinvGamma_ij, BdagGamma_ij)
        real(kind=dp_t), dimension(nspecies,nspecies), intent(in)  :: BdagGamma_ij
        real(kind=dp_t), dimension(nspecies,nspecies), intent(out) :: BinvGamma_ij  
        BinvGamma_ij = BdagGamma_ij
     end subroutine 

    end subroutine compute_BinvGamma_2d

end module convert_variables_module
