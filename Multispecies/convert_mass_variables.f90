module convert_mass_variables_module

  use multifab_module
  use define_bc_module
  use ml_layout_module
  use probin_common_module
  use probin_multispecies_module
  use matrix_utilities 
  use F95_LAPACK

  implicit none

  private

  public :: convert_cons_to_prim,compute_chi,compute_rhoWchi,compute_Lonsager, &
            compute_molconc_rhotot_local,compute_chi_local,compute_Lonsager_local, &
            correct_rho_with_drho,correct_rho_with_drho_local

contains
  
  subroutine correct_rho_with_drho(mla,rho,drho,the_bc_level)

   type(ml_layout), intent(in   )  :: mla
   type(multifab) , intent(inout)  :: rho(:) 
   type(multifab) , intent(inout)  :: drho(:) 
   type(bc_level) , intent(in   )  :: the_bc_level(:)

   ! local variables
   integer :: lo(rho(1)%dim), hi(rho(1)%dim)
   integer :: n,i,ng,dm,nlevs
 
   ! pointer for rho(nspecies), rho_tot(1), molarconc(nspecies) 
   real(kind=dp_t), pointer        :: dp(:,:,:,:)   ! for rho    
   real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for drho

   dm    = mla%dim     ! dimensionality
   ng    = rho(1)%ng   ! number of ghost cells 
   nlevs = mla%nlevel  ! number of levels 
 
    ! loop over all boxes 
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp => dataptr(rho(n),i)
          dp1 => dataptr(drho(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          
          select case(dm)
          case (2)
             call correct_rho_with_drho_2d(dp(:,:,1,:),dp1(:,:,1,:),ng,lo,hi) 
          case (3)
             call correct_rho_with_drho_3d(dp(:,:,:,:),dp1(:,:,:,:),ng,lo,hi) 
          end select
       end do
    end do

  end subroutine correct_rho_with_drho

  subroutine correct_rho_with_drho_2d(rho,drho,ng,lo,hi)
 
    integer          :: lo(2), hi(2), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)       ! density- last dim for #species
    real(kind=dp_t)  :: drho(lo(1)-ng:,lo(2)-ng:,:)      ! total density in each cell 
        
    ! local variables
    integer          :: i,j
    
    ! for specific box, now start loops over alloted cells    
    do j=lo(2)-ng, hi(2)+ng
       do i=lo(1)-ng, hi(1)+ng
         
         call correct_rho_with_drho_local(rho(i,j,:),drho(i,j,:))

       end do
    end do
 
  end subroutine correct_rho_with_drho_2d

  subroutine correct_rho_with_drho_3d(rho,drho,ng,lo,hi)
 
    integer          :: lo(3), hi(3), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)       ! density- last dim for #species
    real(kind=dp_t)  :: drho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)     ! total density in each cell 
    
    ! local variables
    integer          :: i,j,k
    
    ! for specific box, now start loops over alloted cells    
    do k=lo(3)-ng, hi(3)+ng
       do j=lo(2)-ng, hi(2)+ng
          do i=lo(1)-ng, hi(1)+ng

             call correct_rho_with_drho_local(rho(i,j,k,:),drho(i,j,k,:))

          end do
       end do
    end do
 
  end subroutine correct_rho_with_drho_3d

  subroutine correct_rho_with_drho_local(rho,drho)
 
    real(kind=dp_t), intent(inout) :: rho(nspecies)    ! density- last dim for #species
    real(kind=dp_t), intent(out)   :: drho(nspecies)   ! total density in each cell 
    
    ! local variables
    integer          :: row
    real(kind=dp_t)  :: rho_tot_local

    rho_tot_local = sum(rho)  ! total rho in the cell
   
    ! change 0 with tolerance to prevent division by zero in case species
    ! density, molar concentration or total density = 0.
    do row=1, nspecies
       if(rho(row) .lt. fraction_tolerance*rho_tot_local) then
           drho(row) = fraction_tolerance*rho_tot_local
       else
           drho(row) = 0.0d0
       end if
    end do
 
    ! modify rho  
    rho = rho + drho 

  end subroutine correct_rho_with_drho_local 

  subroutine convert_cons_to_prim(mla,rho,rho_tot,molarconc,molmtot,molmass,the_bc_level)
   
   type(ml_layout), intent(in   )  :: mla
   type(multifab) , intent(in   )  :: rho(:) 
   type(multifab) , intent(inout)  :: rho_tot(:) 
   type(multifab) , intent(inout)  :: molarconc(:) 
   type(multifab) , intent(inout)  :: molmtot(:) 
   real(kind=dp_t), intent(in   )  :: molmass(nspecies)
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
                           dp2(:,:,1,:),molmass,dp3(:,:,1,1),ng,lo,hi) 
          case (3)
             call compute_molconc_rhotot_3d(dp(:,:,:,:),dp1(:,:,:,1),&
                           dp2(:,:,:,:),molmass,dp3(:,:,:,1),ng,lo,hi) 
          end select
       end do
    end do

  end subroutine convert_cons_to_prim

  subroutine compute_molconc_rhotot_2d(rho,rho_tot,molarconc,molmass,molmtot,ng,lo,hi)
 
    integer          :: lo(2), hi(2), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)       ! density- last dim for #species
    real(kind=dp_t)  :: rho_tot(lo(1)-ng:,lo(2)-ng:)     ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,:) ! molar concentration
    real(kind=dp_t)  :: molmass(nspecies)                ! species molar mass 
    real(kind=dp_t)  :: molmtot(lo(1)-ng:,lo(2)-ng:)     ! total molar mass 
        
    ! local variables
    integer          :: i,j
    
    ! for specific box, now start loops over alloted cells    
    do j=lo(2)-ng, hi(2)+ng
       do i=lo(1)-ng, hi(1)+ng
         
         call compute_molconc_rhotot_local(rho(i,j,:),rho_tot(i,j),&
              molarconc(i,j,:),molmass,molmtot(i,j))

       end do
    end do
 
  end subroutine compute_molconc_rhotot_2d

  subroutine compute_molconc_rhotot_3d(rho,rho_tot,molarconc,molmass,molmtot,ng,lo,hi)
 
    integer          :: lo(3), hi(3), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)       ! density- last dim for #species
    real(kind=dp_t)  :: rho_tot(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)     ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:) ! molar concentration
    real(kind=dp_t)  :: molmass(nspecies)                                 ! species molar mass 
    real(kind=dp_t)  :: molmtot(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)     ! total molar mass 
    
    ! local variables
    integer          :: i,j,k
    
    ! for specific box, now start loops over alloted cells    
    do k=lo(3)-ng, hi(3)+ng
       do j=lo(2)-ng, hi(2)+ng
          do i=lo(1)-ng, hi(1)+ng

             call compute_molconc_rhotot_local(rho(i,j,k,:),rho_tot(i,j,k),&
                  molarconc(i,j,k,:),molmass,molmtot(i,j,k))

          end do
       end do
    end do
 
  end subroutine compute_molconc_rhotot_3d

  subroutine compute_molconc_rhotot_local(rho,rho_tot,molarconc,molmass,molmtot)
 
    real(kind=dp_t), intent(in)   :: rho(nspecies)       ! density- last dim for #species
    real(kind=dp_t), intent(out)  :: rho_tot             ! total density in each cell 
    real(kind=dp_t), intent(out)  :: molarconc(nspecies) ! molar concentration
    real(kind=dp_t), intent(in)   :: molmass(nspecies)   ! species molar mass 
    real(kind=dp_t), intent(out)  :: molmtot             ! total molar mass 
    
    ! local variables
    integer          :: n
    real(kind=dp_t), dimension(nspecies) :: W            ! mass fraction w_i = rho_i/rho 
    real(kind=dp_t)  :: Sum_woverm, rho_tot_local

    ! calculate total density inside each cell
    rho_tot_local=0.d0 
    do n=1, nspecies  
       rho_tot_local = rho_tot_local + rho(n)
    end do         
    rho_tot = rho_tot_local

    ! calculate mass fraction and total molar mass (1/m=Sum(w_i/m_i))
    Sum_woverm=0.d0
    do n=1, nspecies  
       W(n) = rho(n)/rho_tot
       Sum_woverm = Sum_woverm + W(n)/molmass(n)
    end do
    molmtot = 1.0d0/Sum_woverm 

    ! calculate molar concentrations in each cell (x_i=m*w_i/m_i) 
    do n=1, nspecies 
       molarconc(n) = molmtot*W(n)/molmass(n)
    end do
    
  end subroutine compute_molconc_rhotot_local 

  subroutine compute_chi(mla,rho,rho_tot,molarconc,molmass,chi,D_MS,Temp,zeta_by_Temp,the_bc_level)
   
    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: rho(:)
    type(multifab) , intent(in   )  :: rho_tot(:) 
    type(multifab) , intent(in   )  :: molarconc(:)
    real(kind=dp_t), intent(in   )  :: molmass(nspecies)   ! species molar mass  
    type(multifab) , intent(inout)  :: chi(:) 
    type(multifab) , intent(in   )  :: D_MS(:) 
    type(multifab) , intent(in   )  :: Temp(:) 
    type(multifab) , intent(inout)  :: zeta_by_Temp(:) 
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
    real(kind=dp_t), pointer        :: dp6(:,:,:,:)  ! for Temp
    real(kind=dp_t), pointer        :: dp7(:,:,:,:)  ! for zeta_by_Temp

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
          dp6 => dataptr(Temp(n), i)
          dp7 => dataptr(zeta_by_Temp(n), i)
          lo  =  lwb(get_box(rho(n), i))
          hi  =  upb(get_box(rho(n), i))
          
          select case(dm)
          case (2)
             call compute_chi_2d(dp(:,:,1,:),dp1(:,:,1,1),dp2(:,:,1,:),molmass,dp3(:,:,1,:),&
                  dp5(:,:,1,:),dp6(:,:,1,1),dp7(:,:,1,:),ng,lo,hi) 
          case (3)
             call compute_chi_3d(dp(:,:,:,:),dp1(:,:,:,1),dp2(:,:,:,:),molmass,dp3(:,:,:,:),&
                  dp5(:,:,:,:),dp6(:,:,:,1),dp7(:,:,:,:),ng,lo,hi) 
          end select
       end do
    end do

  end subroutine compute_chi
 
  subroutine compute_chi_2d(rho,rho_tot,molarconc,molmass,chi,D_MS,Temp,zeta_by_Temp,ng,lo,hi)

    integer          :: lo(2), hi(2), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)         ! density; last dimension for species
    real(kind=dp_t)  :: rho_tot(lo(1)-ng:,lo(2)-ng:)       ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,:)   ! molar concentration 
    real(kind=dp_t)  :: molmass(nspecies)                  ! species molar mass 
    real(kind=dp_t)  :: chi(lo(1)-ng:,lo(2)-ng:,:)         ! last dimension for nspecies^2
    real(kind=dp_t)  :: D_MS(lo(1)-ng:,lo(2)-ng:,:)        ! MS diff-coeffs 
    real(kind=dp_t)  :: Temp(lo(1)-ng:,lo(2)-ng:)          ! Temperature 
    real(kind=dp_t)  :: zeta_by_Temp(lo(1)-ng:,lo(2)-ng:,:)

    ! local variables
    integer          :: i,j,row,column
    real(kind=dp_t), dimension(nspecies,nspecies) :: chilocal

    ! for specific box, now start loops over alloted cells 
    do j=lo(2)-ng,hi(2)+ng
       do i=lo(1)-ng,hi(1)+ng
    
          call compute_chi_local(rho(i,j,:),rho_tot(i,j),molarconc(i,j,:),molmass,&
                                 chi(i,j,:),D_MS(i,j,:),Temp(i,j),zeta_by_Temp(i,j,:))

          ! print chi for one cell 
          if(.false.) then 
          if(i.eq.7 .and. j.eq.14) then
            if(use_lapack) then 
              print*, 'print chi via inverse/p-inverse'
            else 
              print*, 'print chi via iterative methods'
            end if
            call set_Xij(chilocal, chi(i,j,:)) 
            do row=1, nspecies
               do column=1, nspecies
                  print*, chilocal(row, column)
               end do
               print*, ''
            end do
          end if
          end if

       end do
    end do

  end subroutine compute_chi_2d

  subroutine compute_chi_3d(rho,rho_tot,molarconc,molmass,chi,D_MS,Temp,zeta_by_Temp,ng,lo,hi)
   
    integer          :: lo(3), hi(3), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)         ! density; last dimension for species
    real(kind=dp_t)  :: rho_tot(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)       ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)   ! molar concentration; 
    real(kind=dp_t)  :: molmass(nspecies)                            ! species molar mass 
    real(kind=dp_t)  :: chi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)         ! last dimension for nspecies^2
    real(kind=dp_t)  :: D_MS(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)        ! SM diffusion constants 
    real(kind=dp_t)  :: Temp(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)          ! Temperature 
    real(kind=dp_t)  :: zeta_by_Temp(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    
    ! local variables
    integer          :: i,j,k

    ! for specific box, now start loops over alloted cells 
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
       
             call compute_chi_local(rho(i,j,k,:),rho_tot(i,j,k),molarconc(i,j,k,:),molmass,&
                                    chi(i,j,k,:),D_MS(i,j,k,:),Temp(i,j,k),zeta_by_Temp(i,j,k,:))

          end do
       end do
    end do
   
  end subroutine compute_chi_3d

  subroutine compute_chi_local(rho,rho_tot,molarconc,molmass,chi,D_MS,Temp,zeta_by_Temp)
    
    real(kind=dp_t), intent(inout)  :: rho(nspecies)         
    real(kind=dp_t), intent(in)     :: rho_tot               
    real(kind=dp_t), intent(inout)  :: molarconc(nspecies) 
    real(kind=dp_t), intent(in)     :: molmass(nspecies)   ! species molar mass   
    real(kind=dp_t), intent(out)    :: chi(nspecies,nspecies)   
    real(kind=dp_t), intent(in)     :: D_MS(nspecies,nspecies) 
    real(kind=dp_t), intent(in)     :: Temp
    real(kind=dp_t), intent(inout)  :: zeta_by_Temp(nspecies)

    ! local variables
    integer                         :: row,column
    real(kind=dp_t)                 :: Sum_knoti   

    ! vectors and matrices to be used by LAPACK 
    real(kind=dp_t), dimension(nspecies,nspecies) :: Lambda
    real(kind=dp_t), dimension(nspecies)          :: W
      
    ! free up memory  
    Lambda        = 0.d0         
    W             = 0.d0

    ! compute Lambda_ij matrix and massfraction W_i = rho_i/rho; molarconc is 
    ! expressed in terms of molmtot,mi,rhotot etc. 
    do row=1, nspecies  
       do column=1, row-1
          Lambda(row, column) = -molarconc(row)*molarconc(column)/D_MS(row,column)
          Lambda(column, row) = Lambda(row, column) 
       end do
       W(row) = rho(row)/rho_tot
    end do

    ! compute Lambda_ii
    do row=1, nspecies
       Sum_knoti = 0.d0
       do column=1, nspecies
          if(column.ne.row) then
             Sum_knoti = Sum_knoti - Lambda(row,column)
          end if
          Lambda(row,row) = Sum_knoti
       end do
    end do

    ! compute zeta_by_Temp for thermodiffusion
    if(is_nonisothermal) then
       do row=1, nspecies
          Sum_knoti = 0.d0
          do column=1, nspecies
             if(column.ne.row) then
                Sum_knoti = Sum_knoti + Lambda(row,column)*(DT_in(row)-DT_in(column))
             end if
             zeta_by_Temp(row) = Sum_knoti/Temp
          end do
       end do
    end if
    
    ! compute chi either selecting inverse/pseudoinverse or iterative methods 
    if(use_lapack) then
       call compute_chi_lapack(Lambda(:,:),chi(:,:),W(:))
    else
       call Dbar2chi_iterative(nspecies,100,D_MS(:,:),molmass(:),molarconc(:),chi(:,:)) 
    end if

  end subroutine compute_chi_local

  subroutine compute_chi_lapack(Lambda,chi,W)
    
    real(kind=dp_t)  :: Lambda(nspecies,nspecies)
    real(kind=dp_t)  :: chi(nspecies,nspecies)
    real(kind=dp_t)  :: W(nspecies)
 
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
    end do
 
    ! calculate Lambda + alpha*W*WT (Equation 6) 
    do row=1, nspecies
       do column=1, nspecies
          chilocal(row,column) = alpha*W(row)*W(column) + Lambda(row,column)
       end do
    end do

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

    ! populate chi with B^(-1)
    chi = chilocal   
 
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
       end do
       
       if(S(row).gt.fraction_tolerance) then 
          Sdag(row,row) = 1.0d0/S(row)
       else
          Sdag(row,row) = 0.0d0
       end if 
    end do

    ! compute chi = V*Sdag*UT, the pseudoinverse of chilocal 
    chi = matmul(V, matmul(Sdag, UT))

    end select
    !===============================================================          
 
    ! compute chi as equation (6) 
    do row=1, nspecies
       do column=1, nspecies
          chi(row,column) = chi(row,column) - 1.0d0/alpha
       end do
    end do
          
  end subroutine compute_chi_lapack

  subroutine compute_Lonsager(mla,rho,rho_tot,molarconc,molmass,molmtot,chi,Gama,Lonsager,the_bc_level)
 
    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: rho(:)
    type(multifab) , intent(in   )  :: rho_tot(:) 
    type(multifab) , intent(in   )  :: molarconc(:) 
    type(multifab) , intent(in   )  :: molmtot(:) 
    real(kind=dp_t), intent(in   )  :: molmass(nspecies)
    type(multifab) , intent(in   )  :: chi(:) 
    type(multifab) , intent(in   )  :: Gama(:) 
    type(multifab) , intent(inout)  :: Lonsager(:)
    type(bc_level) , intent(in   )  :: the_bc_level(:)

    ! local variables
    integer :: lo(rho(1)%dim), hi(rho(1)%dim)
    integer :: n,i,ng,dm,nlevs
 
    ! pointer for rho(nspecies), molarconc(nspecies) 
    real(kind=dp_t), pointer        :: dp(:,:,:,:)   ! for rho    
    real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for rho_tot
    real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for molarconc
    real(kind=dp_t), pointer        :: dp3(:,:,:,:)  ! for molmtot
    real(kind=dp_t), pointer        :: dp4(:,:,:,:)  ! for chi
    real(kind=dp_t), pointer        :: dp5(:,:,:,:)  ! for Gama
    real(kind=dp_t), pointer        :: dp6(:,:,:,:)  ! for Lonsager 

    dm = mla%dim        ! dimensionality
    ng = rho(1)%ng      ! number of ghost cells 
    nlevs = mla%nlevel  ! number of levels 
 
    ! loop over all boxes 
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp => dataptr(rho(n), i)
          dp1 => dataptr(rho_tot(n), i)
          dp2 => dataptr(molarconc(n),i)
          dp3 => dataptr(molmtot(n),i)
          dp4 => dataptr(chi(n), i)
          dp5 => dataptr(Gama(n), i)
          dp6 => dataptr(Lonsager(n), i)
          lo  =  lwb(get_box(rho(n), i))
          hi  =  upb(get_box(rho(n), i))
          
          select case(dm)
          case (2)
             call compute_Lonsager_2d(dp(:,:,1,:),dp1(:,:,1,1),dp2(:,:,1,:),molmass,&
                           dp3(:,:,1,1),dp4(:,:,1,:),dp5(:,:,1,:),dp6(:,:,1,:),ng,lo,hi) 
          case (3)
             call compute_Lonsager_3d(dp(:,:,:,:),dp1(:,:,:,1),dp2(:,:,:,:),molmass,&
                           dp3(:,:,:,1),dp4(:,:,:,:),dp5(:,:,:,:),dp6(:,:,:,:),ng,lo,hi) 
          end select
       end do
    end do
 
  end subroutine compute_Lonsager
  
  subroutine compute_Lonsager_2d(rho,rho_tot,molarconc,molmass,molmtot,chi,Gama,Lonsager,ng,lo,hi)
  
    integer          :: lo(2), hi(2), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)         ! density; last dimension for species
    real(kind=dp_t)  :: rho_tot(lo(1)-ng:,lo(2)-ng:)       ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,:)   ! molar concentration
    real(kind=dp_t)  :: molmass(nspecies)                  ! species molar mass 
    real(kind=dp_t)  :: molmtot(lo(1)-ng:,lo(2)-ng:)       ! total molar mass 
    real(kind=dp_t)  :: chi(lo(1)-ng:,lo(2)-ng:,:)         ! last dimension for nspecies^2
    real(kind=dp_t)  :: Gama(lo(1)-ng:,lo(2)-ng:,:)        ! non-ideality coefficient 
    real(kind=dp_t)  :: Lonsager(lo(1)-ng:,lo(2)-ng:,:) ! last dimension for nspecies^2

    ! local variables
    integer          :: i,j,row,column
    real(kind=dp_t), dimension(nspecies,nspecies) :: Lonsager_local 
 
    ! for specific box, now start loops over alloted cells 
    do j=lo(2)-ng,hi(2)+ng
       do i=lo(1)-ng,hi(1)+ng
        
          call compute_Lonsager_local(rho(i,j,:),rho_tot(i,j),molarconc(i,j,:),&
                          molmass,molmtot(i,j),chi(i,j,:),Gama(i,j,:),Lonsager(i,j,:))

          if(.false.) then
          if(i.eq.7 .and. j.eq.14) then
          call set_Xij(Lonsager_local, Lonsager(i,j,:)) 
            write(*,*), "Lonsager"
            do row=1, nspecies
               write(*,*), Lonsager_local(row,:) 
               !write(*,*), sum(Lonsager_local(row,:))
            end do
            do row=1, nspecies
               !write(*,*), sum(Lonsager_local(:,row))
            end do
          end if
          end if

       end do
    end do

  end subroutine compute_Lonsager_2d

  subroutine compute_Lonsager_3d(rho,rho_tot,molarconc,molmass,molmtot,chi,Gama,Lonsager,ng,lo,hi)

    integer          :: lo(3), hi(3), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)         ! density; last dimension for species
    real(kind=dp_t)  :: rho_tot(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)       ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:) ! molar concentration
    real(kind=dp_t)  :: molmass(nspecies)                                 ! species molar mass 
    real(kind=dp_t)  :: molmtot(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)     ! total molar mass 
    real(kind=dp_t)  :: chi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)         ! last dimension for nspecies^2
    real(kind=dp_t)  :: Gama(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)        ! non-ideality coefficient 
    real(kind=dp_t)  :: Lonsager(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:) ! last dimension for nspecies^2
    
    ! local variables
    integer          :: i,j,k

    ! for specific box, now start loops over alloted cells 
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
       
             call compute_Lonsager_local(rho(i,j,k,:),rho_tot(i,j,k),molarconc(i,j,k,:),&
                          molmass,molmtot(i,j,k),chi(i,j,k,:),Gama(i,j,k,:),Lonsager(i,j,k,:))
              
         end do
      end do
    end do
   
  end subroutine compute_Lonsager_3d

subroutine compute_Lonsager_local(rho,rho_tot,molarconc,molmass,molmtot,chi,Gama,Lonsager)
   
    real(kind=dp_t), intent(in)   :: rho(nspecies)            
    real(kind=dp_t), intent(in)   :: rho_tot                  
    real(kind=dp_t), intent(in)   :: molarconc(nspecies)      ! molar concentration
    real(kind=dp_t), intent(in)   :: molmass(nspecies)        ! species molar mass 
    real(kind=dp_t), intent(in)   :: molmtot                  ! total molar mass 
    real(kind=dp_t), intent(in)   :: chi(nspecies,nspecies)   ! rank conversion done 
    real(kind=dp_t), intent(in)   :: Gama(nspecies,nspecies)        
    real(kind=dp_t), intent(out)  :: Lonsager(nspecies,nspecies) 

    ! local variables
    integer                              :: row,column,info
    real(kind=dp_t), dimension(nspecies) :: W 
    real(kind=dp_t)                      :: rcond 
  
    ! compute massfraction W_i = rho_i/rho; 
    do row=1, nspecies  
       W(row) = rho(row)/rho_tot
    end do

    ! compute Onsager matrix L
    do column=1, nspecies
       do row=1, nspecies
          ! Donev: This will need to be modified when we go to variable temperature
          !Lonsager(row, column) = rho_tot*rho_tot*Temp*W(row)*chi(row,column)*W(column)/Press
          Lonsager(row, column) = molmtot*rho_tot*W(row)*chi(row,column)*W(column)/k_B
       end do
    end do

    ! compute cell-centered Cholesky factor of Lonsager
    if(use_lapack) then
       
       call dpotrf_f95(Lonsager,'L', rcond, 'I', info)
    
       ! remove all upper-triangular entries and NXN entry that lapack doesn't set to zero 
       do row=1, nspecies
          do column=row+1, nspecies
             Lonsager(row, column) = 0.0d0          
          end do
       end do    
       Lonsager(nspecies, nspecies) = 0.0d0          
    
    else
       call choldc(Lonsager,nspecies)   
    end if

  end subroutine compute_Lonsager_local

  subroutine compute_rhoWchi(mla,rho,rho_tot,molarconc,molmass,molmtot,chi,rhoWchi,the_bc_level)
 
    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: rho(:)
    type(multifab) , intent(in   )  :: rho_tot(:) 
    type(multifab) , intent(in   )  :: molarconc(:) 
    type(multifab) , intent(in   )  :: molmtot(:) 
    real(kind=dp_t), intent(in   )  :: molmass(nspecies)
    type(multifab) , intent(in   )  :: chi(:) 
    type(multifab) , intent(inout)  :: rhoWchi(:) 
    type(bc_level) , intent(in   )  :: the_bc_level(:)

    ! local variables
    integer :: lo(rho(1)%dim), hi(rho(1)%dim)
    integer :: n,i,ng,dm,nlevs
 
    ! pointer for rho(nspecies), molarconc(nspecies) 
    real(kind=dp_t), pointer        :: dp(:,:,:,:)   ! for rho    
    real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for rho_tot
    real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for molarconc
    real(kind=dp_t), pointer        :: dp3(:,:,:,:)  ! for molmtot
    real(kind=dp_t), pointer        :: dp4(:,:,:,:)  ! for chi
    real(kind=dp_t), pointer        :: dp5(:,:,:,:)  ! for rhoWchi

    dm = mla%dim        ! dimensionality
    ng = rho(1)%ng      ! number of ghost cells 
    nlevs = mla%nlevel  ! number of levels 

 
    ! loop over all boxes 
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp => dataptr(rho(n), i)
          dp1 => dataptr(rho_tot(n), i)
          dp2 => dataptr(molarconc(n),i)
          dp3 => dataptr(molmtot(n),i)
          dp4 => dataptr(chi(n), i)
          dp5 => dataptr(rhoWchi(n), i)
          lo  =  lwb(get_box(rho(n), i))
          hi  =  upb(get_box(rho(n), i))
          
          select case(dm)
          case (2)
             call compute_rhoWchi_2d(dp(:,:,1,:),dp1(:,:,1,1),dp2(:,:,1,:),molmass,&
                           dp3(:,:,1,1),dp4(:,:,1,:),dp5(:,:,1,:),ng,lo,hi) 
          case (3)
             call compute_rhoWchi_3d(dp(:,:,:,:),dp1(:,:,:,1),dp2(:,:,:,:),molmass,&
                           dp3(:,:,:,1),dp4(:,:,:,:),dp5(:,:,:,:),ng,lo,hi) 
          end select
       end do
    end do
 
  end subroutine compute_rhoWchi
  
  subroutine compute_rhoWchi_2d(rho,rho_tot,molarconc,molmass,molmtot,chi,rhoWchi,ng,lo,hi)
  
    integer          :: lo(2), hi(2), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)         ! density; last dimension for species
    real(kind=dp_t)  :: rho_tot(lo(1)-ng:,lo(2)-ng:)       ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,:)   ! molar concentration
    real(kind=dp_t)  :: molmass(nspecies)                  ! species molar mass 
    real(kind=dp_t)  :: molmtot(lo(1)-ng:,lo(2)-ng:)       ! total molar mass 
    real(kind=dp_t)  :: chi(lo(1)-ng:,lo(2)-ng:,:)         ! last dimension for nspecies^2
    real(kind=dp_t)  :: rhoWchi(lo(1)-ng:,lo(2)-ng:,:) ! last dimension for nspecies^2

    ! local variables
    integer          :: i,j,row,column
    real(kind=dp_t), dimension(nspecies,nspecies) :: rhoWchiloc 
  
    ! for specific box, now start loops over alloted cells 
    do j=lo(2)-ng,hi(2)+ng
       do i=lo(1)-ng,hi(1)+ng
        
          call compute_rhoWchi_local(rho(i,j,:),rho_tot(i,j),molarconc(i,j,:),&
                          molmass,molmtot(i,j),chi(i,j,:),rhoWchi(i,j,:))

          if(.false.) then
          if(i.eq.7 .and. j.eq.14) then
          call set_Xij(rhoWchiloc, rhoWchi(i,j,:))
             do row=1, nspecies
                do column=1, nspecies
                   print*, rhoWchiloc(row,column) 
                end do
                print*, '' 
             end do
          end if
          end if

       end do
    end do

  end subroutine compute_rhoWchi_2d

  subroutine compute_rhoWchi_3d(rho,rho_tot,molarconc,molmass,molmtot,chi,rhoWchi,ng,lo,hi)

    integer          :: lo(3), hi(3), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)         ! density; last dimension for species
    real(kind=dp_t)  :: rho_tot(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)       ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)   ! molar concentration
    real(kind=dp_t)  :: molmass(nspecies)                            ! species molar mass 
    real(kind=dp_t)  :: molmtot(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)       ! total molar mass 
    real(kind=dp_t)  :: chi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)         ! last dimension for nspecies^2
    real(kind=dp_t)  :: rhoWchi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:) ! last dimension for nspecies^2
    
    ! local variables
    integer          :: i,j,k

    ! for specific box, now start loops over alloted cells 
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
       
             call compute_rhoWchi_local(rho(i,j,k,:),rho_tot(i,j,k),molarconc(i,j,k,:),&
                          molmass,molmtot(i,j,k),chi(i,j,k,:),rhoWchi(i,j,k,:))
              
         end do
      end do
    end do
   
  end subroutine compute_rhoWchi_3d
  
  subroutine compute_rhoWchi_local(rho,rho_tot,molarconc,molmass,molmtot,chi,rhoWchi)
   
    real(kind=dp_t), intent(in)   :: rho(nspecies)            
    real(kind=dp_t), intent(in)   :: rho_tot                  
    real(kind=dp_t), intent(in)   :: molarconc(nspecies) ! molar concentration
    real(kind=dp_t), intent(in)   :: molmass(nspecies)   ! species molar mass 
    real(kind=dp_t), intent(in)   :: molmtot             ! total molar mass 
    real(kind=dp_t), intent(in)   :: chi(nspecies,nspecies)   ! rank conversion done 
    real(kind=dp_t), intent(out)  :: rhoWchi(nspecies,nspecies) 
 
    ! local variables
    integer                              :: row,column
    real(kind=dp_t), dimension(nspecies) :: W !,chiw 

    ! compute massfraction W_i = rho_i/rho; 
    do row=1, nspecies  
       W(row) = rho(row)/rho_tot
    end do

    ! populate -rho*W*chi = -rho_i*chi
    do row=1, nspecies
       do column=1, nspecies
          rhoWchi(row,column) = -rho(row)*chi(row,column)  
       end do
    end do
    
  end subroutine compute_rhoWchi_local

  subroutine set_Xij(Xout_ij, Xin_ij)
        
    real(kind=dp_t), dimension(nspecies,nspecies), intent(in)  :: Xin_ij
    real(kind=dp_t), dimension(nspecies,nspecies), intent(out) :: Xout_ij  
   
    ! reshape array into matrix without doing index algebra
    Xout_ij = Xin_ij
  
  end subroutine set_Xij 
  
end module convert_mass_variables_module
