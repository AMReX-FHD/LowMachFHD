module mass_flux_utilities_module 

  use bl_space
  use multifab_module
  use ml_layout_module
  use probin_common_module, only: k_B, molmass, rhobar, molmass, nspecies
  use probin_multispecies_module, only: use_lapack, fraction_tolerance, &
                                        is_ideal_mixture, inverse_type, is_nonisothermal, &
                                        chi_iterations, avg_type
  use matrix_utilities 
  use compute_mixture_properties_module
  !use F95_LAPACK ! Donev: Disabled LAPACK so this builds more easily on different systems

  implicit none

  private

  public :: correct_rho_with_drho, &
            compute_molconc_molmtot, &
            compute_rhotot, &
            compute_Gama, &
            compute_chi, &
            compute_zeta_by_Temp, &
            compute_rhoWchi, &
            compute_sqrtLonsager_fc, &
            compute_baro_coef

contains
  
  subroutine correct_rho_with_drho(mla,rho,drho)

   type(ml_layout), intent(in   )  :: mla
   type(multifab) , intent(inout)  :: rho(:) 
   type(multifab) , intent(inout)  :: drho(:) 

   ! local variables
   integer :: lo(rho(1)%dim), hi(rho(1)%dim)
   integer :: n,i,ng_1,ng_2,dm,nlevs
 
   ! pointer for rho(nspecies), rhotot(1), molarconc(nspecies) 
   real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for rho    
   real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for drho

   type(mfiter) :: mfi
   type(box) :: tilebox
   integer :: tlo(mla%dim), thi(mla%dim)

   type(bl_prof_timer), save :: bpt

   call build(bpt,"correct_rho_with_drho")

   dm    = mla%dim     ! dimensionality
   ng_1  = rho(1)%ng   ! number of ghost cells 
   ng_2  = drho(1)%ng
   nlevs = mla%nlevel  ! number of levels 
 
   !$omp parallel private(n,i,mfi,tilebox,tlo,thi,dp1,dp2,lo,hi)

    ! loop over all boxes 
    do n=1,nlevs
       call mfiter_build(mfi, rho(n), tiling=.true.)

       do while (more_tile(mfi))
          i = get_fab_index(mfi)

          tilebox = get_growntilebox(mfi,rho(n)%ng)
          tlo = lwb(tilebox)
          thi = upb(tilebox)

!       do i=1,nfabs(rho(n))
          dp1 => dataptr(rho(n),i)
          dp2 => dataptr(drho(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))          
          select case(dm)
          case (2)
             call correct_rho_with_drho_2d(dp1(:,:,1,:),dp2(:,:,1,:),ng_1,ng_2,lo,hi,tlo,thi) 
          case (3)
             call correct_rho_with_drho_3d(dp1(:,:,:,:),dp2(:,:,:,:),ng_1,ng_2,lo,hi,tlo,thi) 
          end select
       end do
    end do
    !$omp end parallel

    call destroy(bpt)

  end subroutine correct_rho_with_drho

  subroutine correct_rho_with_drho_2d(rho,drho,ng_1,ng_2,glo,ghi,tlo,thi)
 
    integer          :: glo(2), ghi(2), ng_1, ng_2, tlo(2), thi(2)
    real(kind=dp_t)  ::  rho(glo(1)-ng_1:,glo(2)-ng_1:,:) ! density- last dim for #species
    real(kind=dp_t)  :: drho(glo(1)-ng_2:,glo(2)-ng_2:,:) ! total density in each cell 
        
    ! local variables
    integer          :: i,j
    
    ! for specific box, now start loops over alloted cells    
    do j=tlo(2), thi(2)
       do i=tlo(1), thi(1)
         
         call correct_rho_with_drho_local(rho(i,j,:),drho(i,j,:))

       end do
    end do
 
  end subroutine correct_rho_with_drho_2d

  subroutine correct_rho_with_drho_3d(rho,drho,ng_1,ng_2,glo,ghi,tlo,thi)
 
    integer          :: glo(3), ghi(3), ng_1, ng_2,tlo(:),thi(:)
    real(kind=dp_t)  ::  rho(glo(1)-ng_1:,glo(2)-ng_1:,glo(3)-ng_1:,:) ! density- last dim for #species
    real(kind=dp_t)  :: drho(glo(1)-ng_2:,glo(2)-ng_2:,glo(3)-ng_2:,:) ! total density in each cell 
    
    ! local variables
    integer          :: i,j,k
    
    ! for specific box, now start loops over alloted cells    
    do k=tlo(3), thi(3)
       do j=tlo(2), thi(2)
          do i=tlo(1), thi(1)

             call correct_rho_with_drho_local(rho(i,j,k,:),drho(i,j,k,:))

          end do
       end do
    end do
 
  end subroutine correct_rho_with_drho_3d

  subroutine correct_rho_with_drho_local(rho,drho)
 
    real(kind=dp_t), intent(inout) :: rho(nspecies)    ! density- last dim for #species
    real(kind=dp_t), intent(out)   :: drho(nspecies)   ! total density in each cell 
    
    ! local variables
    real(kind=dp_t)  :: rhotot_local
    real(kind=dp_t)  :: rho_old(nspecies)

    ! make a copy of the input rho's
    rho_old = rho
    
    ! total rho in the cell
    rhotot_local = sum(rho)

    ! make sure each rho is greater than +epsilon 
    ! to prevent division by zero in case species
    ! density, molar concentration or total density = 0.
    rho(1:nspecies) = max(rho(1:nspecies), fraction_tolerance*rhotot_local) 

    ! keep track of how much rho changed so we can reset it back later
    drho = rho - rho_old

  end subroutine correct_rho_with_drho_local 

  subroutine compute_molconc_molmtot(mla,rho,rhotot,molarconc,molmtot)
   
   type(ml_layout), intent(in   )  :: mla
   type(multifab) , intent(in   )  :: rho(:) 
   type(multifab) , intent(in   )  :: rhotot(:) 
   type(multifab) , intent(inout)  :: molarconc(:) 
   type(multifab) , intent(inout)  :: molmtot(:) 

   ! local variables
   integer :: lo(rho(1)%dim), hi(rho(1)%dim)
   integer :: n,i,dm,nlevs,ng_1,ng_2,ng_3,ng_4
 
   ! pointer for rho(nspecies), rhotot(1), molarconc(nspecies) 
   real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for rho    
   real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for rhotot
   real(kind=dp_t), pointer        :: dp3(:,:,:,:)  ! for molar concentrations
   real(kind=dp_t), pointer        :: dp4(:,:,:,:)  ! for total molar concentrations

   type(mfiter) :: mfi
   type(box) :: tilebox
   integer :: tlo(mla%dim), thi(mla%dim)

   type(bl_prof_timer), save :: bpt

   call build(bpt,"compute_molconc_molmtot")

   dm    = mla%dim     ! dimensionality
   ng_1  = rho(1)%ng   ! number of ghost cells 
   ng_2  = rhotot(1)%ng
   ng_3  = molarconc(1)%ng
   ng_4  = molmtot(1)%ng
   nlevs = mla%nlevel  ! number of levels 

   if (ng_3 .ne. ng_4) then
      call bl_error('compute_molconc_molmtot: ng for molarconc and molmass differ')
   end if
 
   !$omp parallel private(n,i,mfi,tilebox,tlo,thi,dp1,dp2,dp3,dp4,lo,hi)

    ! loop over all boxes 
    do n=1,nlevs
       call mfiter_build(mfi, rho(n), tiling=.true.)

       do while (more_tile(mfi))
          i = get_fab_index(mfi)

          tilebox = get_growntilebox(mfi,molarconc(n)%ng)
          tlo = lwb(tilebox)
          thi = upb(tilebox)

 !      do i=1,nfabs(rho(n))
          dp1 => dataptr(rho(n),i)
          dp2 => dataptr(rhotot(n),i)
          dp3 => dataptr(molarconc(n),i)
          dp4 => dataptr(molmtot(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          select case(dm)
          case (2)
             call compute_molconc_molmtot_2d(dp1(:,:,1,:),dp2(:,:,1,1), &
                                             dp3(:,:,1,:),dp4(:,:,1,1), &
                                             ng_1,ng_2,ng_3,ng_4,lo,hi,tlo,thi) 
          case (3)
             call compute_molconc_molmtot_3d(dp1(:,:,:,:),dp2(:,:,:,1), &
                                             dp3(:,:,:,:),dp4(:,:,:,1), &
                                             ng_1,ng_2,ng_3,ng_4,lo,hi,tlo,thi) 
          end select
       end do
    end do
    !$omp end parallel

    call destroy(bpt)

  end subroutine compute_molconc_molmtot

  subroutine compute_molconc_molmtot_2d(rho,rhotot,molarconc,molmtot, &
                                        ng_1,ng_2,ng_3,ng_4,glo,ghi,tlo,thi)
 
    integer          :: glo(2), ghi(2), ng_1, ng_2, ng_3, ng_4
    integer          :: tlo(2), thi(2)
    real(kind=dp_t)  ::       rho(glo(1)-ng_1:,glo(2)-ng_1:,:)       ! density- last dim for #species
    real(kind=dp_t)  ::    rhotot(glo(1)-ng_2:,glo(2)-ng_2:)     ! total density in each cell 
    real(kind=dp_t)  :: molarconc(glo(1)-ng_3:,glo(2)-ng_3:,:) ! molar concentration
    real(kind=dp_t)  ::   molmtot(glo(1)-ng_4:,glo(2)-ng_4:)     ! total molar mass 
        
    ! local variables
    integer          :: i,j
    
    ! for specific box, now start loops over alloted cells    
    do j=tlo(2), thi(2)
       do i=tlo(1), thi(1)
         
         call compute_molconc_molmtot_local(rho(i,j,:),rhotot(i,j), &
                                            molarconc(i,j,:),molmtot(i,j))

       end do
    end do
 
  end subroutine compute_molconc_molmtot_2d

  subroutine compute_molconc_molmtot_3d(rho,rhotot,molarconc,molmtot, &
                                               ng_1,ng_2,ng_3,ng_4,glo,ghi,tlo,thi)
 
    integer          :: glo(3), ghi(3), ng_1, ng_2, ng_3, ng_4
    integer          :: tlo(3), thi(3)
    real(kind=dp_t)  ::       rho(glo(1)-ng_1:,glo(2)-ng_1:,glo(3)-ng_1:,:) ! density- last dim for #species
    real(kind=dp_t)  ::    rhotot(glo(1)-ng_2:,glo(2)-ng_2:,glo(3)-ng_2:)   ! total density in each cell 
    real(kind=dp_t)  :: molarconc(glo(1)-ng_3:,glo(2)-ng_3:,glo(3)-ng_3:,:) ! molar concentration
    real(kind=dp_t)  ::   molmtot(glo(1)-ng_4:,glo(2)-ng_4:,glo(3)-ng_4:)   ! total molar mass 
    
    ! local variables
    integer          :: i,j,k
    
    ! for specific box, now start loops over alloted cells    
    do k=tlo(3), thi(3)
       do j=tlo(2), thi(2)
          do i=tlo(1), thi(1)

             call compute_molconc_molmtot_local(rho(i,j,k,:),rhotot(i,j,k),&
                                                molarconc(i,j,k,:),molmtot(i,j,k))

          end do
       end do
    end do
 
  end subroutine compute_molconc_molmtot_3d

  subroutine compute_molconc_molmtot_local(rho,rhotot,molarconc,molmtot)
 
    real(kind=dp_t), intent(in)   :: rho(nspecies)       ! density- last dim for #species
    real(kind=dp_t), intent(in)   :: rhotot              ! total density in each cell 
    real(kind=dp_t), intent(out)  :: molarconc(nspecies) ! molar concentration
    real(kind=dp_t), intent(out)  :: molmtot             ! total molar mass 
    
    ! local variables
    integer          :: n
    real(kind=dp_t), dimension(nspecies) :: W            ! mass fraction w_i = rho_i/rho 
    real(kind=dp_t)  :: Sum_woverm

    ! calculate mass fraction and total molar mass (1/m=Sum(w_i/m_i))
    Sum_woverm=0.d0
    do n=1, nspecies  
       W(n) = rho(n)/rhotot
       Sum_woverm = Sum_woverm + W(n)/molmass(n)
    end do
    molmtot = 1.0d0/Sum_woverm 

    ! calculate molar concentrations in each cell (x_i=m*w_i/m_i) 
    do n=1, nspecies 
       molarconc(n) = molmtot*W(n)/molmass(n)
    end do
    
  end subroutine compute_molconc_molmtot_local 

  subroutine compute_rhotot(mla,rho,rhotot,ghost_cells_in)
   
    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: rho(:) 
    type(multifab) , intent(inout)  :: rhotot(:) 
    logical, intent(in), optional   :: ghost_cells_in

    ! local variables
    integer :: n,nlevs,comp,ng

    logical ghost_cells

    type(bl_prof_timer), save :: bpt

    call build(bpt, "compute_rhotot")

    nlevs = mla%nlevel  ! number of levels 

    ghost_cells = .false.
    if (present(ghost_cells_in)) ghost_cells = ghost_cells_in

    ng=0
    if (ghost_cells) then
       ng=rhotot(1)%ng
    end if

    do n=1,nlevs
       call multifab_setval(rhotot(n),0.d0,all=.true.)
       do comp=1,nspecies
          call multifab_plus_plus_c(rhotot(n),1,rho(n),comp,1,ng)
       end do
    end do

  end subroutine compute_rhotot

  subroutine compute_Gama(mla,molarconc,Hessian,Gama)

    type(ml_layout), intent(in   )  :: mla
    type(multifab),  intent(in   )  :: molarconc(:) 
    type(multifab),  intent(in   )  :: Hessian(:)    ! Hessian matrix
    type(multifab),  intent(inout)  :: Gama(:)       ! Non-ideality coeficient 
 
    ! local variables
    integer :: lo(molarconc(1)%dim), hi(molarconc(1)%dim)
    integer :: n,i,dm,nlevs,ng_2,ng_4,ng_5

    ! assign pointers for multifabs to be passed
    real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for molarconc
    real(kind=dp_t), pointer        :: dp4(:,:,:,:)  ! for Hessian 
    real(kind=dp_t), pointer        :: dp5(:,:,:,:)  ! for Gama 

    type(mfiter) :: mfi
    type(box) :: tilebox
    integer :: tlo(mla%dim), thi(mla%dim)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "compute_Gama")

    dm    = mla%dim     ! dimensionality
    ng_2  = molarconc(1)%ng
    ng_4  = Hessian(1)%ng
    ng_5  = Gama(1)%ng
    nlevs = mla%nlevel  ! number of levels 
  
   !$omp parallel private(n,i,mfi,tilebox,tlo,thi) &
   !$omp private(dp2,dp4,dp5,lo,hi)

    ! loop over all boxes 
    do n=1,nlevs
       call mfiter_build(mfi, Gama(n), tiling=.true.)

       do while (more_tile(mfi))
          i = get_fab_index(mfi)

          tilebox = get_growntilebox(mfi,Gama(n)%ng)
          tlo = lwb(tilebox)
          thi = upb(tilebox)
          dp2 => dataptr(molarconc(n),i)
          dp4 => dataptr(Hessian(n),i)
          dp5 => dataptr(Gama(n),i)
          lo = lwb(get_box(Gama(n),i))
          hi = upb(get_box(Gama(n),i))          
          select case(dm)
          case (2)
             call compute_Gama_2d(dp2(:,:,1,:),dp4(:,:,1,:),dp5(:,:,1,:), &
                                  ng_2,ng_4,ng_5,lo,hi,tlo,thi) 
          case (3)
             call compute_Gama_3d(dp2(:,:,:,:),dp4(:,:,:,:),dp5(:,:,:,:), &
                                  ng_2,ng_4,ng_5,lo,hi,tlo,thi) 
          end select
       end do
    end do
    !$omp end parallel

    call destroy(bpt)

  end subroutine compute_Gama
  
  subroutine compute_Gama_2d(molarconc,Hessian,Gama, &
                             ng_2,ng_4,ng_5,glo,ghi,tlo,thi)

    integer          :: glo(2), ghi(2), ng_2,ng_4,ng_5
    integer          :: tlo(2), thi(2)
    real(kind=dp_t)  :: molarconc(glo(1)-ng_2:,glo(2)-ng_2:,:) ! molar concentration 
    real(kind=dp_t)  ::   Hessian(glo(1)-ng_4:,glo(2)-ng_4:,:) ! last dimension for nspecies^2
    real(kind=dp_t)  ::      Gama(glo(1)-ng_5:,glo(2)-ng_5:,:) ! last dimension for nspecies^2

    ! local varialbes
    integer          :: i,j

    ! for specific box, now start loops over alloted cells 
    do j=tlo(2),thi(2)
       do i=tlo(1),thi(1)
       
          call compute_Gama_local(molarconc(i,j,:),Hessian(i,j,:),Gama(i,j,:))

       end do
    end do
   
  end subroutine compute_Gama_2d

  subroutine compute_Gama_3d(molarconc,Hessian,Gama, &
                             ng_2,ng_4,ng_5,glo,ghi,tlo,thi)
 
    integer          :: glo(3), ghi(3), ng_2,ng_4,ng_5
    integer          :: tlo(3), thi(3)
    real(kind=dp_t)  :: molarconc(glo(1)-ng_2:,glo(2)-ng_2:,glo(3)-ng_2:,:) ! molar concentration; 
    real(kind=dp_t)  ::   Hessian(glo(1)-ng_4:,glo(2)-ng_4:,glo(3)-ng_4:,:) ! last dimension for nspecies^2
    real(kind=dp_t)  ::      Gama(glo(1)-ng_5:,glo(2)-ng_5:,glo(3)-ng_5:,:) ! last dimension for nspecies^2

    ! local varialbes
    integer          :: i,j,k

    ! for specific box, now start loops over alloted cells 
    do k=tlo(3),thi(3)
       do j=tlo(2),thi(2)
          do i=tlo(1),thi(1)

             call compute_Gama_local(molarconc(i,j,k,:),Hessian(i,j,k,:),Gama(i,j,k,:))
          end do
       end do
    end do
   
  end subroutine compute_Gama_3d

  subroutine compute_Gama_local(molarconc,Hessian,Gama)
   
    real(kind=dp_t), intent(in)   :: molarconc(nspecies)
    real(kind=dp_t), intent(in)   :: Hessian(nspecies,nspecies)
    real(kind=dp_t), intent(out)  :: Gama(nspecies,nspecies)
 
    ! local variables
    integer                                       :: row,column
    real(kind=dp_t), dimension(nspecies,nspecies) :: I, X_xxT

    ! Identity matrix
    I = 0.d0
    do row=1,nspecies
       I(row, row) = 1.d0        
    end do

    ! populate X_xxT
    if (is_ideal_mixture) then
       X_xxT = 0.d0
    else
       do row=1, nspecies  
          ! diagonal entries
          X_xxT(row,row) = molarconc(row) - molarconc(row)**2 
          do column=1, row-1
             ! off-diagnoal entries
             X_xxT(row,column)   = -molarconc(row)*molarconc(column)  ! form x*transpose(x) off diagonals 
             X_xxT(column, row)  = X_xxT(row, column)                 ! symmetric
          end do
       end do
    end if
  
    ! compute Gama 
    Gama = I + matmul(X_xxT, Hessian)     
 
  end subroutine compute_Gama_local

  subroutine compute_chi(mla,rho,rhotot,molarconc,chi,D_bar)
   
    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: rho(:)
    type(multifab) , intent(in   )  :: rhotot(:) 
    type(multifab) , intent(in   )  :: molarconc(:)
    type(multifab) , intent(inout)  :: chi(:) 
    type(multifab) , intent(in   )  :: D_bar(:) 

    ! local variables
    integer :: lo(rho(1)%dim), hi(rho(1)%dim)
    integer :: n,i,dm,nlevs,ng_0,ng_1,ng_2,ng_3,ng_4
 
    ! pointer for rho(nspecies), molarconc(nspecies) 
    real(kind=dp_t), pointer        :: dp0(:,:,:,:)  ! for rho    
    real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for rhotot
    real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for molarconc
    real(kind=dp_t), pointer        :: dp3(:,:,:,:)  ! for chi
    real(kind=dp_t), pointer        :: dp4(:,:,:,:)  ! for D_bar

    type(mfiter) :: mfi
    type(box) :: tilebox
    integer :: tlo(mla%dim), thi(mla%dim)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"compute_chi")

    dm = mla%dim        ! dimensionality
    ng_0 = rho(1)%ng    ! number of ghost cells 
    ng_1 = rhotot(1)%ng
    ng_2 = molarconc(1)%ng
    ng_3 = chi(1)%ng
    ng_4 = D_bar(1)%ng
    nlevs = mla%nlevel  ! number of levels 
 
    !$omp parallel private(n,i,mfi,tilebox,tlo,thi) &
    !$omp private(dp0,dp1,dp2,dp3,dp4,lo,hi)

    ! loop over all boxes 
    do n=1,nlevs
       call mfiter_build(mfi, rho(n), tiling=.true.)

       do while (more_tile(mfi))
          i = get_fab_index(mfi)

          tilebox = get_growntilebox(mfi,chi(n)%ng)
          tlo = lwb(tilebox)
          thi = upb(tilebox)

!       do i=1,nfabs(rho(n))
          dp0  => dataptr(rho(n), i)
          dp1 => dataptr(rhotot(n), i)
          dp2 => dataptr(molarconc(n), i)
          dp3 => dataptr(chi(n), i)
          dp4 => dataptr(D_bar(n), i)
          lo  =  lwb(get_box(rho(n), i))
          hi  =  upb(get_box(rho(n), i))          
          select case(dm)
          case (2)
             call compute_chi_2d(dp0(:,:,1,:),dp1(:,:,1,1),dp2(:,:,1,:),dp3(:,:,1,:),&
                                 dp4(:,:,1,:),&
                                 ng_0,ng_1,ng_2,ng_3,ng_4,lo,hi,tlo,thi) 
          case (3)
             call compute_chi_3d(dp0(:,:,:,:),dp1(:,:,:,1),dp2(:,:,:,:),dp3(:,:,:,:),&
                                 dp4(:,:,:,:),&
                                 ng_0,ng_1,ng_2,ng_3,ng_4,lo,hi,tlo,thi) 
          end select
       end do
    end do
    !$omp end parallel

    call destroy(bpt)

  end subroutine compute_chi
 
  subroutine compute_chi_2d(rho,rhotot,molarconc,chi,D_bar, &
                            ng_0,ng_1,ng_2,ng_3,ng_4,glo,ghi,tlo,thi)

    integer          :: glo(2), ghi(2), ng_0,ng_1,ng_2,ng_3,ng_4
    integer          :: tlo(2), thi(2)
    real(kind=dp_t)  ::          rho(glo(1)-ng_0:,glo(2)-ng_0:,:) ! density; last dimension for species
    real(kind=dp_t)  ::       rhotot(glo(1)-ng_1:,glo(2)-ng_1:)   ! total density in each cell 
    real(kind=dp_t)  ::    molarconc(glo(1)-ng_2:,glo(2)-ng_2:,:) ! molar concentration 
    real(kind=dp_t)  ::          chi(glo(1)-ng_3:,glo(2)-ng_3:,:) ! last dimension for nspecies^2
    real(kind=dp_t)  ::        D_bar(glo(1)-ng_4:,glo(2)-ng_4:,:) ! MS diff-coefs 

    ! local variables
    integer          :: i,j

    ! for specific box, now start loops over alloted cells 
    do j=tlo(2),thi(2)
       do i=tlo(1),thi(1)
    
          call compute_chi_local(rho(i,j,:),rhotot(i,j),molarconc(i,j,:),&
                                 chi(i,j,:),D_bar(i,j,:))

       end do
    end do

  end subroutine compute_chi_2d

  subroutine compute_chi_3d(rho,rhotot,molarconc,chi,D_bar, &
                            ng_0,ng_1,ng_2,ng_3,ng_4,glo,ghi,tlo,thi)
   
    integer          :: glo(3), ghi(3), ng_0,ng_1,ng_2,ng_3,ng_4
    integer          :: tlo(3), thi(3)
    real(kind=dp_t)  ::          rho(glo(1)-ng_0:,glo(2)-ng_0:,glo(3)-ng_0:,:) ! density; last dimension for species
    real(kind=dp_t)  ::       rhotot(glo(1)-ng_1:,glo(2)-ng_1:,glo(3)-ng_1:)   ! total density in each cell 
    real(kind=dp_t)  ::    molarconc(glo(1)-ng_2:,glo(2)-ng_2:,glo(3)-ng_2:,:) ! molar concentration; 
    real(kind=dp_t)  ::          chi(glo(1)-ng_3:,glo(2)-ng_3:,glo(3)-ng_3:,:) ! last dimension for nspecies^2
    real(kind=dp_t)  ::        D_bar(glo(1)-ng_4:,glo(2)-ng_4:,glo(3)-ng_4:,:) ! SM diffusion constants 
    
    ! local variables
    integer          :: i,j,k

    ! for specific box, now start loops over alloted cells 
    do k=tlo(3),thi(3)
       do j=tlo(2),thi(2)
          do i=tlo(1),thi(1)
       
             call compute_chi_local(rho(i,j,k,:),rhotot(i,j,k),molarconc(i,j,k,:),&
                                    chi(i,j,k,:),D_bar(i,j,k,:))

          end do
       end do
    end do
   
  end subroutine compute_chi_3d

  subroutine compute_chi_local(rho,rhotot,molarconc,chi,D_bar)
    
    real(kind=dp_t), intent(in   ) :: rho(nspecies)         
    real(kind=dp_t), intent(in   ) :: rhotot               
    real(kind=dp_t), intent(in   ) :: molarconc(nspecies) 
    real(kind=dp_t), intent(inout) :: chi(nspecies,nspecies)   
    real(kind=dp_t), intent(in   ) :: D_bar(nspecies,nspecies) 

    ! local variables
    integer                         :: row,column
    real(kind=dp_t)                 :: Sum_knoti   

    ! vectors and matrices to be used by LAPACK 
    real(kind=dp_t), dimension(nspecies,nspecies) :: Lambda
    real(kind=dp_t), dimension(nspecies)          :: W

    ! compute chi either selecting inverse/pseudoinverse or iterative methods 
    if (use_lapack) then

       ! compute Lambda_ij matrix and massfraction W_i = rho_i/rho; molarconc is 
       ! expressed in terms of molmtot,mi,rhotot etc. 
       do row=1, nspecies  
          do column=1, row-1
             Lambda(row, column) = -molarconc(row)*molarconc(column)/D_bar(row,column)
             Lambda(column, row) = Lambda(row, column) 
          end do
          W(row) = rho(row)/rhotot
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

       call compute_chi_lapack(Lambda(:,:),chi(:,:),W(:))

    else

       call Dbar2chi_iterative(chi_iterations,D_bar(:,:),molarconc(:),chi(:,:)) 

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
    !stop "LAPACK95 dget? disabled"

    ! populate chi with B^(-1)
    chi = chilocal   
 
    case(2) 
    !==========================================================
    ! Using pseudoinverse 
    !==========================================================

    ! SVD decomposition of chilocal = U * S * VTranspose; note that chilocal 
    ! is changed. also V=(VT)T, UT = (U)T are needed for pseudoinverse of chilocal.
    !stop "LAPACK95 la_gesvd disabled"
    call la_gesvd(chilocal, S, U, VT)
    V = transpose(VT)
    UT = transpose(U)
   
    ! populate diagonal matrix Sdag = 1/S with diagonal=0 below a chosen tolerance
    do row=1, nspecies
       do column=1,nspecies
          Sdag(row,column) = 0.0d0
       end do
       
       if(S(row).gt.fraction_tolerance*sum(S)) then 
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

  subroutine compute_zeta_by_Temp(mla,molarconc,D_bar,D_therm,Temp,zeta_by_Temp)
   
    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: molarconc(:)
    type(multifab) , intent(in   )  :: D_bar(:) 
    type(multifab) , intent(in   )  :: D_therm(:) 
    type(multifab) , intent(in   )  :: Temp(:) 
    type(multifab) , intent(inout)  :: zeta_by_Temp(:) 

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: n,i,dm,nlevs,ng_2,ng_4,ng_5,ng_6,ng_7

    real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for molarconc
    real(kind=dp_t), pointer        :: dp4(:,:,:,:)  ! for D_bar
    real(kind=dp_t), pointer        :: dp5(:,:,:,:)  ! for Temp
    real(kind=dp_t), pointer        :: dp6(:,:,:,:)  ! for zeta_by_Temp
    real(kind=dp_t), pointer        :: dp7(:,:,:,:)  ! for D_therm

    type(mfiter) :: mfi
    type(box) :: tilebox
    integer :: tlo(mla%dim), thi(mla%dim)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"compute_zeta_by_Temp")

    dm = mla%dim        ! dimensionality
    ng_2 = molarconc(1)%ng
    ng_4 = D_bar(1)%ng
    ng_5 = Temp(1)%ng
    ng_6 = zeta_by_Temp(1)%ng
    ng_7 = D_therm(1)%ng
    nlevs = mla%nlevel  ! number of levels 
 
    !$omp parallel private(n,i,mfi,tilebox,tlo,thi) &
    !$omp private(dp0,dp1,dp2,dp4,dp5,dp6,dp7,lo,hi)

    ! loop over all boxes 
    do n=1,nlevs
       call mfiter_build(mfi, zeta_by_Temp(n), tiling=.true.)

       do while (more_tile(mfi))
          i = get_fab_index(mfi)

          tilebox = get_growntilebox(mfi,zeta_by_Temp(n)%ng)
          tlo = lwb(tilebox)
          thi = upb(tilebox)

!       do i=1,nfabs(rho(n))
          dp2 => dataptr(molarconc(n), i)
          dp4 => dataptr(D_bar(n), i)
          dp5 => dataptr(Temp(n), i)
          dp6 => dataptr(zeta_by_Temp(n), i)
          dp7 => dataptr(D_therm(n), i)
          lo  =  lwb(get_box(zeta_by_Temp(n), i))
          hi  =  upb(get_box(zeta_by_Temp(n), i))          
          select case(dm)
          case (2)
             call compute_zeta_by_Temp_2d(dp2(:,:,1,:),&
                                 dp4(:,:,1,:),dp5(:,:,1,1),dp6(:,:,1,:),dp7(:,:,1,:),&
                                 ng_2,ng_4,ng_5,ng_6,ng_7,lo,hi,tlo,thi) 
          case (3)
             call compute_zeta_by_Temp_3d(dp2(:,:,:,:),&
                                 dp4(:,:,:,:),dp5(:,:,:,1),dp6(:,:,:,:),dp7(:,:,:,:),&
                                 ng_2,ng_4,ng_5,ng_6,ng_7,lo,hi,tlo,thi) 
          end select
       end do
    end do
    !$omp end parallel

    call destroy(bpt)

  end subroutine compute_zeta_by_Temp
 
  subroutine compute_zeta_by_Temp_2d(molarconc,D_bar,Temp,zeta_by_Temp,D_therm, &
                                     ng_2,ng_4,ng_5,ng_6,ng_7,glo,ghi,tlo,thi)

    integer          :: glo(2), ghi(2), ng_2,ng_4,ng_5,ng_6,ng_7
    integer          :: tlo(2), thi(2)
    real(kind=dp_t)  ::    molarconc(glo(1)-ng_2:,glo(2)-ng_2:,:) ! molar concentration 
    real(kind=dp_t)  ::        D_bar(glo(1)-ng_4:,glo(2)-ng_4:,:) ! MS diff-coefs 
    real(kind=dp_t)  ::         Temp(glo(1)-ng_5:,glo(2)-ng_5:)   ! Temperature 
    real(kind=dp_t)  :: zeta_by_Temp(glo(1)-ng_6:,glo(2)-ng_6:,:) ! zeta/T
    real(kind=dp_t)  ::      D_therm(glo(1)-ng_7:,glo(2)-ng_7:,:) ! thermo diff-coefs 

    ! local variables
    integer          :: i,j

    ! for specific box, now start loops over alloted cells 
    do j=tlo(2),thi(2)
       do i=tlo(1),thi(1)
    
          call compute_zeta_by_Temp_local(molarconc(i,j,:),&
                                          D_bar(i,j,:),Temp(i,j),zeta_by_Temp(i,j,:),&
                                          D_therm(i,j,:))

       end do
    end do

  end subroutine compute_zeta_by_Temp_2d

  subroutine compute_zeta_by_Temp_3d(molarconc,D_bar,Temp,zeta_by_Temp,D_therm, &
                                     ng_2,ng_4,ng_5,ng_6,ng_7,glo,ghi,tlo,thi)
   
    integer          :: glo(3), ghi(3), ng_2,ng_4,ng_5,ng_6,ng_7
    integer          :: tlo(3), thi(3)
    real(kind=dp_t)  ::    molarconc(glo(1)-ng_2:,glo(2)-ng_2:,glo(3)-ng_2:,:) ! molar concentration; 
    real(kind=dp_t)  ::        D_bar(glo(1)-ng_4:,glo(2)-ng_4:,glo(3)-ng_4:,:) ! SM diffusion constants 
    real(kind=dp_t)  ::         Temp(glo(1)-ng_5:,glo(2)-ng_5:,glo(3)-ng_5:)   ! Temperature 
    real(kind=dp_t)  :: zeta_by_Temp(glo(1)-ng_6:,glo(2)-ng_6:,glo(3)-ng_6:,:) ! zeta/T
    real(kind=dp_t)  ::      D_therm(glo(1)-ng_7:,glo(2)-ng_7:,glo(3)-ng_7:,:) ! thermo diffusion constants 
    
    ! local variables
    integer          :: i,j,k

    ! for specific box, now start loops over alloted cells 
    do k=tlo(3),thi(3)
       do j=tlo(2),thi(2)
          do i=tlo(1),thi(1)
       
             call compute_zeta_by_Temp_local(molarconc(i,j,k,:),&
                                             D_bar(i,j,k,:),Temp(i,j,k),zeta_by_Temp(i,j,k,:),&
                                             D_therm(i,j,k,:))

          end do
       end do
    end do
   
  end subroutine compute_zeta_by_Temp_3d

  subroutine compute_zeta_by_Temp_local(molarconc,D_bar,Temp,zeta_by_Temp,D_therm)
    
    real(kind=dp_t), intent(in   ) :: molarconc(nspecies) 
    real(kind=dp_t), intent(in   ) :: D_bar(nspecies,nspecies) 
    real(kind=dp_t), intent(in   ) :: Temp
    real(kind=dp_t), intent(inout) :: zeta_by_Temp(nspecies)
    real(kind=dp_t), intent(in   ) :: D_therm(nspecies)

    ! local variables
    integer                         :: row,column
    real(kind=dp_t)                 :: Sum_knoti   

    ! vectors and matrices to be used by LAPACK 
    real(kind=dp_t), dimension(nspecies,nspecies) :: Lambda

    ! compute zeta_by_Temp for thermodiffusion
    if(is_nonisothermal) then

       ! compute Lambda_ij matrix; molarconc is 
       ! expressed in terms of molmtot,mi,rhotot etc. 
       do row=1, nspecies  
          do column=1, row-1
             Lambda(row, column) = -molarconc(row)*molarconc(column)/D_bar(row,column)
             Lambda(column, row) = Lambda(row, column) 
          end do
       end do

       do row=1, nspecies
          Sum_knoti = 0.d0
          do column=1, nspecies
             if(column.ne.row) then
                Sum_knoti = Sum_knoti + Lambda(row,column)*(D_therm(row)-D_therm(column))
             end if
             zeta_by_Temp(row) = Sum_knoti/Temp
          end do
       end do

    end if

  end subroutine compute_zeta_by_Temp_local


  subroutine compute_sqrtLonsager_fc(mla,rho,rhotot,sqrtLonsager_fc,dx)
 
    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: rho(:)
    type(multifab) , intent(in   )  :: rhotot(:) 
    type(multifab) , intent(inout)  :: sqrtLonsager_fc(:,:)
    real(kind=dp_t), intent(in   )   :: dx(:,:)

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: n,i,dm,nlevs,ng_0,ng_1,ng_2
 
    ! pointer for rho(nspecies), molarconc(nspecies) 
    real(kind=dp_t), pointer        :: dp0(:,:,:,:)  ! for rho    
    real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for rhotot
    real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for sqrtLonsager_x
    real(kind=dp_t), pointer        :: dp3(:,:,:,:)  ! for sqrtLonsager_y
    real(kind=dp_t), pointer        :: dp4(:,:,:,:)  ! for sqrtLonsager_z

    type(bl_prof_timer), save :: bpt

    call build(bpt, "compute_sqrtLonsager_fc")

    nlevs = mla%nlevel  ! number of levels 
    dm = mla%dim        ! dimensionality
    ng_0 = rho(1)%ng    ! number of ghost cells 
    ng_1 = rhotot(1)%ng
    ng_2 = sqrtLonsager_fc(1,1)%ng

    ! loop over all boxes 
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp0 => dataptr(rho(n), i)
          dp1 => dataptr(rhotot(n), i)
          dp2 => dataptr(sqrtLonsager_fc(n,1), i)
          dp3 => dataptr(sqrtLonsager_fc(n,2), i)
          lo  =  lwb(get_box(rho(n), i))
          hi  =  upb(get_box(rho(n), i))         
          select case(dm)
          case (2)
             call compute_sqrtLonsager_2d(dp0(:,:,1,:),dp1(:,:,1,1), &
                                          dp2(:,:,1,:),dp3(:,:,1,:),&
                                          ng_0,ng_1,ng_2,lo,hi,dx(1,:)) 
          case (3)
             dp4 => dataptr(sqrtLonsager_fc(n,3), i)
             call compute_sqrtLonsager_3d(dp0(:,:,:,:),dp1(:,:,:,1), &
                                          dp2(:,:,:,:),dp3(:,:,:,:),dp4(:,:,:,:), &
                                          ng_0,ng_1,ng_2,lo,hi,dx(1,:))
          end select
       end do
    end do
 
    call destroy(bpt)

  end subroutine compute_sqrtLonsager_fc
  
  subroutine compute_sqrtLonsager_2d(rho,rhotot,sqrtLonsager_x,sqrtLonsager_y, &
                                     ng_0,ng_1,ng_2,lo,hi,dx)
  
    integer          :: lo(2), hi(2), ng_0, ng_1, ng_2
    real(kind=dp_t)  ::            rho(lo(1)-ng_0:,lo(2)-ng_0:,:) ! density; last dimension for species
    real(kind=dp_t)  ::         rhotot(lo(1)-ng_1:,lo(2)-ng_1:)   ! total density in each cell 
    real(kind=dp_t)  :: sqrtLonsager_x(lo(1)-ng_2:,lo(2)-ng_2:,:) ! last dimension for nspecies^2
    real(kind=dp_t)  :: sqrtLonsager_y(lo(1)-ng_2:,lo(2)-ng_2:,:) ! last dimension for nspecies^2
    real(kind=dp_t)  :: dx(:)

    ! local variables
    integer         :: i,j
    real(kind=dp_t) :: rhoav(nspecies)

    ! x-faces
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)+1
          call compute_nonnegative_rho_av(rho(i-1,j,:), rho(i,j,:), rhoav, dx)
          call compute_sqrtLonsager_local(rhoav,sum(rhoav),sqrtLonsager_x(i,j,:))
    end do
    end do

    ! y-faces
    do j=lo(2),hi(2)+1
    do i=lo(1),hi(1)
          call compute_nonnegative_rho_av(rho(i,j-1,:), rho(i,j,:), rhoav, dx)
          call compute_sqrtLonsager_local(rhoav,sum(rhoav),sqrtLonsager_y(i,j,:))
    end do
    end do

  end subroutine compute_sqrtLonsager_2d

  subroutine compute_sqrtLonsager_3d(rho,rhotot,sqrtLonsager_x,sqrtLonsager_y,sqrtLonsager_z, &
                                     ng_0,ng_1,ng_2,lo,hi,dx)

    integer          :: lo(3), hi(3), ng_0, ng_1, ng_2
    real(kind=dp_t)  ::            rho(lo(1)-ng_0:,lo(2)-ng_0:,lo(3)-ng_0:,:) ! density; last dimension for species
    real(kind=dp_t)  ::         rhotot(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)   ! total density in each cell 
    real(kind=dp_t)  :: sqrtLonsager_x(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_2:,:) ! last dimension for nspecies^2
    real(kind=dp_t)  :: sqrtLonsager_y(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_2:,:) ! last dimension for nspecies^2
    real(kind=dp_t)  :: sqrtLonsager_z(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_2:,:) ! last dimension for nspecies^2
    real(kind=dp_t)  :: dx(:)

    ! local variables
    integer         :: i,j,k
    real(kind=dp_t) :: rhoav(nspecies)

    ! x-faces
    do k=lo(3),hi(3)
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)+1
          call compute_nonnegative_rho_av(rho(i-1,j,k,:), rho(i,j,k,:), rhoav, dx)
          call compute_sqrtLonsager_local(rhoav,sum(rhoav),sqrtLonsager_x(i,j,k,:))
    end do
    end do
    end do

    ! y-faces
    do k=lo(3),hi(3)
    do j=lo(2),hi(2)+1
    do i=lo(1),hi(1)
          call compute_nonnegative_rho_av(rho(i,j-1,k,:), rho(i,j,k,:), rhoav, dx)
          call compute_sqrtLonsager_local(rhoav,sum(rhoav),sqrtLonsager_y(i,j,k,:))
    end do
    end do
    end do

    ! z-faces
    do k=lo(3),hi(3)+1
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)
          call compute_nonnegative_rho_av(rho(i,j,k-1,:), rho(i,j,k,:), rhoav, dx)
          call compute_sqrtLonsager_local(rhoav,sum(rhoav),sqrtLonsager_z(i,j,k,:))
    end do
    end do
    end do
   
  end subroutine compute_sqrtLonsager_3d

  subroutine compute_nonnegative_rho_av(rho1, rho2, rhoav, dx)
    real(kind=dp_t), intent(in   ) :: rho1(nspecies), rho2(nspecies) ! Densities in two neighboring cells
    real(kind=dp_t), intent(  out) :: rhoav(nspecies)                ! Face-centered average  
    real(kind=dp_t), intent(in   ) :: dx(:)

    real(kind=dp_t) :: dv, value1, value2, tmp1, tmp2
    integer :: comp

    ! cell volume
    dv = product(dx(1:MAX_SPACEDIM))

    do comp=1,nspecies
      value1 = rho1(comp)/molmass(comp) ! Convert to number density
      value2 = rho2(comp)/molmass(comp)

      select case(avg_type)
      case(1) ! Arithmetic with a C0-smoothed Heaviside
        if ( (value1 .le. 0.d0) .or. (value2 .le. 0.d0) ) then
          rhoav(comp)=0.d0
        else
          tmp1=min(dv*value1,1.d0)
          tmp2=min(dv*value2,1.d0)
          rhoav(comp)=molmass(comp)*(value1+value2)/2.d0*tmp1*tmp2
        end if
      case(2) ! Geometric
        rhoav(comp)=molmass(comp)*sqrt(max(value1,0.d0)*max(value2,0.d0))
      case(3) ! Harmonic
        ! What we want here is the harmonic mean of max(value1,0) and max(value2,0)
        ! Where we define the result to be zero if either one is zero
        ! But numerically we want to avoid here division by zero
        if ( (value1 .le. 10.d0*tiny(1.d0)) .or. (value2 .le. 10.d0*tiny(1.d0)) ) then
          rhoav(comp)=0.d0
        else
          rhoav(comp)=molmass(comp)*2.d0/(1.d0/value1+1.d0/value2)
        end if
      case(10) ! Arithmetic with (discontinuous) Heaviside
        if ( (value1 .le. 0.d0) .or. (value2 .le. 0.d0) ) then
          rhoav(comp)=0.d0
        else
          rhoav(comp)=molmass(comp)*(value1+value2)/2.d0
        end if
      case(11) ! Arithmetic with C1-smoothed Heaviside
        if ( (value1 .le. 0.d0) .or. (value2 .le. 0.d0) ) then
          rhoav(comp)=0.d0
        else
          tmp1=dv*value1
          if (tmp1<1.d0) then
            tmp1=(3.d0-2.d0*tmp1)*tmp1**2
          else
            tmp1=1.d0
          end if
          tmp2=dv*value2
          if (tmp2<1.d0) then
            tmp2=(3.d0-2.d0*tmp2)*tmp2**2
          else
            tmp2=1.d0
          end if
          rhoav(comp)=molmass(comp)*(value1+value2)/2.d0*tmp1*tmp2
        endif
      case(12) ! Arithmetic with C2-smoothed Heaviside
        if ( (value1 .le. 0.d0) .or. (value2 .le. 0.d0) ) then
          rhoav(comp)=0.d0
        else
          tmp1=dv*value1
          if (tmp1<1.d0) then
            tmp1=(10.d0-15.d0*tmp1+6.d0*tmp1**2)*tmp1**3
          else
          tmp1=1.d0
          end if
          tmp2=dv*value2
          if (tmp2<1.d0) then
            tmp2=(10.d0-15.d0*tmp2+6.d0*tmp2**2)*tmp2**3
          else
            tmp2=1.d0
          end if
          rhoav(comp)=molmass(comp)*(value1+value2)/2.d0*tmp1*tmp2
        endif
      case default
        call bl_error("compute_nonnegative_rho_av: invalid avg_type")
      end select
    end do

  end subroutine compute_nonnegative_rho_av

  ! This routine must be called with non-negative densities, i.e., after calling compute_nonnegative_rho_av
  subroutine compute_sqrtLonsager_local(rho,rhotot,sqrtLonsager)
   
    real(kind=dp_t), intent(in)   :: rho(nspecies)            
    real(kind=dp_t), intent(in)   :: rhotot
    real(kind=dp_t), intent(out)  :: sqrtLonsager(nspecies,nspecies) 

    ! local variables
    integer                              :: row,column,info
    real(kind=dp_t), dimension(nspecies) :: W 
    real(kind=dp_t)                      :: rcond 

    real(kind=dp_t) :: molarconc(nspecies)
    real(kind=dp_t) :: molmtot
    real(kind=dp_t) :: chi(nspecies,nspecies)
    real(kind=dp_t) :: D_bar(nspecies,nspecies)

    real(kind=dp_t) :: rho_tmp(nspecies)
    real(kind=dp_t) :: rhotot_tmp

  
    type(bl_prof_timer), save :: bpt

    call build(bpt,"compute_sqrtLonsager_local")

    ! make copies of densities and adjust so they contain positive values
    rho_tmp(1:nspecies) = max(rho(1:nspecies), fraction_tolerance*rhotot)
    rhotot_tmp = sum(rho_tmp(1:nspecies))

    ! compute molarconc and molmtot
    call compute_molconc_molmtot_local(rho_tmp,rhotot_tmp,molarconc,molmtot)

    ! compute D_bar
    call compute_D_bar_local(rho_tmp,rhotot_tmp,D_bar)

    ! compute chi
    call compute_chi_local(rho_tmp,rhotot_tmp,molarconc,chi,D_bar)

    ! Restore original values (which have been made non-negative already by compute_nonnegative_rho_av)
    rho_tmp(1:nspecies) = rho(1:nspecies)
    rhotot_tmp = sum(rho_tmp(1:nspecies))

    ! compute massfraction W_i using these non-negative values *without* the correction of order fraction_tolerance
    ! This ensures that if W_i=0 there is no stochastic flux generated at all for species i
    do row=1, nspecies  
       W(row) = rho_tmp(row)/rhotot_tmp
    end do

    ! compute Onsager matrix L (store in sqrtLonsager)
    do column=1, nspecies
       do row=1, nspecies
          sqrtLonsager(row, column) = molmtot*rhotot*W(row)*chi(row,column)*W(column)/k_B
       end do
    end do

    ! compute cell-centered Cholesky factor, sqrtLonsager
    if(use_lapack) then
       
       call dpotrf_f95(sqrtLonsager,'L', rcond, 'I', info)
       !stop "LAPACK95 dpotrf_f95 disabled"
    
       ! remove all upper-triangular entries and NXN entry that lapack doesn't set to zero 
       do row=1, nspecies
          do column=row+1, nspecies
             sqrtLonsager(row, column) = 0.0d0          
          end do
       end do    
       sqrtLonsager(nspecies, nspecies) = 0.0d0          
    
    else
       call choldc(sqrtLonsager,nspecies)   
    end if

    call destroy(bpt)

  end subroutine compute_sqrtLonsager_local

  subroutine compute_rhoWchi(mla,rho,chi,rhoWchi)
 
    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: rho(:)
    type(multifab) , intent(in   )  :: chi(:) 
    type(multifab) , intent(inout)  :: rhoWchi(:) 

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: n,i,dm,nlevs,ng_1,ng_3,ng_4
 
    ! pointer for rho(nspecies), molarconc(nspecies) 
    real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for rho    
    real(kind=dp_t), pointer        :: dp3(:,:,:,:)  ! for chi
    real(kind=dp_t), pointer        :: dp4(:,:,:,:)  ! for rhoWchi

    type(mfiter) :: mfi
    type(box) :: tilebox
    integer :: tlo(mla%dim), thi(mla%dim)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"compute_rhoWchi")

    dm = mla%dim        ! dimensionality
    ng_1 = rho(1)%ng    ! number of ghost cells 
    ng_3 = chi(1)%ng
    ng_4 = rhoWchi(1)%ng
    nlevs = mla%nlevel  ! number of levels 
 
    !$omp parallel private(n,i,mfi,tilebox,tlo,thi) &
    !$omp private(dp1,dp3,dp4,lo,hi)

    ! loop over all boxes 
    do n=1,nlevs
       call mfiter_build(mfi, rho(n), tiling=.true.)

       do while (more_tile(mfi))
          i = get_fab_index(mfi)

          tilebox = get_growntilebox(mfi,rhoWchi(n)%ng)
          tlo = lwb(tilebox)
          thi = upb(tilebox)

!       do i=1,nfabs(rho(n))
          dp1 => dataptr(rho(n), i)
          dp3 => dataptr(chi(n), i)
          dp4 => dataptr(rhoWchi(n), i)
          lo  =  lwb(get_box(rho(n), i))
          hi  =  upb(get_box(rho(n), i))
          
          select case(dm)
          case (2)
             call compute_rhoWchi_2d(dp1(:,:,1,:),dp3(:,:,1,:),dp4(:,:,1,:), &
                                     ng_1,ng_3,ng_4,lo,hi,tlo,thi) 
          case (3)
             call compute_rhoWchi_3d(dp1(:,:,:,:),dp3(:,:,:,:),dp4(:,:,:,:), &
                                     ng_1,ng_3,ng_4,lo,hi,tlo,thi) 
          end select
       end do
    end do
    !$omp end parallel

    call destroy(bpt)

  end subroutine compute_rhoWchi
  
  subroutine compute_rhoWchi_2d(rho,chi,rhoWchi,ng_1,ng_3,ng_4,glo,ghi,tlo,thi)
  
    integer          :: glo(2), ghi(2), ng_1,ng_3,ng_4,tlo(2),thi(2)
    real(kind=dp_t)  ::     rho(glo(1)-ng_1:,glo(2)-ng_1:,:) ! density; last dimension for species
    real(kind=dp_t)  ::     chi(glo(1)-ng_3:,glo(2)-ng_3:,:) ! last dimension for nspecies^2
    real(kind=dp_t)  :: rhoWchi(glo(1)-ng_4:,glo(2)-ng_4:,:) ! last dimension for nspecies^2

    ! local variables
    integer          :: i,j
  
    ! for specific box, now start loops over alloted cells 
    do j=tlo(2),thi(2)
       do i=tlo(1),thi(1)
        
          call compute_rhoWchi_local(rho(i,j,:),chi(i,j,:),rhoWchi(i,j,:))

       end do
    end do

  end subroutine compute_rhoWchi_2d

  subroutine compute_rhoWchi_3d(rho,chi,rhoWchi,ng_1,ng_3,ng_4,glo,ghi,tlo,thi)

    integer          :: glo(3), ghi(3), ng_1,ng_3,ng_4,tlo(3),thi(3)
    real(kind=dp_t)  ::     rho(glo(1)-ng_1:,glo(2)-ng_1:,glo(3)-ng_1:,:) ! density; last dimension for species
    real(kind=dp_t)  ::     chi(glo(1)-ng_3:,glo(2)-ng_3:,glo(3)-ng_3:,:) ! last dimension for nspecies^2
    real(kind=dp_t)  :: rhoWchi(glo(1)-ng_4:,glo(2)-ng_4:,glo(3)-ng_4:,:) ! last dimension for nspecies^2
    
    ! local variables
    integer          :: i,j,k

    ! for specific box, now start loops over alloted cells 
    do k=tlo(3),thi(3)
       do j=tlo(2),thi(2)
          do i=tlo(1),thi(1)
       
             call compute_rhoWchi_local(rho(i,j,k,:),chi(i,j,k,:),rhoWchi(i,j,k,:))
              
         end do
      end do
    end do
   
  end subroutine compute_rhoWchi_3d
  
  subroutine compute_rhoWchi_local(rho,chi,rhoWchi)
   
    real(kind=dp_t), intent(in)   :: rho(nspecies)            
    real(kind=dp_t), intent(in)   :: chi(nspecies,nspecies)   ! rank conversion done 
    real(kind=dp_t), intent(out)  :: rhoWchi(nspecies,nspecies) 
 
    ! local variables
    integer                              :: row,column

    ! populate rho*W*chi = rho_i*chi
    do row=1, nspecies
       do column=1, nspecies
          rhoWchi(row,column) = rho(row)*chi(row,column)  
       end do
    end do

  end subroutine compute_rhoWchi_local

  subroutine compute_baro_coef(mla,baro_coef,rho,rhotot,Temp)
 
    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(inout)  :: baro_coef(:)
    type(multifab) , intent(in   )  :: rho(:)
    type(multifab) , intent(in   )  :: rhotot(:) 
    type(multifab) , intent(in   )  :: Temp(:) 

    ! local
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: n,i,dm,nlevs,ng_1,ng_2,ng_3,ng_4

    real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for baro_coef
    real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for rho
    real(kind=dp_t), pointer        :: dp3(:,:,:,:)  ! for rhotot
    real(kind=dp_t), pointer        :: dp4(:,:,:,:)  ! for Temp

    type(bl_prof_timer), save :: bpt

    call build(bpt,"compute_baro_coef")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_1 = baro_coef(1)%ng
    ng_2 = rho(1)%ng
    ng_3 = rhotot(1)%ng
    ng_4 = Temp(1)%ng

    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp1 => dataptr(baro_coef(n),i)
          dp2 => dataptr(rho(n),i)
          dp3 => dataptr(rhotot(n),i)
          dp4 => dataptr(Temp(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          select case(dm)
          case (2)
             call compute_baro_coef_2d(dp1(:,:,1,:),dp2(:,:,1,:),dp3(:,:,1,1),dp4(:,:,1,1), &
                                        ng_1,ng_2,ng_3,ng_4,lo,hi) 
          case (3)
             call bl_error("compute_baro_coef_3d not written yet")
          end select
       end do
    end do

    call destroy(bpt)

  end subroutine compute_baro_coef

  subroutine compute_baro_coef_2d(baro_coef,rho,rhotot,Temp,ng_1,ng_2,ng_3,ng_4,lo,hi)
 
    integer        , intent(in   ) :: lo(2), hi(2), ng_1, ng_2, ng_3, ng_4
    real(kind=dp_t), intent(inout) :: baro_coef(lo(1)-ng_1:,lo(2)-ng_1:,:)
    real(kind=dp_t), intent(in   ) ::       rho(lo(1)-ng_2:,lo(2)-ng_2:,:)
    real(kind=dp_t), intent(in   ) ::    rhotot(lo(1)-ng_3:,lo(2)-ng_3:)
    real(kind=dp_t), intent(in   ) ::      Temp(lo(1)-ng_4:,lo(2)-ng_4:)
    
    ! local variables
    integer :: i,j,comp
    real(kind=dp_t) :: n

    do j=lo(2)-ng_1,hi(2)+ng_1
       do i=lo(1)-ng_1,hi(1)+ng_1

          n = 0.d0
          do comp=1,nspecies
             n = n + rho(i,j,comp)/molmass(comp)
          end do

          do comp=1,nspecies
             baro_coef(i,j,comp) = rho(i,j,comp)/rhobar(comp) / (n*k_B*Temp(i,j))
          end do

       end do
    end do

  end subroutine compute_baro_coef_2d
  
end module mass_flux_utilities_module
