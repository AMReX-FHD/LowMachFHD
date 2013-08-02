module diffusive_flux_module

  use multifab_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use div_and_grad_module
  use probin_multispecies_module
  use ml_layout_module
  use convert_stag_module
  use F95_LAPACK
  
  implicit none

  private

  public :: diffusive_flux

contains
 
  subroutine diffusive_flux(mla,rho,molarconc,BinvGama,Dbar,Gama,flux,dx,the_bc_level) 

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rho(:)
    type(multifab) , intent(inout) :: molarconc(:) ! Donev: Make this intent(in)
    type(multifab) , intent(inout) :: BinvGama(:) ! Cell-centered coefficients
    real(kind=dp_t), intent(in   ) :: Dbar(:,:)
    real(kind=dp_t), intent(in   ) :: Gama(:,:)
    type(multifab) , intent(inout) :: flux(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local variables
    integer :: lo(rho(1)%dim), hi(rho(1)%dim)
    integer :: n,i,ng,dm,nlevs
    
    ! for the face centered B^(-1)*Gama
    type(multifab) :: BinvGama_face(mla%nlevel,mla%dim)
    
    ! pointer for rho(nspecies), molarconc(nspecies), BinvGama(nspecies^2)
    real(kind=dp_t), pointer       :: dp(:,:,:,:)   ! for rho    
    real(kind=dp_t), pointer       :: dp1(:,:,:,:)  ! for mole x
    real(kind=dp_t), pointer       :: dp2(:,:,:,:)  ! for B^(-1)*Gama

    dm    = mla%dim     ! dimensionality
    ng = rho(1)%ng      ! number of ghost cells 
    nlevs = mla%nlevel  ! number of levels 

    ! loop over all boxes 
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp => dataptr(rho(n),i)
          dp1 => dataptr(molarconc(n),i)
          dp2 => dataptr(BinvGama(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          
          select case(dm)
          case (2)
             call cal_rhoxBinvGama_2d(dp(:,:,1,:), dp1(:,:,1,:), dp2(:,:,1,:), & 
                           Dbar(:,:), Gama(:,:), ng, lo, hi, dx(n,:)) 
          !case (3)
          !   call init_rho_3d(dp(:,:,:,:), dp1(:,:,:,:), ng, lo, hi, prob_lo, prob_hi, dx(n,:))
          end select
       end do

       ! filling up ghost cells for two adjacent grids at the same level
       ! including periodic domain boundary ghost cells
       call multifab_fill_boundary(rho(n))
       call multifab_fill_boundary(molarconc(n))
       call multifab_fill_boundary(BinvGama(n))

       ! fill non-periodic domain boundary ghost cells 
       ! Amit: shouldn't the 2nd 1 be scal_bc_comp? And also nspecies^2 in BinvGama?
       call multifab_physbc(rho(n),1,1,nspecies,the_bc_level(n))
       call multifab_physbc(molarconc(n),1,1,nspecies,the_bc_level(n))
       ! Donev: BinvGamma should NOT need values behind a boundary
       ! so don't call multifab_physbc here --- focus on periodic BCs first
       call multifab_physbc(BinvGama(n),1,1,nspecies,the_bc_level(n))
    end do
    
    ! calculate face centrered grad(molarconc) 
    call compute_grad(mla,molarconc,flux,dx,1,scal_bc_comp,1,nspecies,the_bc_level)
   
    ! change B^(-1)*Gama from cell-centered to face-centered
    ! Donev: These have nspecies**2 components, not nspecies!!!
    ! Should be done in a loop
    ! do i=1, nspecies**2
    ! call average_cc_to_face(nlevs,BinvGama,BinvGama_face,i,scal_bc_comp,1,the_bc_level) 
    ! Or, if you add +1 component to the bc_tower, you do:
    ! call average_cc_to_face(nlevs,BinvGama,BinvGama_face,i,scal_bc_comp+nspecies,1,the_bc_level)     
    ! end do
    call average_cc_to_face(nlevs,BinvGama,BinvGama_face,1,scal_bc_comp,nspecies,the_bc_level) 

    ! compute flux as B^(-1)*Gama X grad(molarconc)
    do n=1,nlevs
       do i=1,dm
          ! Donev: NO!!! This needs to do a matrix-vector product
          ! The routine below only muliplies component-by-component
          ! two multifabs that have the same size
          ! You need to debug the output to make sure what you are computing here
          call multifab_mult_mult_c(flux(n,i),1,BinvGama_face(n,i),1,nspecies,0)          
       end do
    end do
    
    !Donev: If grad(temperature)
    !call compute_grad(mla,temperature,flux,dx,1,scal_bc_comp+n_species,n_species+1,1,the_bc_level)
 
    ! destroy B^(-1)*Gama multifab to prevent leakage in memory
    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(BinvGama_face(n,i))
       end do
    end do

  end subroutine diffusive_flux

  
  subroutine cal_rhoxBinvGama_2d(rho, molarconc, BinvGama, Dbar, Gama, ng, lo, hi, dx)

    integer          :: lo(2), hi(2), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)       ! density; last dimension for species
    real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,:) ! molar concentration; Donev: intent(in)
    ! Donev: As I explained in my notes, you can do here:    
    !real(kind=dp_t)  :: BinvGama(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,1:nspecies,1:nspecies) 
    ! This makes BinvGama a rank-4 array, see below
    ! Or, even better, see my comment at the end of this routine
    real(kind=dp_t)  :: BinvGama(lo(1)-ng:,lo(2)-ng:,:)  ! B^(-1)*Gamma; last dimension for nspecies^2
    real(kind=dp_t)  :: Dbar(:,:)                        ! SM diffusion constants 
    real(kind=dp_t)  :: Gama(:,:)                        ! non-ideality coefficient 
    real(kind=dp_t)  :: dx(:)                            ! grid spacing

    ! local variables
    integer          :: i,j,k,n
    integer          :: row,column
    real(kind=dp_t)  :: mtot,rho_cell,alpha
 
    ! dummy matrices for inversion using LAPACK 
    real(kind=dp_t), dimension(nspecies,nspecies) :: Bij,c

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!Calculate BinvGama. This chunk should be copied in 3D code too.!!!! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Amit: mtot doesn't change(?), so we can declare it as a global
    ! variable in the probin_multispecies code. 
    mtot = 0.d0
    do row=1, nspecies  
       mtot = mtot + m_bc(row)
    enddo
  
    ! for specific box, now start loops over alloted cells    
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
         
          ! calculate total density inside one-cell 
          ! Donev: This should be moved to convert_cons_to_prim
          do row=1, nspecies 
          rho_cell=0.d0 
          do row=1, nspecies  
             rho_cell = rho_cell + rho(i,j,row)
          enddo         
             
          ! calculate molar concentrations 
          ! Donev: This should be moved to convert_cons_to_prim
          do row=1, nspecies 
             molarconc(i,j,row) = mtot*rho(i,j,row)/(m_bc(row)*rho_cell)
          enddo
   
          ! calculate Bprime matrix (stored as Bij to save memory)
          Bij=0.d0
          do row=1, nspecies  
             do column=1, row-1
                Bij(row, column) = rho(i,j,row)*mtot**2/(m_bc(row)* &
                       m_bc(column)*Dbar(row, column)*rho_cell**2) 
                Bij(column, row) = rho(i,j,column)*mtot**2/(m_bc(row)* &
                       m_bc(column)*Dbar(column, row)*rho_cell**2) 
             enddo
             
             do column=1, nspecies
                if (column.ne.row) then
                   Bij(row, row) = Bij(row, row) - mtot**2/(m_bc(row)*rho_cell**2)* & 
                              (rho(i,j,column)/(m_bc(column)*Dbar(row,column)))
                endif
             enddo
          enddo

          if(.false.)  then
            if(i.eq.lo(1) .and. j.eq.lo(2)) then  ! optional checkng Bij, Dbar etc            
               do row=1, nspecies
                  do column=1, nspecies 
                     print*, Bij(row, column)
                     print*, Dbar(row, column)
                  enddo
                  print*, ''        
               enddo
            endif 
          endif

          ! adjust parameter alpha to max val of Bij
          alpha = maxval(Bij) ! Donev: this needs maxval(abs(Bij))
  
          ! transform Bprimeij to Bij
          do row=1, nspecies   
             Bij(row, row) = alpha + Bij(row, row)
          enddo
    
          ! calculate A^(-1)*B; result is written to second argument (B) 
          call la_gesvx(A=Bij, B=Gama, X=c) 
             
          ! Donev:
          ! If you do what I said in the comment above, you can just do:
          ! BinvGama(i,j,:,:)=Bij
          ! and not need to do the stuff below   
          ! Or, even better, see comment at the end and do this:
          ! call set_Bij(BinvGama(i,j,:),Bij)
          
          
          ! store B^(-1)*Gamma on each cell via Lexicographic order
          do row=1, nspecies
             do column=1, nspecies 
                BinvGama(i,j,column+nspecies*(row-1)) = Bij(row, column) 
                if(.false.) then 
                   if(i.eq.lo(1) .and. j.eq.lo(2)) then  ! optional check BinvGama for a cell            
                      print*, BinvGama(i,j,row+nspecies*(column-1))
                   endif
                endif
             enddo
          enddo
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
       end do
    end do
    
    ! Donev:
    ! Use contained (internal) subroutine to do the copy without doing index algebra:
    ! contains
    ! subroutine set_Bij(BinvGama_ij, B_ij)
    !    real(kind=dp_t), dimension(nspecies,nspecies), intent(in) :: B_ij
    !    real(kind=dp_t), dimension(nspecies,nspecies), intent(out) :: BinvGamma_ij
    !    BinvGamma_ij = B_ij
    ! end subroutine
 
    end subroutine cal_rhoxBinvGama_2d

end module diffusive_flux_module

