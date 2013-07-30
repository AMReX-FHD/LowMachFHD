module init_module

  use multifab_module
  use define_bc_module
  use multifab_physbc_module
  use probin_common_module
  use probin_multispecies_module
  use F95_LAPACK 
 
  implicit none

  private

  public :: init_rho

contains
  
  !subroutine init_rho(rho,molarconc,BinvGama,dx,prob_lo,prob_hi,the_bc_level)
  subroutine init_rho(rho,BinvGama,dx,prob_lo,prob_hi,the_bc_level)

    type(multifab) , intent(inout) :: rho(:)            
    !type(multifab) , intent(inout) :: molarconc(:)            
    type(multifab) , intent(inout) :: BinvGama(:)   
    real(kind=dp_t), intent(in   ) :: dx(:,:)           
    real(kind=dp_t), intent(in   ) :: prob_lo(rho(1)%dim)
    real(kind=dp_t), intent(in   ) :: prob_hi(rho(1)%dim)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local variables
    integer :: lo(rho(1)%dim), hi(rho(1)%dim)
    integer :: nspecies, dm, ng, i, n, nlevs

    ! pointer for rho(nspecies), molarconc(nspecies), BinvGama(nspecies^2)
    real(kind=dp_t), pointer :: dp(:,:,:,:)       
    !real(kind=dp_t), pointer :: dp1(:,:,:,:)       
    real(kind=dp_t), pointer :: dp2(:,:,:,:) 
 
    dm = rho(1)%dim
    ng = rho(1)%ng
    nlevs = size(rho,1)

    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp => dataptr(rho(n),i)
          !dp1 => dataptr(molarconc(n),i)
          dp2 => dataptr(BinvGama(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          
          select case(dm)
          case (2)
             !call init_rho_2d(dp(:,:,1,:), dp1(:,:,1,:), dp2(:,:,1,:), ng, lo, hi, prob_lo, prob_hi, dx(n,:)) 
             call init_rho_2d(dp(:,:,1,:), dp2(:,:,1,:), ng, lo, hi, prob_lo, prob_hi, dx(n,:)) 
          !case (3)
          !   call init_rho_3d(dp(:,:,:,:), dp1(:,:,:,:), ng, lo, hi, prob_lo, prob_hi, dx(n,:))
          end select
       end do

       ! filling up ghost cells for two adjacent grids at the same level
       ! including periodic domain boundary ghost cells
       call multifab_fill_boundary(rho(n))
       !call multifab_fill_boundary(molarconc(n))
       call multifab_fill_boundary(BinvGama(n))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(rho(n),1,1,nspecies,the_bc_level(n))
       !call multifab_physbc(molarconc(n),1,1,nspecies,the_bc_level(n))
       call multifab_physbc(BinvGama(n),1,1,nspecies**2,the_bc_level(n))
    end do

  end subroutine init_rho

  !subroutine init_rho_2d(rho, molarconc, BinvGama, ng, lo, hi, prob_lo, prob_hi, dx)
  subroutine init_rho_2d(rho, BinvGama, ng, lo, hi, prob_lo, prob_hi, dx)

    integer          :: lo(2), hi(2), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)      ! Last dimension for species
    !real(kind=dp_t)  :: molarconc(lo(1)-ng:,lo(2)-ng:,:)      ! Last dimension for species
    real(kind=dp_t)  :: BinvGama(lo(1)-ng:,lo(2)-ng:,:) ! B^(-1)*Gamma nspecies X nspecies
    real(kind=dp_t)  :: prob_lo(2)
    real(kind=dp_t)  :: prob_hi(2)
    real(kind=dp_t)  :: dx(:)

    ! local varables
    integer          :: i,j,k,n
    integer          :: row,column
    real(kind=dp_t)  :: x,y,m,rho_cell,alpha,L(2)
 
    !for LAPACK 
    real(kind=dp_t), dimension(nspecies,nspecies) :: Bij,Dbar,Gama,c

    L(1:2) = prob_hi(1:2)-prob_lo(1:2) ! Domain length
    
    !populate local quantities: SM diff_coeffs, Gamma which for initial case
    !doesn't change in each and every cell. I don't need to populate D(n,n)
    !which is zero and not used anywhere. Also calculate total mass.  
    n=0; 
    m = 0.d0
    do row=1, nspecies  
       do column=1, row-1
          n=n+1
          Dbar(row, column) = Dbar_bc(n)
          Dbar(column, row) = Dbar(row, column) ! symmetric
          Gama(row, column) = 0.d0       
          Gama(column, row) = Gama(row, column) ! symmetric
       enddo
       Dbar(row, row) = 0.d0   ! set to unit matrix
       Gama(row, row) = 1.d0   ! set to unit matrix
       m = m + m_bc(row)
    enddo
     
    ! Start loops over each cell     
    do j=lo(2),hi(2)
       y = prob_lo(2) + (dble(j)+0.5d0) * dx(2) - 0.5d0
       do i=lo(1),hi(1)
          x = prob_lo(1) + (dble(i)+0.5d0) * dx(1) - 0.5d0
          if ((x-L(1)*0.5d0)**2 + (y-L(2)*0.5d0)**2 & 
               .lt. L(1)*L(2)*0.1d0) then
             rho(i,j,1)           = c_bc(1,1)
             rho(i,j,2:nspecies)  = c_bc(1,2) 
            
             rho_cell=0.d0 
             do row=1, nspecies  ! calculate total density inside cell
                rho_cell = rho_cell + rho(i,j,row)
             enddo         
            
             ! calculate molar concentrations 
             !do row=1, nspecies 
             !   molarconc(i,j,row) = m*rho(i,j,row)/(m_bc(row)*rho_cell)
             !enddo
    
             !!!!!!!!!! Calculate Bprimeij, then Bij and then Bij^(-1)*Gamma !!!!!!!!!!!
             !!!!!!!!!! May be one more subroutine to make as it's called twice!!!!!!!!!!!
             Bij=0.d0
             do row=1, nspecies  
                do column=1, row-1
                   Bij(row, column) = rho(i,j,row)*m**2/(m_bc(row)* &
                            m_bc(column)*Dbar(row, column)*rho_cell**2) 
                   Bij(column, row) = rho(i,j,column)*m**2/(m_bc(row)* &
                            m_bc(column)*Dbar(column, row)*rho_cell**2) 
                enddo

                
                do column=1, nspecies
                   if (column.ne.row) then
                     Bij(row, row) = Bij(row, row) - m**2/(m_bc(row)*rho_cell**2)* & 
                                        (rho(i,j,column)/(m_bc(column)*Dbar(row,column)))
                   endif
                enddo
             enddo

             if(.false.)  then
             if(i.eq.lo(1) .and. j.eq.lo(2)) then  ! checking B matrix for one cell            
                do row=1, nspecies
                   do column=1, nspecies 
                      print*, Bij(row, column)
                   enddo
                   print*, ''        
                enddo
             endif 
             endif

             alpha = maxval(Bij) ! adjust parameter alpha to max val of Bij
  
             do row=1, nspecies   !transform Bprimeij to Bij
                Bij(row, row) = alpha + Bij(row, row)
             enddo
    
             call la_gesvx(A=Bij, B=Gama, X=c) ! calculate A^(-1)*B; result is written to B 
             
             !store B^(-1)*Gamma on each cell via Lexicographic projection
             do row=1, nspecies
                do column=1, nspecies 
                   BinvGama(i,j,row+nspecies*(column-1)) = Bij(row, column) 
                   if(.false.) then 
                      if(i.eq.lo(1) .and. j.eq.lo(2)) then  ! checking B matrix for one cell            
                         print*, BinvGama(i,j,row+nspecies*(column-1))
                      endif
                   endif
                enddo
             enddo
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 
          else
             rho(i,j,1)           = c_bc(2,1)
             rho(i,j,2:nspecies)  = c_bc(2,2)
          
             ! Amit: non-understandable how to deal with zero densities
             !!!!!!!!!! Calculate Bprimeij, then Bij and then Bij^(-1)*Gamma !!!!!!!!!!!
             Bij=0.d0
             do row=1, nspecies  
                do column=1, row-1
                   Bij(row, column) = rho(i,j,row)*m**2/(m_bc(row)* &
                            m_bc(column)*Dbar(row, column)*rho_cell**2) 
                   Bij(column, row) = rho(i,j,column)*m**2/(m_bc(row)* &
                            m_bc(column)*Dbar(column, row)*rho_cell**2) 
                enddo

                
                do column=1, nspecies
                   if (column.ne.row) then
                     Bij(row, row) = Bij(row, row) - m**2/(m_bc(row)*rho_cell**2)* & 
                                        (rho(i,j,column)/(m_bc(column)*Dbar(row,column)))
                   endif
                enddo
             enddo

             alpha = maxval(Bij) ! adjust parameter alpha to max val of Bij
  
             do row=1, nspecies   !transform Bprimeij to Bij
                Bij(row, row) = alpha + Bij(row, row)
             enddo
    
             call la_gesvx(A=Bij, B=Gama, X=c) ! calculate A^(-1)*B; result is written to B 
             
             !store B^(-1)*Gamma on each cell via Lexicographic projection
             do row=1, nspecies
                do column=1, nspecies 
                   BinvGama(i,j,row+nspecies*(column-1)) = Bij(row, column) 
                enddo
             enddo
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
                 
          endif
       end do
    end do
 
    end subroutine init_rho_2d

    subroutine init_rho_3d(rho, diff_coeffs, ng, lo, hi, prob_lo, prob_hi, dx)

    integer         :: lo(3), hi(3), ng
    real(kind=dp_t) :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:) ! Last dimension for species 
    real(kind=dp_t) :: diff_coeffs(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:) 
    real(kind=dp_t) :: prob_lo(3)
    real(kind=dp_t) :: prob_hi(3)
    real(kind=dp_t) :: dx(:)
 
    ! local variables
    integer          :: i,j,k
    real(kind=dp_t) :: x,y,z,L(3)

    L(1:3) = prob_hi(1:3)-prob_lo(1:3) ! Domain length
    
    !$omp parallel private(i,j,k,x,y,z)
    do k=lo(3),hi(3)
       z = prob_lo(3) + (dble(k)+0.5d0) * dx(3) - 0.5d0
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0) * dx(2) - 0.5d0
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0) * dx(1) - 0.5d0
             if ((x-L(1)*0.5d0)**2 + (y-L(2)*0.5d0)**2 + & 
                 (z-L(3)*0.5d0)**2 .lt. L(1)*L(2)*L(3)*0.001d0) then
                 rho(i,j,k,1)           = c_bc(1,1)
                 rho(i,j,k,2:nspecies)  = c_bc(1,2)
                 diff_coeffs(i,j,k,1:nspecies) = d_bc(1:nspecies)
             else
                 rho(i,j,k,1)           = c_bc(2,1)
                 rho(i,j,k,2:nspecies)  = c_bc(2,2)
                 diff_coeffs(i,j,k,1:nspecies) = d_bc(1:nspecies)
             endif
          end do
       end do
    end do
    !$omp end parallel do
 
  end subroutine init_rho_3d

end module init_module
