module init_module

  use multifab_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use probin_common_module
  use probin_multispecies_module
  use F95_LAPACK 
 
  implicit none

  private

  public :: init_rho

contains
  
  subroutine init_rho(rho,Dbar,Gama,mass,dx,prob_lo,prob_hi,the_bc_level,rho_part_bc_comp)

    type(multifab) , intent(inout) :: rho(:)            
    real(kind=dp_t), intent(inout) :: Dbar(:,:)           
    real(kind=dp_t), intent(inout) :: Gama(:,:)           
    real(kind=dp_t), intent(inout) :: mass(:)           
    real(kind=dp_t), intent(in   ) :: dx(:,:)           
    real(kind=dp_t), intent(in   ) :: prob_lo(rho(1)%dim)
    real(kind=dp_t), intent(in   ) :: prob_hi(rho(1)%dim)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    integer,         intent(in   ) :: rho_part_bc_comp
 
    ! local variables
    integer :: lo(rho(1)%dim), hi(rho(1)%dim)
    integer :: dm, ng, i, n, nlevs
    real(kind=dp_t), pointer :: dp(:,:,:,:)   ! pointer for rho (last dim: nspecies)   

    dm = rho(1)%dim
    ng = rho(1)%ng
    nlevs = size(rho,1)

    ! populate SM Dbar matrix, Gama, masses etc
    call populate_DbarGama(Dbar(:,:), Gama(:,:), mass(:)) 

    ! looping over boxes 
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp => dataptr(rho(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          !print*, lo, hi 
          
          select case(dm)
          case (2)
             call init_rho_2d(dp(:,:,1,:), ng, lo, hi, prob_lo, prob_hi, dx(n,:))
          !case (3)
          !   call init_rho_3d(dp(:,:,:,:), dp1(:,:,:,:), ng, lo, hi, prob_lo, prob_hi, dx(n,:))
          end select
       end do

       ! filling up ghost cells for two adjacent grids at the same level
       ! including periodic domain boundary ghost cells
       call multifab_fill_boundary(rho(n))

       ! fill non-periodic domain boundary ghost cells
       ! Amit: changed here.
       call multifab_physbc(rho(n),1,rho_part_bc_comp,nspecies,the_bc_level(n),dx(n,:),.false.)
    end do

  end subroutine init_rho

  subroutine populate_DbarGama(Dbar, Gama, mass)
  
    real(kind=dp_t)  :: Dbar(:,:)   ! SM diffusion constants 
    real(kind=dp_t)  :: Gama(:,:)   ! non-ideality coefficient 
    real(kind=dp_t)  :: mass(:)     ! mass of each species
    integer          :: n,row,column
 
    ! populate Dbar,Gama & mass, which for initial case doesn't change in each cell. 
    n=0; 
    do row=1, nspecies  
       do column=1, row-1
          n=n+1
          Dbar(row, column) = Dbar_bc(n)
          Dbar(column, row) = Dbar(row, column) ! symmetric
          Gama(row, column) = 0.d0       
          Gama(column, row) = Gama(row, column) ! symmetric
       enddo
       Dbar(row, row) = 0.d0     ! self-diffusion is zero
       Gama(row, row) = 1.d0     ! set to unit matrix for time being
       mass(row)      = m_bc(row) ! set to unit matrix for time being
       !mtot = mtot + m_bc(row)   ! total mass of the species
    enddo

  end subroutine populate_DbarGama

  subroutine init_rho_2d(rho, ng, lo, hi, prob_lo, prob_hi, dx)

    integer          :: lo(2), hi(2), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)  ! last dimension for species
    real(kind=dp_t)  :: prob_lo(2)
    real(kind=dp_t)  :: prob_hi(2)
    real(kind=dp_t)  :: dx(:)
 
    ! local varables
    integer          :: i,j
    real(kind=dp_t)  :: x,y,L(2)
 
    L(1:2) = prob_hi(1:2)-prob_lo(1:2) ! Domain length
    
    ! for specific box, now start loop over alloted cells     
    
    ! select problem type, 1=bubble, 2=constant gradient
    select case(init_type) 
    case(1) 
      do j=lo(2),hi(2)
         y = prob_lo(2) + (dble(j)+0.5d0) * dx(2) - 0.5d0
         do i=lo(1),hi(1)
            x = prob_lo(1) + (dble(i)+0.5d0) * dx(1) - 0.5d0
       
            if ((x-L(1)*0.5d0)**2 + (y-L(2)*0.5d0)**2 & 
                 .lt. L(1)*L(2)*0.1d0) then
                 rho(i,j,1:nspecies) = c_bc(1,1:nspecies)
            else
                 rho(i,j,1:nspecies) = c_bc(2,1:nspecies)
            endif
    
         end do
      end do
  
    case(2) 
      do j=lo(2),hi(2)
         y = prob_lo(2) + (dble(j)+0.5d0) * dx(2) - 0.5d0
         do i=lo(1),hi(1)
            x = prob_lo(1) + (dble(i)+0.5d0) * dx(1) - 0.5d0
        
            rho(i,j,1:nspecies) = sin(x*y)
            !if(i.eq.4 .and. j.eq.4)
              ! print*, x, y, rho(i,j,1:nspecies)
            !endif
         end do
            print*, ''
      end do
    
    end select
   
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
    real(kind=dp_t)  :: x,y,z,L(3)

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
                 rho(i,j,k,1:nspecies) = c_bc(1,1:nspecies)
             else
                 rho(i,j,k,1:nspecies) = c_bc(2,1:nspecies)
             endif
          
          end do
       end do
    end do
    !$omp end parallel do
 
  end subroutine init_rho_3d

end module init_module
