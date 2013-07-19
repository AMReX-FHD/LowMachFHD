module init_module

  use multifab_module
  use define_bc_module
  use multifab_physbc_module
  use probin_common_module
  use probin_multispecies_module
  
  implicit none

  private

  public :: init_rho

contains
  
  subroutine init_rho(rho,diff_coeffs,dx,prob_lo,prob_hi,the_bc_level)

    type(multifab) , intent(inout) :: rho(:)            
    type(multifab) , intent(inout) :: diff_coeffs(:)   
    real(kind=dp_t), intent(in   ) :: dx(:,:)           
    real(kind=dp_t), intent(in   ) :: prob_lo(rho(1)%dim)
    real(kind=dp_t), intent(in   ) :: prob_hi(rho(1)%dim)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local variables
    integer :: lo(rho(1)%dim), hi(rho(1)%dim)
    integer :: nspecies, dm, ng, i, n, nlevs

    ! pointer for rho & diff_coeffs, last index for nspecies
    real(kind=dp_t), pointer :: dp(:,:,:,:)       
    real(kind=dp_t), pointer :: dp1(:,:,:,:) 
 
    dm = rho(1)%dim
    ng = rho(1)%ng
    nlevs = size(rho,1)

    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp => dataptr(rho(n),i)
          dp1 => dataptr(diff_coeffs(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          
          select case(dm)
          case (2)
             call init_rho_2d(dp(:,:,1,:), dp1(:,:,1,:), ng, lo, hi, prob_lo, prob_hi, dx(n,:)) 
          case (3)
             call init_rho_3d(dp(:,:,:,:), ng, lo, hi, prob_lo, prob_hi, dx(n,:))
          end select
       end do

       ! filling up ghost cells for two adjacent grids at the same level
       ! including periodic domain boundary ghost cells
       call multifab_fill_boundary(rho(n))
       call multifab_fill_boundary(diff_coeffs(n))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(rho(n),1,1,nspecies,the_bc_level(n))
       call multifab_physbc(diff_coeffs(n),1,1,nspecies,the_bc_level(n))
    end do

  end subroutine init_rho

  subroutine init_rho_2d(rho, diff_coeffs, ng, lo, hi, prob_lo, prob_hi, dx)

    integer          :: lo(2), hi(2), ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)         ! Last dimension for species
    real(kind=dp_t)  :: diff_coeffs(lo(1)-ng:,lo(2)-ng:,:) !to hold D matrix
    real(kind=dp_t)  :: prob_lo(2)
    real(kind=dp_t)  :: prob_hi(2)
    real(kind=dp_t)  :: dx(:)
 
    ! local varables
    integer          :: i,j
    real(kind=dp_t)  :: x,y,L(2)
    
    L(1:2) = prob_hi(1:2)-prob_lo(1:2) ! Domain length
 
    do j=lo(2),hi(2)
       y = prob_lo(2) + (dble(j)+0.5d0) * dx(2) - 0.5d0
       do i=lo(1),hi(1)
          x = prob_lo(1) + (dble(i)+0.5d0) * dx(1) - 0.5d0
          if ((x-L(1)*0.5d0)*(x-L(1)*0.5d0) + (y-L(2)*0.5d0)*(y-L(2)*0.5d0) & 
               .lt. L(1)*L(2)*0.1d0) then
             rho(i,j,1)           = c_bc(1,1)
             rho(i,j,2:nspecies)  = c_bc(1,2)
             diff_coeffs(i,j,1:nspecies) = d_bc(1:nspecies)
          else
             rho(i,j,1)           = c_bc(2,1)
             rho(i,j,2:nspecies)  = c_bc(2,2)
             diff_coeffs(i,j,1:nspecies) = d_bc(1:nspecies)
          endif 
       end do
    end do
 
    end subroutine init_rho_2d

    subroutine init_rho_3d(rho, ng, lo, hi, prob_lo, prob_hi, dx)

    integer          :: lo(3), hi(3), ng
    real(kind=dp_t) :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    !real(kind=dp_t) :: diff_coeffs(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:,:)
    ! Example:
    ! real(dp_t) :: diff_coeffs(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:,:) ! Last dimensions are (1:nspecies,1:nspecies)
    ! You can use D = diff_coeffs(i,j,k,1:nspecies,1:nspecies)
    real(kind=dp_t) :: prob_lo(3)
    real(kind=dp_t) :: prob_hi(3)
    real(kind=dp_t) :: dx(:)
 
    ! local variables
    integer          :: i,j,k
    real(kind=dp_t) :: x,y,z,L(3)

    L(1:3) = prob_hi(1:3)-prob_lo(1:3) ! Domain length
    
    !$omp parallel private(i,j,k,x,y,z)
    do k=lo(3),hi(3)
       z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)
             if ((x*x + y*y + z*z) .lt. 0.5d0) then
                rho(i,j,k,1)           = 0.5d0
                rho(i,j,k,2:nspecies)  = 0.3d0
                !diff_coeffs(i,j,k,:,:) = 1.d0
             else
                rho(i,j,k,1)           = 0.0d0
                rho(i,j,k,2:nspecies)  = 0.0d0
                !diff_coeffs(i,j,k,:,:) = 1.d0
             endif
          end do
       end do
    end do
    !$omp end parallel do
 
  end subroutine init_rho_3d

end module init_module
