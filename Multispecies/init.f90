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
  
  subroutine init_rho(rho,dx,prob_lo,the_bc_level)

    type(multifab) , intent(inout) :: rho(:)   !1st=nlevs 
    real(kind=dp_t), intent(in   ) :: dx(:,:)  !1st=nlevs, 2nd=dim         
    real(kind=dp_t), intent(in   ) :: prob_lo(rho(1)%dim)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local variables
    integer :: lo(rho(1)%dim), hi(rho(1)%dim)
    integer :: nspecies, dm, ng, i, n, nlevs

    real(kind=dp_t), pointer :: dp(:,:,:,:)
 
    ! Amit: how to populate nspecies using probin_multispecies_module, call?
    !nspecies= xxxxxx%nspecies 
    dm = rho(1)%dim
    ng = rho(1)%ng
    nlevs = size(rho,1)

    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp => dataptr(rho(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          
          select case(dm)
          case (2)
             call init_rho_2d(dp(:,:,1,1:nspecies), ng, lo, hi, prob_lo, dx(n,1)) 
          case (3)
             call init_rho_3d(dp(:,:,:,1:nspecies), ng, lo, hi, prob_lo, dx(n,1))
          end select
       end do

       ! filling up ghost cells for two adjacent grids at the same level
       ! includeing periodic domain boundary ghost cells
       call multifab_fill_boundary(rho(n))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(rho(n),1,1,nspecies,the_bc_level(n))
    end do

  end subroutine init_rho

  subroutine init_rho_2d(rho, ng, lo, hi, prob_lo, dx)

    integer          :: lo(2), hi(2), ng
    real(kind=dp_t) :: rho(lo(1)-ng:,lo(2)-ng:,:) ! Last dimension is 1:nspecies
    !real(kind=dp_t) :: diff_coeffs(lo(1)-ng:,lo(2)-ng:,:,:) !to hold D matrix
    real(kind=dp_t) :: prob_lo(2)
    real(kind=dp_t) :: dx
 
    ! local varables
    integer          :: i,j
    real(kind=dp_t) :: x,y,r2
 
    ! Amit: Caution: remember to change multifab_physbc accordingly
    do j=lo(2),hi(2)
       y = prob_lo(2) + (dble(j)+0.5d0) * dx
       do i=lo(1),hi(1)
          x = prob_lo(1) + (dble(i)+0.5d0) * dx
          if ((x*x + y*y) .lt. 0.5d0) then
             rho(i,j,2:nspecies)  = 0.3d0
             !diff_coeffs(i,j,:,:) = 1.d0
          else
             rho(i,j,1)           = 0.1d0
             !diff_coeffs(i,j,:,:) = 1.d0
          endif 
       end do
    end do

    end subroutine init_rho_2d

    subroutine init_rho_3d(rho, ng, lo, hi, prob_lo, dx)

    integer          :: lo(3), hi(3), ng
    real(kind=dp_t) :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    !real(kind=dp_t) :: diff_coeffs(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:,:)
    ! Example:
    ! real(dp_t) :: diff_coffs(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:,:) ! Last dimensions are (1:nspecies,1:nspecies)
    ! You can use D = diff_coeffs(i,j,k,1:nspecies,1:nspecies)
    real(kind=dp_t) :: prob_lo(3)
    real(kind=dp_t) :: dx
 
    ! local varables
    integer          :: i,j,k
    real(kind=dp_t) :: x,y,z,r2

    ! Amit: New omp parallel implementation as learned in summer school
    !$omp parallel default(none), reduction(+:rho) private(i,j,k,x,y,z,r2)
    !$omp do 
    do k=lo(3),hi(3)
       z = prob_lo(3) + (dble(k)+0.5d0) * dx
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0) * dx
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0) * dx
             if ((x*x + y*y + z*z) .lt. 0.5d0) then
                rho(i,j,k,2:nspecies)  = 0.3d0
                !diff_coeffs(i,j,k,:,:) = 1.d0
             else
                rho(i,j,k,1)           = 0.1d0
                !diff_coeffs(i,j,k,:,:) = 1.d0
             endif
          end do
       end do
    end do
    !$omp end do
    !$omp end parallel
 
  end subroutine init_rho_3d

end module init_module
