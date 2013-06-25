module init_module

  use multifab_module
  use define_bc_module
  use multifab_physbc_module
  ! Donev: Variables you need here:
  use probin_common_module
  !use probin_multispecies_module
  
  ! Donev:
  ! replace "double precision" with "real(dp_t)" everywhere in all codes

  implicit none

  private

  public :: init_rho

contains
  
  subroutine init_rho(rho,dx,prob_lo,the_bc_tower)

    type(multifab) , intent(inout) :: rho(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(in   ) :: prob_lo(rho(1)%dim)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: lo(rho(1)%dim), hi(rho(1)%dim)
    integer :: nsp, dm, ng, i, n, nlevs

    real(kind=dp_t), pointer :: dp(:,:,:,:)
 
    nsp= 4 ! Amit: have to figureout way to pass this from main.f90.
    dm = rho(1)%dim
    ng = rho(1)%ng

    nlevs = size(rho,1)

    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp => dataptr(rho(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          ! Donev: This is wrong, it only handles the last species
          ! You need 1:nsp, not nsp as the last index.
          ! Fix this
          ! Amit : Have done this earlier but somehow the 4-dimensional rho
          ! is not getting passed correctly neither with : or 1:nsp.   
 
          select case(dm)
          case (2)
             call init_rho_2d(dp(:,:,1,1:nsp), ng, lo, hi, prob_lo, dx(n,1)) 
          case (3)
             call init_rho_3d(dp(:,:,:,1:nsp), ng, lo, hi, prob_lo, dx(n,1))
          end select
       end do

       ! filling up ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(rho(n))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(rho(n),1,1,nsp,the_bc_tower%bc_tower_array(n))
    end do

  end subroutine init_rho

  subroutine init_rho_2d(rho, ng, lo, hi, prob_lo, dx)

    integer          :: lo(2), hi(2), ng
    ! Donev: These are higher-dimensional arrays
    double precision :: rho(lo(1)-ng:,lo(2)-ng:,:) ! Last dimension is 1:nspecies
    double precision :: prob_lo(2)
    double precision :: dx
 
    ! local varables
    integer          :: i,j
    double precision :: x,y,r2

    ! Donev: You need to set values for rho(i,j,1:nspecies)
    do j=lo(2),hi(2)
       y = prob_lo(2) + (dble(j)+0.5d0) * dx
       do i=lo(1),hi(1)
          x = prob_lo(1) + (dble(i)+0.5d0) * dx
          if ((x*x + y*y) .lt. 0.5d0) then
             rho(i,j) = 0.4d0
             if((x*x + y*y) .lt. 0.45d0) then
             rho(i,j) = 0.3d0
             else if((x*x + y*y) .lt. 0.35d0) then
             rho(i,j) = 0.25d0
             endif
          else
             rho(i,j) = 0.d0
          endif 
       end do
    end do
    ! Amit: Currently species are initialized in concetric circles with density .4, .3, 
    ! .25 and 0. The outerone is 0 now because have to manually put that number in 
    ! multifab_physbc for ghost cells and boundary condition. I'll copy that
    ! file in this folder because we have to keep track of the numbers there.

    end subroutine init_rho_2d

    subroutine init_rho_3d(rho, ng, lo, hi, prob_lo, dx)

    integer          :: lo(3), hi(3), ng
    ! Donev: Make this a rank-4 array
    double precision :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    ! Example:
    ! real(dp_t) :: diff_coffs(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:,:) ! Last dimensions are (1:nspecies,1:nspecies)
    ! You can use D = diff_coeffs(i,j,k,1:nspecies,1:nspecies)
    double precision :: prob_lo(3)
    double precision :: dx
 
    ! local varables
    integer          :: i,j,k
    double precision :: x,y,z,r2

    !$omp parallel default(none), reduction(+:rho) private(i,j,k,x,y,z,r2)
    !$omp do 
    do k=lo(3),hi(3)
       z = prob_lo(3) + (dble(k)+0.5d0) * dx
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0) * dx
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0) * dx
             if ((x*x + y*y + z*z) .lt. 0.5d0) then
                rho(i,j,k) = 0.4d0
                if((x*x + y*y + z*z) .lt. 0.45d0) then
                rho(i,j,k) = 0.3d0
                else if((x*x + y*y + z*z) .lt. 0.35d0) then
                rho(i,j,k) = 0.25d0
                endif
             else
                rho(i,j,k) = 0.d0
             endif
          end do
       end do
    end do
    !$omp end do
    !$omp end parallel
    !Amit: New omp parallel implementation as learned in summer school
 
  end subroutine init_rho_3d

end module init_module
