module init_rho_module

  use multifab_module
  use define_bc_module
  use multifab_physbc_module

  implicit none

  private

  public :: init_rho

contains
  
  subroutine init_rho(rho,dx,prob_lo,the_bc_tower)

    type(multifab) , intent(inout) :: rho
    real(kind=dp_t), intent(in   ) :: dx
    real(kind=dp_t), intent(in   ) :: prob_lo(rho%dim)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: lo(rho%dim), hi(rho%dim)
    integer :: nsp, dm, ng, i

    real(kind=dp_t), pointer :: dp(:,:,:,:)
 
    nsp= 4 ! Amit: what is the way to pass from main.f90 ?
    dm = rho%dim
    ng = rho%ng

    do i=1,nfabs(rho)
       dp => dataptr(rho,i)
       lo = lwb(get_box(rho,i))
       hi = upb(get_box(rho,i))
       ! Donev: This is wrong, it only handles the last species
       ! You need 1:nsp, not nsp as the last index.
       ! Fix this
       select case(dm)
       case (2)
          call init_rho_2d(dp(:,:,1,nsp), ng, lo, hi, prob_lo, dx) 
       case (3)
          call init_rho_3d(dp(:,:,:,nsp), ng, lo, hi, prob_lo, dx)
    ! Amit: 4th index is replaced with nsp instead of 1 to include 
    ! all species for both case (2) and (3).
       end select
    end do

    ! filling up ghost cells for two adjacent grids at the same level
    ! this includes periodic domain boundary ghost cells
    call multifab_fill_boundary(rho)

    ! fill non-periodic domain boundary ghost cells
    ! Amit: Confusing 1,1, entry. 
    call multifab_physbc(rho,1,1,nsp,the_bc_tower%bc_tower_array(1))

  end subroutine init_rho

  subroutine init_rho_2d(rho, ng, lo, hi, prob_lo, dx)

    integer          :: lo(2), hi(2), ng
    double precision :: rho(lo(1)-ng:,lo(2)-ng:)
    double precision :: prob_lo(2)
    double precision :: dx
 
    ! local varables
    integer          :: i,j
    double precision :: x,y,r2

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
    ! multifab_physbc for ghost cells and boundary condition.
    end subroutine init_rho_2d

    subroutine init_rho_3d(rho, ng, lo, hi, prob_lo, dx)

    integer          :: lo(3), hi(3), ng
    double precision :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    double precision :: prob_lo(3)
    double precision :: dx
 
    ! local varables
    integer          :: i,j,k
    double precision :: x,y,z,r2

    !$omp parallel do private(i,j,k,x,y,z,r2)
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
    !$omp end parallel do

  end subroutine init_rho_3d

end module init_rho_module
