module eos_check_module

  use multifab_module
  use ml_layout_module
  use convert_rhoc_to_c_module
  use probin_common_module, only: rhobar, nspecies, algorithm_type, rho0, rho_eos_form

  implicit none

  private

  public :: eos_check, compute_rhotot_eos, compute_rho_eos

contains

  subroutine eos_check(mla,rho)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: rho(:)

    ! local
    integer i,n,dm,nlevs,ng_r
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: rp(:,:,:,:)

    real(kind=dp_t) :: eos_error, eos_error_grid, eos_error_proc

    type(bl_prof_timer), save :: bpt

    call build(bpt, "eos_check")

    nlevs = mla%nlevel
    dm = mla%dim

    ng_r = rho(1)%ng

    eos_error_proc = -1.d20
    do n=1,nlevs
       do i=1,nfabs(rho(n))
         rp => dataptr(rho(n), i)
         lo =  lwb(get_box(rho(n), i))
         hi =  upb(get_box(rho(n), i))
         eos_error_grid = -1.d20
         select case (dm)
         case (2)
            call eos_check_2d(rp(:,:,1,:),ng_r,eos_error_grid,lo,hi)
         case (3)
            call eos_check_3d(rp(:,:,:,:),ng_r,eos_error_grid,lo,hi)
         end select
         eos_error_proc = max(eos_error_grid, eos_error_proc)
      end do
   end do

   ! This sets eos_error to be the max of eos_error_proc over all processors.
   call parallel_reduce(eos_error, eos_error_proc, MPI_MAX)

   if (parallel_IOProcessor()) then
      print*,"EOS ERROR in L1 norm: ",eos_error
      print*,""
   end if

   call destroy(bpt)

  end subroutine eos_check

  subroutine eos_check_2d(rho,ng_r,eos_error,lo,hi)

    integer        , intent(in   ) :: lo(:), hi(:), ng_r
    real(kind=dp_t), intent(in   ) :: rho(lo(1)-ng_r:,lo(2)-ng_r:,:)
    real(kind=dp_t), intent(inout) :: eos_error

    ! local
    integer :: i,j,n

    real(kind=dp_t) :: sum,error

    if (algorithm_type .eq. 6) then

       do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          sum = 0.d0
          do n=1,nspecies
             sum = sum + rho(i,j,n)
          end do
          error = abs(rho0 - sum)
          eos_error = max(eos_error,error)

       end do
       end do

    else

       do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          sum = 0.d0
          do n=1,nspecies
             sum = sum + rho(i,j,n)/rhobar(n)
          end do
          error = abs(1.d0 - sum)
          eos_error = max(eos_error,error)

       end do
       end do

    end if

  end subroutine eos_check_2d

  subroutine eos_check_3d(rho,ng_r,eos_error,lo,hi)

    integer        , intent(in   ) :: lo(:), hi(:), ng_r
    real(kind=dp_t), intent(in   ) :: rho(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:,:)
    real(kind=dp_t), intent(inout) :: eos_error

    ! local
    integer :: i,j,k,n

    real(kind=dp_t) :: sum,error

    if (algorithm_type .eq. 6) then

       do k=lo(3),hi(3)
       do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          sum = 0.d0
          do n=1,nspecies
             sum = sum + rho(i,j,k,n)
          end do
          error = abs(rho0 - sum)
          eos_error = max(eos_error,error)
             
       end do
       end do
       end do

    else

       do k=lo(3),hi(3)
       do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          sum = 0.d0
          do n=1,nspecies
             sum = sum + rho(i,j,k,n)/rhobar(n)
          end do
          error = abs(1.d0 - sum)
          eos_error = max(eos_error,error)
             
       end do
       end do
       end do

    end if

  end subroutine eos_check_3d
  

  subroutine compute_rho_eos(rho, rhotot)
    real(kind=dp_t), intent(in ) :: rho(nspecies)
    real(kind=dp_t), intent(out) :: rhotot
    
    real(kind=dp_t) :: c(nspecies)

    ! compute rhotot via sum
    rhotot = sum(rho)
    ! compute concentrations
    c=rho/rhotot  
   
    if(rho_eos_form==2) then ! Linearized EOS of incompressible components    
       rhotot = rho0  - rho0* sum((rho0/rhobar(1:nspecies-1)-1.0d0)*c(1:nspecies-1))
    else ! Nonlinear EOS of incompressible components    
       rhotot = 1.0d0 / sum(c/rhobar(1:nspecies))       
    end if   
  
  end subroutine compute_rho_eos


  subroutine compute_rhotot_eos(mla,rho,rhotot_eos,rhotot_eos_comp)

    ! compute rhotot_eos = sum(c_i/rhobar_i) in the valid region

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rho(:)
    type(multifab) , intent(inout) :: rhotot_eos(:)
    integer        , intent(in   ) :: rhotot_eos_comp

    ! local
    integer i,n,dm,nlevs,ng_r,ng_rr
    real(kind=dp_t), pointer :: rp(:,:,:,:), rpp(:,:,:,:)
    integer :: lo(mla%dim),hi(mla%dim)

    type(bl_prof_timer),save :: bpt

    call build(bpt,"compute_rhotot_eos")

    nlevs = mla%nlevel
    dm = mla%dim
    ng_r = rho(1)%ng
    ng_rr = rhotot_eos(1)%ng

    do n=1,nlevs
       do i=1,nfabs(rho(n))
         rp => dataptr(rho(n), i)
         rpp => dataptr(rhotot_eos(n), i)
         lo =  lwb(get_box(rho(n), i))
         hi =  upb(get_box(rho(n), i))
         select case (dm)
         case (2)
            call compute_rhotot_eos_2d(rp(:,:,1,:),ng_r,rpp(:,:,1,rhotot_eos_comp),ng_rr,lo,hi)
         case (3)
            call compute_rhotot_eos_3d(rp(:,:,:,:),ng_r,rpp(:,:,:,rhotot_eos_comp),ng_rr,lo,hi)
         end select
         
      end do
    end do

    call destroy(bpt)

  end subroutine compute_rhotot_eos

  subroutine compute_rhotot_eos_2d(rho,ng_r,rhotot_eos,ng_rr,lo,hi)

    integer        , intent(in   ) :: lo(:), hi(:), ng_r, ng_rr
    real(kind=dp_t), intent(in   ) :: rho(lo(1)-ng_r:,lo(2)-ng_r:,:)
    real(kind=dp_t), intent(inout) :: rhotot_eos(lo(1)-ng_rr:,lo(2)-ng_rr:)

    ! local
    integer :: i,j
    
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)    
       call compute_rho_eos(rho(i,j,:), rhotot_eos(i,j))
    end do
    end do

  end subroutine compute_rhotot_eos_2d

  subroutine compute_rhotot_eos_3d(rho,ng_r,rhotot_eos,ng_rr,lo,hi)

    integer        , intent(in   ) :: lo(:), hi(:), ng_r, ng_rr
    real(kind=dp_t), intent(in   ) :: rho(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:,:)
    real(kind=dp_t), intent(inout) :: rhotot_eos(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:)

    ! local
    integer :: i,j,k

    do k=lo(3),hi(3)
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)    
       call compute_rho_eos(rho(i,j,k,:), rhotot_eos(i,j,k))
    end do
    end do
    end do

  end subroutine compute_rhotot_eos_3d
  
end module eos_check_module
