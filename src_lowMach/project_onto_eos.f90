module project_onto_eos_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use probin_lowmach_module, only: rhobar
  use probin_gmres_module  , only: mg_rel_tol



  implicit none

  private

  public :: project_onto_eos

contains

  subroutine project_onto_eos(mla,snew)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: snew(:)

    ! local
    integer i,n,nlevs,dm,ng_s
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: sp(:,:,:,:)

    nlevs = mla%nlevel
    dm    = mla%dim

    ng_s = snew(1)%ng

    do n=1,nlevs
       do i=1,nfabs(snew(n))
          sp  => dataptr(snew(n), i)
          lo = lwb(get_box(snew(n), i))
          hi = upb(get_box(snew(n), i))
          select case (dm)
          case (2)
             call project_onto_eos_2d(sp(:,:,1,:), ng_s, lo, hi)
          case (3)
             call project_onto_eos_3d(sp(:,:,:,:), ng_s, lo, hi)
          end select
       end do
    end do

  end subroutine project_onto_eos

  subroutine project_onto_eos_2d(snew,ng_s,lo,hi)

    ! L2 Projection onto EOS Constraint

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    real(kind=dp_t), intent(inout) :: snew(lo(1)-ng_s:,lo(2)-ng_s:,:)

    real(kind(1.0d0)) :: A,B,sum, r_tmp(3)

    integer i,j,ncell

    ! NOTE: comp=1 corresponds to rho
    !       comp=2 corresponds to rho1 = rho*c
    !       We do not store rho2 = rho*(1-c)

    A = rhobar(1)**2 / (rhobar(1)**2 + rhobar(2)**2)
    B = rhobar(1)*rhobar(2) / (rhobar(1)**2 + rhobar(2)**2)
    ncell = (hi(1)-lo(1)+1)*(hi(2)-lo(2)+1)
    sum = 0.d0

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          ! sum = sum_ij (A*rho1 - B*rho2) - sum_ij (rho1)
          !     = sum_ij (A*rho1 - B*(rho-rho1)) - sum_ij (rho1)
          sum = sum + (A-1.d0)*snew(i,j,2) - B*(snew(i,j,1)-snew(i,j,2))

       end do
    end do

    sum = sum / dble(ncell)

    r_tmp(1:2)=0.0d0 ! Check the L2 norm of the change
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          ! rho1 = A*rho1 - B*rho2 - sum
          !      = A*rho1 - B*(rho-rho1) - sum
          r_tmp(3) = snew(i,j,2)
          snew(i,j,2) = A*snew(i,j,2) - B*(snew(i,j,1)-snew(i,j,2)) - sum
          r_tmp(1) = r_tmp(1) + (r_tmp(3)-snew(i,j,2))**2

          ! rho = rho1 + rho2 
          !     = rho1 + (1-rho1/rhobar1)*rhobar2
          r_tmp(3) = snew(i,j,1)
          snew(i,j,1) = snew(i,j,2) + (1.d0-snew(i,j,2)/rhobar(1))*rhobar(2)
          r_tmp(2) = r_tmp(2) + (r_tmp(3)-snew(i,j,1))**2

       end do
    end do
    
    r_tmp(1:2) = sqrt(r_tmp(1:2)/dble(ncell))
    r_tmp(1:2) = r_tmp(1:2) / sqrt(rhobar(1)**2 + rhobar(2)**2) ! Relative change
    if(any( r_tmp(1:2) > 1000*mg_rel_tol)) then      
       call bl_warn('EOS adjustment exceeded Poisson solver tolerance')
       print*,r_tmp(1:2),mg_rel_tol
    end if

  end subroutine project_onto_eos_2d

  subroutine project_onto_eos_3d(snew,ng_s,lo,hi)

    ! L2 Projection onto EOS Constraint

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    real(kind=dp_t), intent(inout) :: snew(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)

    real(kind(1.0d0)) :: A,B,sum, r_tmp(3)

    integer i,j,k,ncell

    ! NOTE: comp=1 corresponds to rho
    !       comp=2 corresponds to rho1 = rho*c
    !       We do not store rho2 = rho*(1-c)

    A = rhobar(1)**2 / (rhobar(1)**2 + rhobar(2)**2)
    B = rhobar(1)*rhobar(2) / (rhobar(1)**2 + rhobar(2)**2)
    ncell = (hi(1)-lo(1)+1)*(hi(2)-lo(2)+1)*(hi(3)-lo(3)+1)
    sum = 0.d0
    
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ! sum = sum_ij (A*rho1 - B*rho2)
             !     = sum_ij (A*rho1 - B*(rho-rho1))
             sum = sum + (A-1.d0)*snew(i,j,k,2) - B*(snew(i,j,k,1)-snew(i,j,k,2))

          end do
       end do
    end do

    sum = sum / dble(ncell)

    r_tmp(1:2)=0.0d0 ! Check the L2 norm of the change
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ! rho1 = A*rho1 - B*rho2 - sum
             !      = A*rho1 - B*(rho-rho1) - sum
             r_tmp(3) = snew(i,j,k,2)
             snew(i,j,k,2) = A*snew(i,j,k,2) &
                  - B*(snew(i,j,k,1)-snew(i,j,k,2)) - sum
             r_tmp(1) = r_tmp(1) + (r_tmp(3)-snew(i,j,k,2))**2

             ! rho = rho1 + rho2 
             !     = rho1 + (1-rho1/rhobar1)*rhobar2
             r_tmp(3) = snew(i,j,k,1)
             snew(i,j,k,1) = snew(i,j,k,2) + (1.d0-snew(i,j,k,2)/rhobar(1))*rhobar(2)
             r_tmp(2) = r_tmp(2) + (r_tmp(3)-snew(i,j,k,1))**2

          end do
       end do
    end do

    r_tmp(1:2) = sqrt(r_tmp(1:2)/dble(ncell))
    r_tmp(1:2) = r_tmp(1:2) / sqrt(rhobar(1)**2 + rhobar(2)**2) ! Relative change
    if(any( r_tmp(1:2) > 1000*mg_rel_tol)) then      
       call bl_warn('EOS adjustment exceeded Poisson solver tolerance')
       print*,r_tmp(1:2),mg_rel_tol
    end if

  end subroutine project_onto_eos_3d

end module project_onto_eos_module
