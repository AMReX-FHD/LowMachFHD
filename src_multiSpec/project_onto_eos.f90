module project_onto_eos_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use probin_common_module, only: rhobar
  use probin_gmres_module, only: mg_rel_tol
  use probin_multispecies_module, only: nspecies



  implicit none

  private

  public :: project_onto_eos

contains

  subroutine project_onto_eos(mla,rho)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rho(:)

    ! local
    integer i,n,nlevs,dm,ng_r
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: sp(:,:,:,:)

    nlevs = mla%nlevel
    dm    = mla%dim

    ng_r = rho(1)%ng

    do n=1,nlevs
       do i=1,nfabs(rho(n))
          sp  => dataptr(rho(n), i)
          lo = lwb(get_box(rho(n), i))
          hi = upb(get_box(rho(n), i))
          select case (dm)
          case (2)
             call project_onto_eos_2d(sp(:,:,1,:), ng_r, lo, hi)
          case (3)
             call project_onto_eos_3d(sp(:,:,:,:), ng_r, lo, hi)
          end select
       end do
    end do

  end subroutine project_onto_eos

  subroutine project_onto_eos_2d(rho,ng_r,lo,hi)

    ! L2 Projection onto EOS Constraint

    integer        , intent(in   ) :: lo(:),hi(:),ng_r
    real(kind=dp_t), intent(inout) :: rho(lo(1)-ng_r:,lo(2)-ng_r:,:)

    real(kind=dp_t) :: sum,sum2,sum_sq,A(nspecies),r_tmp(3),tmp1

    integer i,j,l,m,ncell

    ! number of cells on the grid
    ncell = (hi(1)-lo(1)+1)*(hi(2)-lo(2)+1)

    sum_sq = 0.d0
    do l=1,nspecies
       sum_sq = sum_sq + rhobar(l)**2
    end do

    do l=1,nspecies
       do m=1,nspecies
          A(m) = rhobar(l)*rhobar(m) / sum_sq
       end do

       sum = 0.d0

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ! sum = sum_ij (A(i)*rhol - sum_{m!=l} A(m)*rhom) - sum_ij (rhol)
             !     = sum_ij (A(1)-1.0)*rhol - sum_{m!=l} A(m)*rhom
             sum2 = 0.d0
             do m=1,nspecies
                if (m .ne. l) then
                   sum2 = sum2 + A(m)*rho(i,j,m)
                end if
             end do
             sum = sum + (A(l)-1.d0)*rho(i,j,l) - sum2

          end do
       end do

       sum = sum / dble(ncell)

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             if (l .ne. nspecies) then

                sum2 = 0.d0
                do m=1,nspecies
                   if (m .ne. l) then
                      sum2 = sum2 + A(m)*rho(i,j,m)
                   end if
                end do
                tmp1 = rho(i,j,l)
                rho(i,j,l) = A(l)*rho(i,j,l) - sum2 - sum

                if (abs(rho(i,j,l)-tmp1) .gt. 1.d-9) then
                   print*,'eos big adjustment',i,j,abs(rho(i,j,l)-tmp1)
                end if

             else

                sum2 = 0
                do m=1,nspecies-1
                   sum2 = sum2 + rho(i,j,m)/rhobar(m)
                end do
                tmp1 = rho(i,j,l)
                rho(i,j,l) = (1.d0 - sum2)*rhobar(l)

                if (abs(rho(i,j,l)-tmp1) .gt. 1.d-9) then
                   print*,'eos big adjustment',i,j,abs(rho(i,j,l)-tmp1)
                end if

             end if

          end do
       end do

    end do ! end loop over species

!    r_tmp(1:2)=0.0d0 ! Check the L2 norm of the change
!    do j=lo(2),hi(2)
!       do i=lo(1),hi(1)

          ! rho1 = A*rho1 - B*rho2 - sum
          !      = A*rho1 - B*(rho-rho1) - sum
!          r_tmp(3) = rho(i,j,2)
!          rho(i,j,2) = A*rho(i,j,2) - B*(rho(i,j,1)-rho(i,j,2)) - sum
!          r_tmp(1) = r_tmp(1) + (r_tmp(3)-rho(i,j,2))**2

          ! rho = rho1 + rho2 
          !     = rho1 + (1-rho1/rhobar1)*rhobar2
!          r_tmp(3) = rho(i,j,1)
!          rho(i,j,1) = rho(i,j,2) + (1.d0-rho(i,j,2)/rhobar(1))*rhobar(2)
!          r_tmp(2) = r_tmp(2) + (r_tmp(3)-rho(i,j,1))**2

!       end do
!    end do
    
!    r_tmp(1:2) = sqrt(r_tmp(1:2)/dble(ncell))
!    r_tmp(1:2) = r_tmp(1:2) / sqrt(rhobar(1)**2 + rhobar(2)**2) ! Relative change
!    if(any( r_tmp(1:2) > 1000*mg_rel_tol)) then      
!       call bl_warn('EOS adjustment exceeded Poisson solver tolerance')
!       print*,r_tmp(1:2),mg_rel_tol
!    end if

  end subroutine project_onto_eos_2d

  subroutine project_onto_eos_3d(rho,ng_r,lo,hi)

    ! L2 Projection onto EOS Constraint

    integer        , intent(in   ) :: lo(:),hi(:),ng_r
    real(kind=dp_t), intent(inout) :: rho(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:,:)

  end subroutine project_onto_eos_3d

end module project_onto_eos_module
