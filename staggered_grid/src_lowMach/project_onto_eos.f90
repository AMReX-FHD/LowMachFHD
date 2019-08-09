module project_onto_eos_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use probin_common_module, only: rhobar, nspecies, algorithm_type, rho0, total_volume
  use probin_gmres_module, only: gmres_rel_tol

  implicit none

  private

  public :: project_onto_eos

contains

  subroutine project_onto_eos(mla,rho)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rho(:)

    ! local
    integer :: i,l,n,nlevs,dm,ng_r,ng_t
    integer :: lo(mla%dim),hi(mla%dim)
    real(kind=dp_t) :: rhobar_sq

    real(kind=dp_t) :: sum_spec(nspecies), sum_spec_proc(nspecies)
    real(kind=dp_t) :: sum_change(nspecies), sum_change_proc(nspecies)
    
    type(multifab) :: rho_tilde(mla%nlevel)
    
    real(kind=dp_t), pointer :: dp1(:,:,:,:)
    real(kind=dp_t), pointer :: dp2(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"project_onto_eos")

    nlevs = mla%nlevel
    dm    = mla%dim

    do n=1,nlevs
       call multifab_build(rho_tilde(n),mla%la(n),nspecies,0)
    end do
    
    ng_r = rho(1)%ng
    ng_t = rho_tilde(1)%ng

    sum_spec(:) = 0.d0
    sum_spec_proc(:) = 0.d0
    sum_change(:) = 0.d0
    sum_change_proc(:) = 0.d0

    ! boussinesq algorithms set the final species to be rho0-rho(1:nspecies-1)
    if (algorithm_type .eq. 4 .or. algorithm_type .eq. 6) then
       
       do n=1,nlevs
          do i=1,nfabs(rho(n))
             dp1 => dataptr(rho(n), i)
             dp2 => dataptr(rho_tilde(n), i)
             lo = lwb(get_box(rho(n), i))
             hi = upb(get_box(rho(n), i))
             select case (dm)
             case (2)
                call project_onto_eos_bousq_2d(dp1(:,:,1,:), ng_r, lo, hi)
             case (3)
                call project_onto_eos_bousq_3d(dp1(:,:,:,:), ng_r, lo, hi)
             end select
          end do
       end do

    else

       ! compute sum_spec on each processor
       do n=1,nlevs
          do i=1,nfabs(rho(n))
             dp1 => dataptr(rho(n), i)
             dp2 => dataptr(rho_tilde(n), i)
             lo = lwb(get_box(rho(n), i))
             hi = upb(get_box(rho(n), i))
             select case (dm)
             case (2)
                call project_onto_eos1_2d(dp1(:,:,1,:), ng_r, dp2(:,:,1,:), ng_t, &
                                          sum_spec_proc, lo, hi)
             case (3)
                call project_onto_eos1_3d(dp1(:,:,:,:), ng_r, dp2(:,:,:,:), ng_t, &
                                          sum_spec_proc, lo, hi)
             end select
          end do
       end do


       ! collect sum_spec(:)
       call parallel_reduce(sum_spec(1:nspecies), sum_spec_proc(1:nspecies), MPI_SUM)

       ! divide sum_spec(:) by total volume
       sum_spec = sum_spec / dble(total_volume)

       ! update rho and compute sum_change on each processor
       do n=1,nlevs
          do i=1,nfabs(rho(n))
             dp1 => dataptr(rho(n), i)
             dp2 => dataptr(rho_tilde(n), i)
             lo = lwb(get_box(rho(n), i))
             hi = upb(get_box(rho(n), i))
             select case (dm)
             case (2)
                call project_onto_eos2_2d(dp1(:,:,1,:), ng_r, dp2(:,:,1,:), ng_t, &
                                          sum_spec, sum_change_proc, lo, hi)
             case (3)
                call project_onto_eos2_3d(dp1(:,:,:,:), ng_r, dp2(:,:,:,:), ng_t, &
                                          sum_spec, sum_change_proc, lo, hi)
             end select
          end do
       end do

       ! collect sum_change
       call parallel_reduce(sum_change(1:nspecies), sum_change_proc(1:nspecies), MPI_SUM)
       
       ! redefine rhobar_sq = sum(rhobar_i^2)
       rhobar_sq = 0.d0
       do l=1,nspecies
          rhobar_sq = rhobar_sq + rhobar(l)**2
       end do
       sum_change(:) = sqrt(sum_change(:)/dble(total_volume))
       sum_change(:) = sum_change(:) / sqrt(rhobar_sq)
       if(any( sum_change(1:nspecies) > 1000*gmres_rel_tol)) then      
          call bl_warn('EOS adjustment exceeded GMRES solver tolerance')
          print*,sum_change(1:nspecies),gmres_rel_tol
       end if
       
    end if

    do n=1,nlevs
       call multifab_destroy(rho_tilde(n))
    end do
       
    call destroy(bpt)

  end subroutine project_onto_eos

  subroutine project_onto_eos_bousq_2d(rho,ng_r,lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:),ng_r
    real(kind=dp_t), intent(inout) :: rho(lo(1)-ng_r:,lo(2)-ng_r:,:)

    integer i,j,l,ncell

    ! boussinesq
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)
       rho(i,j,nspecies) = rho0 - sum(rho(i,j,1:nspecies-1))
    end do
    end do

  end subroutine project_onto_eos_bousq_2d

  subroutine project_onto_eos_bousq_3d(rho,ng_r,lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:),ng_r
    real(kind=dp_t), intent(inout) :: rho(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:,:)
    
    integer i,j,k

    ! boussinesq
    do k=lo(3),hi(3)
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)
       rho(i,j,k,nspecies) = rho0 - sum(rho(i,j,k,1:nspecies-1))
    end do
    end do
    end do

  end subroutine project_onto_eos_bousq_3d

  subroutine project_onto_eos1_2d(rho,ng_r,rho_tilde,ng_t,sum_spec,lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:),ng_r,ng_t
    real(kind=dp_t), intent(in   ) :: rho      (lo(1)-ng_r:,lo(2)-ng_r:,:)
    real(kind=dp_t), intent(inout) :: rho_tilde(lo(1)-ng_t:,lo(2)-ng_t:,:)
    real(kind=dp_t), intent(inout) :: sum_spec(nspecies)

    real(kind=dp_t) :: rhobar_sq,delta_eos
    real(kind=dp_t) :: rho_tmp,w(nspecies)
    integer i,j,l,ncell

    ! L2 Projection onto EOS Constraint
    ! number of cells on the grid
    ncell = (hi(1)-lo(1)+1)*(hi(2)-lo(2)+1)

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          ! compute mass fractions, w_i = rho_i/rho
          rho_tmp = 0.d0
          do l=1,nspecies
             rho_tmp = rho_tmp + rho(i,j,l)
          end do
          do l=1,nspecies
             w(l) = rho(i,j,l)/rho_tmp
          end do

          ! rhobar_sq = (sum_i (w_i/rhobar_i^2))^-1
          rhobar_sq = 0.d0
          do l=1,nspecies
             rhobar_sq = rhobar_sq + w(l)/rhobar(l)**2
          end do
          rhobar_sq = 1.d0/rhobar_sq

          ! delta_eos = sum_l (rho_l/rhobar_l) - 1
          delta_eos = -1.d0
          do l=1,nspecies
             delta_eos = delta_eos + rho(i,j,l)/rhobar(l)
          end do

          do l=1,nspecies
             ! rho_tilde_i = rho - w_i*(rhobar_sq/rhobar_i) * delta_eos
             rho_tilde(i,j,l) = rho(i,j,l) - w(l)*(rhobar_sq/rhobar(l))*delta_eos
             ! sum_spec_i = sum (rho_i - rho_tilde_i)
             sum_spec(l) = sum_spec(l) + rho(i,j,l) - rho_tilde(i,j,l)
          end do

       end do
    end do

  end subroutine project_onto_eos1_2d

  subroutine project_onto_eos1_3d(rho,ng_r,rho_tilde,ng_t,sum_spec,lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:),ng_r,ng_t
    real(kind=dp_t), intent(in   ) :: rho      (lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:,:)
    real(kind=dp_t), intent(inout) :: rho_tilde(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:,:)
    real(kind=dp_t), intent(inout) :: sum_spec(nspecies)

    real(kind=dp_t) :: rhobar_sq,delta_eos
    real(kind=dp_t) :: rho_tmp,w(nspecies)
    integer i,j,k,l,ncell

    ! L2 Projection onto EOS Constraint
    ! number of cells on the grid
    ncell = (hi(1)-lo(1)+1)*(hi(2)-lo(2)+1)*(hi(3)-lo(3)+1)

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ! compute mass fractions, w_i = rho_i/rho
             rho_tmp = 0.d0
             do l=1,nspecies
                rho_tmp = rho_tmp + rho(i,j,k,l)
             end do
             do l=1,nspecies
                w(l) = rho(i,j,k,l)/rho_tmp
             end do

             ! rhobar_sq = (sum_i (w_i/rhobar_i^2))^-1
             rhobar_sq = 0.d0
             do l=1,nspecies
                rhobar_sq = rhobar_sq + w(l)/rhobar(l)**2
             end do
             rhobar_sq = 1.d0/rhobar_sq

             ! delta_eos = sum_l (rho_l/rhobar_l) - 1
             delta_eos = -1.d0
             do l=1,nspecies
                delta_eos = delta_eos + rho(i,j,k,l)/rhobar(l)
             end do

             do l=1,nspecies
                ! rho_tilde_i = rho - w_i*(rhobar_sq/rhobar_i) * delta_eos
                rho_tilde(i,j,k,l) = rho(i,j,k,l) - w(l)*(rhobar_sq/rhobar(l))*delta_eos
                ! sum_spec_i = sum (rho_i - rho_tilde_i)
                sum_spec(l) = sum_spec(l) + rho(i,j,k,l) - rho_tilde(i,j,k,l)
             end do

          end do
       end do
    end do

  end subroutine project_onto_eos1_3d

    subroutine project_onto_eos2_2d(rho,ng_r,rho_tilde,ng_t,sum_spec,sum_change,lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:),ng_r,ng_t
    real(kind=dp_t), intent(inout) :: rho      (lo(1)-ng_r:,lo(2)-ng_r:,:)
    real(kind=dp_t), intent(inout) :: rho_tilde(lo(1)-ng_t:,lo(2)-ng_t:,:)
    real(kind=dp_t), intent(in   ) :: sum_spec(nspecies)
    real(kind=dp_t), intent(inout) :: sum_change(nspecies)

    real(kind=dp_t) :: rho_tmp
    integer i,j,l

    do j=lo(2),hi(2)
    do i=lo(1),hi(1)
       do l=1,nspecies
          rho_tmp = rho(i,j,l)
          rho(i,j,l) = rho_tilde(i,j,l) + sum_spec(l)
          sum_change(l) = sum_change(l) + (rho(i,j,l)-rho_tmp)**2
       end do
    end do
    end do

  end subroutine project_onto_eos2_2d

  subroutine project_onto_eos2_3d(rho,ng_r,rho_tilde,ng_t,sum_spec,sum_change,lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:),ng_r,ng_t
    real(kind=dp_t), intent(inout) :: rho      (lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:,:)
    real(kind=dp_t), intent(inout) :: rho_tilde(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:,:)
    real(kind=dp_t), intent(in   ) :: sum_spec(nspecies)
    real(kind=dp_t), intent(inout) :: sum_change(nspecies)

    real(kind=dp_t) :: rho_tmp
    integer i,j,k,l

    do k=lo(3),hi(3)
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)
       do l=1,nspecies
          rho_tmp = rho(i,j,k,l)
          rho(i,j,k,l) = rho_tilde(i,j,k,l) + sum_spec(l)
          sum_change(l) = sum_change(l) + (rho(i,j,k,l)-rho_tmp)**2
       end do
    end do
    end do
    end do

  end subroutine project_onto_eos2_3d

end module project_onto_eos_module
