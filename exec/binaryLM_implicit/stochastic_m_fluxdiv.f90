module stochastic_m_fluxdiv_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use BoxLibRNGs
  use bc_module
  use multifab_physbc_stag_module
  use multifab_fill_random_module
  use multifab_filter_module
  use sum_momenta_module
  use probin_common_module , only: visc_type, visc_coef, variance_coef, k_B, &
                                   stoch_stress_form, filtering_width
  use probin_binarylm_module, only: temperature

  implicit none

  private

  public :: stochastic_m_fluxdiv, fill_m_stochastic, &
       init_m_stochastic, destroy_m_stochastic, add_m_fluctuations

  ! Stochastic fluxes for momentum are generated on:
  ! -cell-centered grid for diagonal components
  ! -node-centered (2D) or edge-centered (3D) grid for off-diagonal components
  type(multifab), allocatable, save :: mflux_cc(:,:), mflux_nd(:,:), mflux_ed(:,:,:)
  
  integer, save :: n_rngs ! how many random number stages
  
contains

  ! Note that here we *increment* stoch_m_force so it must be initialized externally!
  subroutine stochastic_m_fluxdiv(mla,the_bc_level,stoch_m_force,eta,eta_ed,dx,dt,weights)
    
    type(ml_layout), intent(in   ) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(multifab) , intent(inout) :: stoch_m_force(:,:)
    type(multifab) , intent(in   ) :: eta(:)
    type(multifab) , intent(in   ) :: eta_ed(:,:)
    real(dp_t)     , intent(in   ) :: dx(:,:),dt,weights(:)

    ! local
    integer :: n,nlevs,dm,i,comp
    integer :: ng_c,ng_e,ng_y,ng_w,ng_n,ng_f

    real(dp_t) :: variance

    type(multifab) :: mflux_cc_temp(mla%nlevel)
    type(multifab) :: mflux_nd_temp(mla%nlevel)
    type(multifab) :: mflux_ed_temp(mla%nlevel,3)

    logical :: nodal_temp(mla%dim)

    real(kind=dp_t), pointer :: fp(:,:,:,:), dp(:,:,:,:), sp(:,:,:,:)
    real(kind=dp_t), pointer :: fxp(:,:,:,:), fyp(:,:,:,:), fzp(:,:,:,:)
    real(kind=dp_t), pointer :: dxp(:,:,:,:), dyp(:,:,:,:), dzp(:,:,:,:)
    real(kind=dp_t), pointer :: ep1(:,:,:,:), ep2(:,:,:,:), ep3(:,:,:,:)
    integer :: lo(mla%dim), hi(mla%dim)

    nlevs = mla%nlevel
    dm = mla%dim

    do n=1,nlevs

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! create a place to store temporary copy of the random numbers
       ! these copies will be scaled and used to create the stochastic flux divergence
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! we need dm cell-centered fluxes for momentum
       call multifab_build(mflux_cc_temp(n),mla%la(n),dm,max(1,filtering_width))
       if (dm .eq. 2) then
          ! in 2D, we need 2 random fluxes at each node
          nodal_temp = .true.
          call multifab_build(mflux_nd_temp(n),mla%la(n),2,filtering_width,nodal_temp)
       else if (dm .eq. 3) then
          ! in 3D, we need 2 random fluxes at each edge
          nodal_temp(1) = .true.
          nodal_temp(2) = .true.
          nodal_temp(3) = .false.
          call multifab_build(mflux_ed_temp(n,1),mla%la(n),2,filtering_width,nodal_temp)
          nodal_temp(1) = .true.
          nodal_temp(2) = .false.
          nodal_temp(3) = .true.
          call multifab_build(mflux_ed_temp(n,2),mla%la(n),2,filtering_width,nodal_temp)
          nodal_temp(1) = .false.
          nodal_temp(2) = .true.
          nodal_temp(3) = .true.
          call multifab_build(mflux_ed_temp(n,3),mla%la(n),2,filtering_width,nodal_temp)
       end if

       call multifab_setval(mflux_cc_temp(n),0.d0,all=.true.)
       if (dm .eq. 2) then
          call multifab_setval(mflux_nd_temp(n),0.d0,all=.true.)
       else if (dm .eq. 3) then
          call multifab_setval(mflux_ed_temp(n,1),0.d0,all=.true.)
          call multifab_setval(mflux_ed_temp(n,2),0.d0,all=.true.)
          call multifab_setval(mflux_ed_temp(n,3),0.d0,all=.true.)
       end if

       ! add weighted contribution of fluxes
       do comp=1,n_rngs
          call saxpy(mflux_cc_temp(n),weights(comp),mflux_cc(n,comp),all=.true.)
          if (dm .eq. 2) then
             call saxpy(mflux_nd_temp(n),weights(comp),mflux_nd(n,comp),all=.true.)
          else if (dm .eq. 3) then
             call saxpy(mflux_ed_temp(n,1),weights(comp),mflux_ed(n,1,comp),all=.true.)
             call saxpy(mflux_ed_temp(n,2),weights(comp),mflux_ed(n,2,comp),all=.true.)
             call saxpy(mflux_ed_temp(n,3),weights(comp),mflux_ed(n,3,comp),all=.true.)
          end if
       end do

    end do

    ng_c = mflux_cc_temp(1)%ng
    ng_y = eta(1)%ng
    ng_f = stoch_m_force(1,1)%ng

    do n=1,nlevs

       if (visc_type < 0) then
          ! eta varies in space, add its contribution below in an i/j/k loop
          variance = sqrt(variance_coef*2.d0*k_B*temperature          /(product(dx(n,1:dm))*dt))
       else
          ! eta is constant in space, include it here
          variance = sqrt(variance_coef*2.d0*k_B*temperature*visc_coef/(product(dx(n,1:dm))*dt))
       end if


       if (dm .eq. 2) then

          ng_n = mflux_nd_temp(1)%ng
          ng_w = eta_ed(1,1)%ng
          
          do i=1,nfabs(mflux_cc_temp(n))
             
             fp  => dataptr(mflux_cc_temp(n),i)
             sp  => dataptr(mflux_nd_temp(n),i)
             dp  => dataptr(eta(n),i)
             lo = lwb(get_box(mflux_cc_temp(n),i))
             hi = upb(get_box(mflux_cc_temp(n),i))
             ! multiply by variance
             fp = variance*fp
             sp = variance*sp
             ! if eta varies in space, multiply pointwise by sqrt(eta)
             if (visc_type < 0) then
                ep1 => dataptr(eta_ed(n,1),i)
                call mult_by_sqrt_eta_2d(fp(:,:,1,:),ng_c,sp(:,:,1,:),ng_n, &
                                         dp(:,:,1,1),ng_y,ep1(:,:,1,1),ng_w,lo,hi)
             end if
             ! apply boundary conditions
             call mflux_bc_2d(sp(:,:,1,:),ng_n,lo,hi, &
                              the_bc_level(n)%phys_bc_level_array(i,:,:))
          end do

          ! sync up random numbers at boundaries and ghost cells
          call multifab_internal_sync(mflux_nd_temp(n))
          call multifab_fill_boundary(mflux_nd_temp(n))
          call multifab_fill_boundary(mflux_cc_temp(n))
          
          if(filtering_width>0) then
             call multifab_filter(mflux_nd_temp(n), filtering_width)
             call multifab_filter(mflux_cc_temp(n), filtering_width)
             call multifab_fill_boundary(mflux_cc_temp(n)) ! First ghost cell is used in divergence
          end if

       else if (dm .eq. 3) then

          ng_e = mflux_ed_temp(1,1)%ng
          ng_w = eta_ed(1,1)%ng

          do i=1,nfabs(mflux_cc_temp(n))
             fp  => dataptr(mflux_cc_temp(n),i)
             fxp => dataptr(mflux_ed_temp(n,1),i)
             fyp => dataptr(mflux_ed_temp(n,2),i)
             fzp => dataptr(mflux_ed_temp(n,3),i)
             dp => dataptr(eta(n),i)
             lo = lwb(get_box(mflux_cc_temp(n),i))
             hi = upb(get_box(mflux_cc_temp(n),i))
             ! multiply by variance
             fp  = variance*fp
             fxp = variance*fxp
             fyp = variance*fyp
             fzp = variance*fzp
             ! if eta varies in space, multiply pointwise by sqrt(eta)
             if (visc_type < 0) then
                ep1 => dataptr(eta_ed(n,1),i)
                ep2 => dataptr(eta_ed(n,2),i)
                ep3 => dataptr(eta_ed(n,3),i)
                call mult_by_sqrt_eta_3d(fp(:,:,:,:),ng_c, &
                                         fxp(:,:,:,:),fyp(:,:,:,:),fzp(:,:,:,:),ng_e, &
                                         dp(:,:,:,1),ng_y, &
                                         ep1(:,:,:,1),ep2(:,:,:,1),ep3(:,:,:,1),ng_w,lo,hi)
             end if
             ! apply boundary conditions
             call mflux_bc_3d(fxp(:,:,:,:),fyp(:,:,:,:),fzp(:,:,:,:),ng_e,lo,hi, &
                              the_bc_level(n)%phys_bc_level_array(i,:,:))
          end do
          
          ! sync up random numbers at boundaries and ghost cells
          call multifab_internal_sync(mflux_ed_temp(n,1))
          call multifab_internal_sync(mflux_ed_temp(n,2))
          call multifab_internal_sync(mflux_ed_temp(n,3))
          call multifab_fill_boundary(mflux_ed_temp(n,1))
          call multifab_fill_boundary(mflux_ed_temp(n,2))
          call multifab_fill_boundary(mflux_ed_temp(n,3))
          call multifab_fill_boundary(mflux_cc_temp(n))

          if(filtering_width>0) then
             call multifab_filter(mflux_ed_temp(n,1), filtering_width)
             call multifab_filter(mflux_ed_temp(n,2), filtering_width)
             call multifab_filter(mflux_ed_temp(n,3), filtering_width)
             call multifab_filter(mflux_cc_temp(n)  , filtering_width)
             call multifab_fill_boundary(mflux_cc_temp(n)) ! First ghost cell is used in divergence
          end if

       end if

       ! calculate divergence and add to stoch_m_force
       do i=1,nfabs(stoch_m_force(n,1))
          fp => dataptr(mflux_cc_temp(n), i)
          dxp => dataptr(stoch_m_force(n,1),i)
          dyp => dataptr(stoch_m_force(n,2),i)
          lo =  lwb(get_box(stoch_m_force(n,1), i))
          hi =  upb(get_box(stoch_m_force(n,1), i))
          select case (dm)
          case (2)
             sp => dataptr(mflux_nd_temp(n), i)
             call stoch_m_force_2d(fp(:,:,1,:), sp(:,:,1,:), dxp(:,:,1,1), dyp(:,:,1,1), &
                                   ng_c, ng_n, ng_f, dx(n,:), lo, hi)
          case (3)
             dzp => dataptr(stoch_m_force(n,3), i)
             fxp => dataptr(mflux_ed_temp(n,1), i)
             fyp => dataptr(mflux_ed_temp(n,2), i)
             fzp => dataptr(mflux_ed_temp(n,3), i)
             call stoch_m_force_3d(fp(:,:,:,:), fxp(:,:,:,:), fyp(:,:,:,:), fzp(:,:,:,:), &
                                   dxp(:,:,:,1), dyp(:,:,:,1), dzp(:,:,:,1), &
                                   ng_c, ng_e, ng_f, dx(n,:), lo, hi)

          end select
       end do

       do i=1,dm
          call multifab_physbc_domainvel(stoch_m_force(n,i),vel_bc_comp+i-1, &
                                         the_bc_level(n),dx(n,:))
       end do

    end do

    do n=1,nlevs
       call multifab_destroy(mflux_cc_temp(n))
       if (dm .eq. 2) then
          call multifab_destroy(mflux_nd_temp(n))
       else if (dm .eq. 3) then
          call multifab_destroy(mflux_ed_temp(n,1))
          call multifab_destroy(mflux_ed_temp(n,2))
          call multifab_destroy(mflux_ed_temp(n,3))
       end if
    end do

  contains
    
    subroutine mult_by_sqrt_eta_2d(mflux_cc,ng_c,mflux_nd,ng_n,eta,ng_y,eta_nodal,ng_w,lo,hi)

      integer        , intent(in   ) :: lo(:),hi(:),ng_c,ng_n,ng_y,ng_w
      real(kind=dp_t), intent(inout) ::  mflux_cc(lo(1)-ng_c:,lo(2)-ng_c:,:)
      real(kind=dp_t), intent(inout) ::  mflux_nd(lo(1)-ng_n:,lo(2)-ng_n:,:)
      real(kind=dp_t), intent(in   ) ::       eta(lo(1)-ng_y:,lo(2)-ng_y:)
      real(kind=dp_t), intent(in   ) :: eta_nodal(lo(1)-ng_w:,lo(2)-ng_w:)

      ! local
      integer i,j

      do j=lo(2)-1,hi(2)+1
         do i=lo(1)-1,hi(1)+1
            mflux_cc(i,j,:) = mflux_cc(i,j,:) * sqrt(eta(i,j))
         end do
      end do

      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)+1
            mflux_nd(i,j,:) = mflux_nd(i,j,:) * sqrt(eta_nodal(i,j))
         end do
      end do

    end subroutine mult_by_sqrt_eta_2d

    subroutine mult_by_sqrt_eta_3d(mflux_cc,ng_c,mflux_xy,mflux_xz,mflux_yz,ng_e,eta,ng_y, &
         eta_xy,eta_xz,eta_yz,ng_w,lo,hi)

      integer        , intent(in   ) :: lo(:),hi(:),ng_c,ng_e,ng_y,ng_w
      real(kind=dp_t), intent(inout) :: mflux_cc(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:,:)
      real(kind=dp_t), intent(inout) :: mflux_xy(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
      real(kind=dp_t), intent(inout) :: mflux_xz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
      real(kind=dp_t), intent(inout) :: mflux_yz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
      real(kind=dp_t), intent(in   ) ::      eta(lo(1)-ng_y:,lo(2)-ng_y:,lo(3)-ng_y:)
      real(kind=dp_t), intent(in   ) ::   eta_xy(lo(1)-ng_w:,lo(2)-ng_w:,lo(3)-ng_w:)
      real(kind=dp_t), intent(in   ) ::   eta_xz(lo(1)-ng_w:,lo(2)-ng_w:,lo(3)-ng_w:)
      real(kind=dp_t), intent(in   ) ::   eta_yz(lo(1)-ng_w:,lo(2)-ng_w:,lo(3)-ng_w:)

      ! local
      integer i,j,k

      do k=lo(3)-1,hi(3)+1
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
               mflux_cc(i,j,k,:) = mflux_cc(i,j,k,:) * sqrt(eta(i,j,k))
            end do
         end do
      end do

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)+1
               mflux_xy(i,j,k,:) = mflux_xy(i,j,k,:) * sqrt(eta_xy(i,j,k))
            end do
         end do
      end do

      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)+1
               mflux_xz(i,j,k,:) = mflux_xz(i,j,k,:) * sqrt(eta_xz(i,j,k))
            end do
         end do
      end do

      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)
               mflux_yz(i,j,k,:) = mflux_yz(i,j,k,:) * sqrt(eta_yz(i,j,k))
            end do
         end do
      end do

    end subroutine mult_by_sqrt_eta_3d

    subroutine mflux_bc_2d(mflux_nd,ng_n,lo,hi,phys_bc)

      integer        , intent(in   ) :: lo(:),hi(:),ng_n
      real(kind=dp_t), intent(inout) :: mflux_nd(lo(1)-ng_n:,lo(2)-ng_n:,:)
      integer        , intent(in   ) :: phys_bc(:,:)

      ! y-mom fluxes that live on x-domain boundaries
      if (phys_bc(1,1) .eq. NO_SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_RESERVOIR) then
         mflux_nd(lo(1),lo(2):hi(2)+1,:) = sqrt(2.d0)*mflux_nd(lo(1),lo(2):hi(2)+1,:)
      else if (phys_bc(1,1) .eq. SLIP_WALL) then
         mflux_nd(lo(1),lo(2):hi(2)+1,:) = 0.d0
      end if
      if (phys_bc(1,2) .eq. NO_SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_RESERVOIR) then
         mflux_nd(hi(1)+1,lo(2):hi(2)+1,:) = sqrt(2.d0)*mflux_nd(hi(1)+1,lo(2):hi(2)+1,:)
      else if (phys_bc(1,2) .eq. SLIP_WALL) then
         mflux_nd(hi(1)+1,lo(2):hi(2)+1,:) = 0.d0
      end if

      ! x-mom fluxes that live on y-domain boundaries
      if (phys_bc(2,1) .eq. NO_SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_RESERVOIR) then
         mflux_nd(lo(1):hi(1)+1,lo(2),:) = sqrt(2.d0)*mflux_nd(lo(1):hi(1)+1,lo(2),:)
      else if (phys_bc(2,1) .eq. SLIP_WALL) then
         mflux_nd(lo(1):hi(1)+1,lo(2),:) = 0.d0
      end if
      if (phys_bc(2,2) .eq. NO_SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_RESERVOIR) then
         mflux_nd(lo(1):hi(1)+1,hi(2)+1,:) = sqrt(2.d0)*mflux_nd(lo(1):hi(1)+1,hi(2)+1,:)
      else if (phys_bc(2,2) .eq. SLIP_WALL) then
         mflux_nd(lo(1):hi(1)+1,hi(2)+1,:) = 0.d0
      end if

    end subroutine mflux_bc_2d

    subroutine mflux_bc_3d(mflux_xy,mflux_xz,mflux_yz,ng_e,lo,hi,phys_bc)

      integer        , intent(in   ) :: lo(:),hi(:),ng_e
      real(kind=dp_t), intent(inout) :: mflux_xy(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
      real(kind=dp_t), intent(inout) :: mflux_xz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
      real(kind=dp_t), intent(inout) :: mflux_yz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
      integer        , intent(in   ) :: phys_bc(:,:)

      ! y-mom and z-mom fluxes that live on x-domain boundaries
      if (phys_bc(1,1) .eq. NO_SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_RESERVOIR) then
         mflux_xy(lo(1),lo(2):hi(2)+1,lo(3):hi(3),:) = sqrt(2.d0)*mflux_xy(lo(1),lo(2):hi(2)+1,lo(3):hi(3),:)
         mflux_xz(lo(1),lo(2):hi(2),lo(3):hi(3)+1,:) = sqrt(2.d0)*mflux_xz(lo(1),lo(2):hi(2),lo(3):hi(3)+1,:)
      else if (phys_bc(1,1) .eq. SLIP_WALL) then
         mflux_xy(lo(1),lo(2):hi(2)+1,lo(3):hi(3),:) = 0.d0
         mflux_xz(lo(1),lo(2):hi(2),lo(3):hi(3)+1,:) = 0.d0
      end if
      if (phys_bc(1,2) .eq. NO_SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_RESERVOIR) then
         mflux_xy(hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3),:) = sqrt(2.d0)*mflux_xy(hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3),:)
         mflux_xz(hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1,:) = sqrt(2.d0)*mflux_xz(hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1,:)
      else if (phys_bc(1,2) .eq. SLIP_WALL) then
         mflux_xy(hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3),:) = 0.d0
         mflux_xz(hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1,:) = 0.d0
      end if

      ! x-mom and z-mom fluxes that live on y-domain boundaries
      if (phys_bc(2,1) .eq. NO_SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_RESERVOIR) then
         mflux_xy(lo(1):hi(1)+1,lo(2),lo(3):hi(3),:) = sqrt(2.d0)*mflux_xy(lo(1):hi(1)+1,lo(2),lo(3):hi(3),:)
         mflux_yz(lo(1):hi(1),lo(2),lo(3):hi(3)+1,:) = sqrt(2.d0)*mflux_yz(lo(1):hi(1),lo(2),lo(3):hi(3)+1,:)
      else if (phys_bc(2,1) .eq. SLIP_WALL) then
         mflux_xy(lo(1):hi(1)+1,lo(2),lo(3):hi(3),:) = 0.d0
         mflux_yz(lo(1):hi(1),lo(2),lo(3):hi(3)+1,:) = 0.d0
      end if
      if (phys_bc(2,2) .eq. NO_SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_RESERVOIR) then
         mflux_xy(lo(1):hi(1)+1,hi(2)+1,lo(3):hi(3),:) = sqrt(2.d0)*mflux_xy(lo(1):hi(1)+1,hi(2)+1,lo(3):hi(3),:)
         mflux_yz(lo(1):hi(1),hi(2)+1,lo(3):hi(3)+1,:) = sqrt(2.d0)*mflux_yz(lo(1):hi(1),hi(2)+1,lo(3):hi(3)+1,:)
      else if (phys_bc(2,2) .eq. SLIP_WALL) then
         mflux_xy(lo(1):hi(1)+1,hi(2)+1,lo(3):hi(3),:) = 0.d0
         mflux_yz(lo(1):hi(1),hi(2)+1,lo(3):hi(3)+1,:) = 0.d0
      end if

      ! x-mom and y-mom fluxes that live on z-domain boundaries
      if (phys_bc(3,1) .eq. NO_SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_RESERVOIR) then
         mflux_xz(lo(1):hi(1)+1,lo(2):hi(2),lo(3),:) = sqrt(2.d0)*mflux_xz(lo(1):hi(1)+1,lo(2):hi(2),lo(3),:)
         mflux_yz(lo(1):hi(1),lo(2):hi(2)+1,lo(3),:) = sqrt(2.d0)*mflux_yz(lo(1):hi(1),lo(2):hi(2)+1,lo(3),:)
      else if (phys_bc(3,1) .eq. SLIP_WALL) then
         mflux_xz(lo(1):hi(1)+1,lo(2):hi(2),lo(3),:) = 0.d0
         mflux_yz(lo(1):hi(1),lo(2):hi(2)+1,lo(3),:) = 0.d0
      end if
      if (phys_bc(3,2) .eq. NO_SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_RESERVOIR) then
         mflux_xz(lo(1):hi(1)+1,lo(2):hi(2),hi(3)+1,:) = sqrt(2.d0)*mflux_xz(lo(1):hi(1)+1,lo(2):hi(2),hi(3)+1,:)
         mflux_yz(lo(1):hi(1),lo(2):hi(2)+1,hi(3)+1,:) = sqrt(2.d0)*mflux_yz(lo(1):hi(1),lo(2):hi(2)+1,hi(3)+1,:)
      else if (phys_bc(3,2) .eq. SLIP_WALL) then
         mflux_xz(lo(1):hi(1)+1,lo(2):hi(2),hi(3)+1,:) = 0.d0
         mflux_yz(lo(1):hi(1),lo(2):hi(2)+1,hi(3)+1,:) = 0.d0
      end if

    end subroutine mflux_bc_3d

    subroutine stoch_m_force_2d(flux_cc,flux_nd,divx,divy,ng_c,ng_n,ng_f,dx,lo,hi)

      integer        , intent(in   ) :: lo(:),hi(:),ng_c,ng_n,ng_f
      real(kind=dp_t), intent(in   ) :: flux_cc(lo(1)-ng_c:,lo(2)-ng_c:,:)
      real(kind=dp_t), intent(in   ) :: flux_nd(lo(1)-ng_n:,lo(2)-ng_n:,:)
      real(kind=dp_t), intent(inout) ::    divx(lo(1)-ng_f:,lo(2)-ng_f:)
      real(kind=dp_t), intent(inout) ::    divy(lo(1)-ng_f:,lo(2)-ng_f:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      integer :: i,j

      ! divergence on x-faces
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1

            divx(i,j) = divx(i,j) + &
                 (flux_cc(i,j,1) - flux_cc(i-1,j,1)) / dx(1) + &
                 (flux_nd(i,j+1,1) - flux_nd(i,j,1)) / dx(2)

         end do
      end do

      ! divergence on y-faces
      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)

            divy(i,j) = divy(i,j) + &
                 (flux_nd(i+1,j,2) - flux_nd(i,j,2)) / dx(1) + &
                 (flux_cc(i,j,2) - flux_cc(i,j-1,2)) / dx(2)

         end do
      end do


    end subroutine stoch_m_force_2d

    subroutine stoch_m_force_3d(flux_cc,flux_xy,flux_xz,flux_yz,divx,divy,divz, &
                                ng_c,ng_e,ng_f,dx,lo,hi)

      integer        , intent(in   ) :: lo(:),hi(:),ng_c,ng_e,ng_f
      real(kind=dp_t), intent(in   ) :: flux_cc(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:,:)
      real(kind=dp_t), intent(in   ) :: flux_xy(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
      real(kind=dp_t), intent(in   ) :: flux_xz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
      real(kind=dp_t), intent(in   ) :: flux_yz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
      real(kind=dp_t), intent(inout) ::    divx(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
      real(kind=dp_t), intent(inout) ::    divy(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
      real(kind=dp_t), intent(inout) ::    divz(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      ! local
      integer :: i,j,k

      ! divergence on x-faces
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)+1

               divx(i,j,k) = divx(i,j,k) + &
                    (flux_cc(i,j,k,1) - flux_cc(i-1,j,k,1)) / dx(1) + &
                    (flux_xy(i,j+1,k,1) - flux_xy(i,j,k,1)) / dx(2) + &
                    (flux_xz(i,j,k+1,1) - flux_xz(i,j,k,1)) / dx(3)

            end do
         end do
      end do

      ! divergence on y-faces
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)

               divy(i,j,k) = divy(i,j,k) + &
                    (flux_xy(i+1,j,k,2) - flux_xy(i,j,k,2)) / dx(1) + &
                    (flux_cc(i,j,k,2) - flux_cc(i,j-1,k,2)) / dx(2) + &
                    (flux_yz(i,j,k+1,1) - flux_yz(i,j,k,1)) / dx(3)

            end do
         end do
      end do

      ! divergence on z-faces
      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

               divz(i,j,k) = divz(i,j,k) + &
                    (flux_xz(i+1,j,k,2) - flux_xz(i,j,k,2)) / dx(1) + &
                    (flux_yz(i,j+1,k,2) - flux_yz(i,j,k,2)) / dx(2) + &
                    (flux_cc(i,j,k,3) - flux_cc(i,j,k-1,3)) / dx(3)

            end do
         end do
      end do

    end subroutine stoch_m_force_3d

  end subroutine stochastic_m_fluxdiv

  ! fill the stochastic multifabs with random numbers
  subroutine fill_m_stochastic(mla)

    type(ml_layout), intent(in   ) :: mla

    ! local
    integer :: n,nlevs,dm,comp

    nlevs = mla%nlevel
    dm = mla%dim

    do comp=1,n_rngs

       ! Diagonal components of stochastic stress tensor for momentum
       select case(stoch_stress_form)
       case (0) ! Non-symmetric
          call multifab_fill_random(mflux_cc(:,comp))
       case default ! Symmetric
          call multifab_fill_random(mflux_cc(:,comp), variance=2.d0)
       end select

       ! Off-diagonal components of stochastic stress tensor for momentum
       if (dm .eq. 2) then
          ! in 2D, we need 2 random fluxes at each node
          select case(stoch_stress_form)
          case(0) ! Non-symmetric
             call multifab_fill_random(mflux_nd(:,comp))
          case default ! Symmetric
             call multifab_fill_random(mflux_nd(:,comp), comp=1)
             do n=1,nlevs
                call multifab_copy_c(mflux_nd(n,comp),2,mflux_nd(n,comp),1)
             end do
          end select
       else if (dm .eq. 3) then
          ! in 3D, we need 2 random fluxes at each edge
          select case(stoch_stress_form)
          case(0) ! Non-symmetric
             call multifab_fill_random(mflux_ed(:,1,comp))
             call multifab_fill_random(mflux_ed(:,2,comp))
             call multifab_fill_random(mflux_ed(:,3,comp))
          case default ! Symmetric
             call multifab_fill_random(mflux_ed(:,1,comp), comp=1)
             call multifab_fill_random(mflux_ed(:,2,comp), comp=1)
             call multifab_fill_random(mflux_ed(:,3,comp), comp=1)
             do n = 1, nlevs
                call multifab_copy_c(mflux_ed(n,1,comp),2,mflux_ed(n,1,comp),1)
                call multifab_copy_c(mflux_ed(n,2,comp),2,mflux_ed(n,2,comp),1)
                call multifab_copy_c(mflux_ed(n,3,comp),2,mflux_ed(n,3,comp),1)
             end do
          end select
       end if

    end do

  end subroutine fill_m_stochastic

  ! call this once at the beginning of simulation to allocate multifabs
  ! that will hold random numbers
  subroutine init_m_stochastic(mla,n_rngs_in)

    type(ml_layout), intent(in   ) :: mla
    integer        , intent(in   ) :: n_rngs_in

    ! local
    integer :: n,nlevs,dm,comp
    logical :: nodal_temp(mla%dim)
    
    n_rngs = n_rngs_in

    nlevs = mla%nlevel
    dm = mla%dim

    allocate(mflux_cc(mla%nlevel         , n_rngs))
    allocate(mflux_nd(mla%nlevel         , n_rngs))
    allocate(mflux_ed(mla%nlevel, 3      , n_rngs))

    do n=1,nlevs
       do comp=1,n_rngs
          ! we need dm cell-centered fluxes for momentum
          call multifab_build(mflux_cc(n,comp),mla%la(n),dm,max(1,filtering_width))
          if (dm .eq. 2) then
             ! in 2D, we need 2 random fluxes at each node
             nodal_temp = .true.
             call multifab_build(mflux_nd(n,comp),mla%la(n),2,filtering_width,nodal_temp)
          else if (dm .eq. 3) then
             ! in 3D, we need 2 random fluxes at each edge
             nodal_temp(1) = .true.
             nodal_temp(2) = .true.
             nodal_temp(3) = .false.
             call multifab_build(mflux_ed(n,1,comp),mla%la(n),2,filtering_width,nodal_temp)
             nodal_temp(1) = .true.
             nodal_temp(2) = .false.
             nodal_temp(3) = .true.
             call multifab_build(mflux_ed(n,2,comp),mla%la(n),2,filtering_width,nodal_temp)
             nodal_temp(1) = .false.
             nodal_temp(2) = .true.
             nodal_temp(3) = .true.
             call multifab_build(mflux_ed(n,3,comp),mla%la(n),2,filtering_width,nodal_temp)
          end if
       end do ! end loop over n_rngs
    end do ! end loop over nlevs

  end subroutine init_m_stochastic

  ! call this once at the end of simulation to deallocate memory
  subroutine destroy_m_stochastic(mla)

    type(ml_layout), intent(in   ) :: mla

    ! local
    integer :: n,nlevs,dm,comp
    
    nlevs = mla%nlevel
    dm = mla%dim

    do n=1,nlevs
       do comp=1,n_rngs
          call multifab_destroy(mflux_cc(n,comp))
          if (dm .eq. 2) then
             call multifab_destroy(mflux_nd(n,comp))
          else if (dm .eq. 3) then
             call multifab_destroy(mflux_ed(n,1,comp))
             call multifab_destroy(mflux_ed(n,2,comp))
             call multifab_destroy(mflux_ed(n,3,comp))
          end if
       end do
    end do
    
    deallocate(mflux_cc,mflux_nd,mflux_ed)

  end subroutine destroy_m_stochastic

 ! Add equilibrium fluctuations to the momentum (valid and ghost regions)
 subroutine add_m_fluctuations(mla,dx,variance,s_cc,s_face,m_face)

   type(ml_layout), intent(in   ) :: mla
   real(dp_t)     , intent(in   ) :: variance, dx(:,:)
   type(multifab) , intent(in   ) :: s_cc(:), s_face(:,:)
   type(multifab) , intent(inout) :: m_face(:,:)

   ! local
   type(multifab) :: mactemp(mla%nlevel,mla%dim)
   integer :: n,i,dm,nlevs
   real(dp_t) :: av_mom(mla%dim)

   dm = mla%dim
   nlevs = mla%nlevel

   do n=1,nlevs
      do i=1,dm
         call multifab_build_edge(mactemp(n,i),mla%la(n),1,1,i)
      end do
   end do

   ! Generate random numbers first and store them in mactemp
   do n=1,nlevs
      do i=1,dm
         call multifab_fill_random(mactemp(n:n,i), &
              variance=abs(variance)*k_B*temperature/product(dx(n,1:dm)), variance_mfab=s_face(n:n,i))
         call saxpy(m_face(n,i), 1.0_dp_t, mactemp(n,i), all=.true.)
      end do
   end do

   do n=1,nlevs
      do i=1,dm
         ! We need to ensure periodic BCs are obeyed for the random values
         call multifab_internal_sync(m_face(n,i))
      end do
   enddo
   
   if(variance<0) then ! Ensure zero total momentum
      if (parallel_IOProcessor()) then
         write(*,"(A,100G17.9)") "Randomly INITIALized momenta"
      end if
      
      call sum_momenta(mla, m_face, av_mom)
      do i=1,dm
         call setval(mactemp(1,i), -av_mom(i))
         call saxpy(m_face(1,i), 1.0_dp_t, mactemp(1,i))
      end do
   end if

   do n=1,nlevs
      do i=1,dm
         call multifab_destroy(mactemp(n,i))
      end do
   end do

 end subroutine add_m_fluctuations

end module stochastic_m_fluxdiv_module
