module stochastic_rhoc_fluxdiv_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use BoxLibRNGs
  use sum_momenta_module
  use convert_stag_module
  use bc_module
  use multifab_physbc_stag_module
  use multifab_fill_random_module
  use multifab_filter_module
  use probin_common_module , only: visc_type, visc_coef, diff_type, variance_coef, k_B, &
                                   filtering_width, stoch_stress_form, diff_coef, &
                                   rhobar, molmass
  use probin_binarylm_module, only: conc_scal, temperature

  implicit none

  private

  public :: stochastic_rhoc_fluxdiv, fill_rhoc_stochastic, &
       init_rhoc_stochastic, destroy_rhoc_stochastic, add_c_fluctuations

  ! Stochastic fluxes for scalars are face-centered
  type(multifab), allocatable, save :: sflux_fc(:,:,:)
  
  logical, save :: warn_bad_c=.false. ! Should we issue warnings about c<0 or c>1

  integer, save :: n_rngs ! how many random number stages
  
contains

  ! Note that here we *increment* stoch_rhoc_force so it must be initialized externally!
  subroutine stochastic_rhoc_fluxdiv(mla,the_bc_level,stoch_rhoc_force,s_fc, &
                                     chi_fc,dx,dt,vel_bc_n,weights)
    
    type(ml_layout), intent(in   ) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(multifab) , intent(inout) :: stoch_rhoc_force(:)
    type(multifab) , intent(in   ) :: s_fc(:,:)
    type(multifab) , intent(in   ) :: chi_fc(:,:)
    real(dp_t)     , intent(in   ) :: dx(:,:),dt,weights(:)
    type(multifab) , intent(inout) :: vel_bc_n(:,:)

    ! local
    integer :: n,nlevs,i,dm,m,comp
    integer :: ng_x,ng_z,ng_s,ng_f,ng_b

    real(dp_t) :: variance

    type(multifab) :: sflux_fc_temp(mla%nlevel,mla%dim)

    real(kind=dp_t), pointer :: fp(:,:,:,:), sp(:,:,:,:), dp(:,:,:,:)
    real(kind=dp_t), pointer :: fxp(:,:,:,:), fyp(:,:,:,:), fzp(:,:,:,:)
    real(kind=dp_t), pointer :: bxp(:,:,:,:), byp(:,:,:,:), bzp(:,:,:,:)
    integer :: lo(mla%dim), hi(mla%dim)

    nlevs = mla%nlevel
    dm = mla%dim

    do n=1,nlevs

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! create a place to store temporary copy of the random numbers
       ! these copies will be scaled and used to create the stochastic flux divergence
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       do i=1,dm
          ! we need one face-centered flux for each concentration
          call multifab_build_edge(sflux_fc_temp(n,i),mla%la(n),1,sflux_fc(n,i,1)%ng,i)
          call multifab_setval(sflux_fc_temp(n,i),0.d0,all=.true.)
          ! add weighted contribution of fluxes
          do comp=1,n_rngs
             call saxpy(sflux_fc_temp(n,i),weights(comp),sflux_fc(n,i,comp),all=.true.)
          end do
       end do

    end do

    ng_x = sflux_fc_temp(1,1)%ng
    ng_z = chi_fc(1,1)%ng
    ng_s = s_fc(1,1)%ng
    ng_f = stoch_rhoc_force(1)%ng
    ng_b = vel_bc_n(1,1)%ng

    do n=1,nlevs

       if (diff_type < 0) then
          ! chi varies in space, add its contribution below in an i/j/k loop
          variance = sqrt(variance_coef*conc_scal*2.d0          /(product(dx(n,1:dm))*dt))
       else
          ! chi is constant in space, include it here
          variance = sqrt(variance_coef*conc_scal*2.d0*diff_coef/(product(dx(n,1:dm))*dt))
       end if

       do i=1,dm       
          do m=1, nfabs(sflux_fc_temp(n,i))
             fp => dataptr(sflux_fc_temp(n,i),m)
             sp => dataptr(s_fc(n,i),m)
             lo = lwb(get_box(sflux_fc_temp(n,i),m))
             hi = upb(get_box(sflux_fc_temp(n,i),m))
             ! multiply by variance
             fp = variance*fp
             ! If chi varies in space, multiply pointwise by sqrt(chi).
             ! Then, include multiplicative 
             ! sqrt(rho * mu_c^-1 k_b T) = sqrt(rho * c * (1-c) * M)
             ! scaling for random rho*c flux and set rho*c flux on walls to zero
             select case (dm)
             case (2)
                if (diff_type < 0) then
                   dp => dataptr(chi_fc(n,i),m)
                   call mult_by_sqrt_chi_2d(fp(:,:,1,:),ng_x,dp(:,:,1,1),ng_z,i,lo,hi)
                end if
                call scale_rhoc_2d(fp(:,:,1,:),ng_x,sp(:,:,1,1:),ng_s,i,lo,hi, &
                                   the_bc_level(n)%phys_bc_level_array(m,:,:))
             case (3)
                if (diff_type < 0) then
                   dp => dataptr(chi_fc(n,i),m)
                   call mult_by_sqrt_chi_3d(fp(:,:,:,:),ng_x,dp(:,:,:,1),ng_z,i,lo,hi)
                end if
                call scale_rhoc_3d(fp(:,:,:,:),ng_x,sp(:,:,:,1:),ng_s,i,lo,hi, &
                                   the_bc_level(n)%phys_bc_level_array(m,:,:))
             end select
          end do
             
          ! Now sync up the random numbers at the boundaries:
          call multifab_internal_sync(sflux_fc_temp(n,i))
          call multifab_fill_boundary(sflux_fc_temp(n,i))          

          if (filtering_width > 0) then
             call multifab_filter(sflux_fc_temp(n,i), filtering_width)
          end if

       end do

       ! calculate divergence and add to stoch_rhoc_force
       do i=1,nfabs(stoch_rhoc_force(n))
          fp  => dataptr(stoch_rhoc_force(n),i)
          fxp => dataptr(sflux_fc_temp(n,1), i)
          fyp => dataptr(sflux_fc_temp(n,2), i)
          bxp => dataptr(vel_bc_n(n,1), i)
          byp => dataptr(vel_bc_n(n,2), i)
          lo =  lwb(get_box(stoch_rhoc_force(n), i))
          hi =  upb(get_box(stoch_rhoc_force(n), i))
          select case (dm)
          case (2)
             call stoch_rhoc_force_2d(fxp(:,:,1,:), fyp(:,:,1,:), ng_x, &
                                   fp(:,:,1,1:), ng_f, &
                                   bxp(:,:,1,1), byp(:,:,1,1), ng_b, &
                                   dx(n,:),lo,hi, &
                                   the_bc_level(n)%adv_bc_level_array(i,:,:,:))
          case (3)
             fzp => dataptr(sflux_fc_temp(n,3), i)
             bzp => dataptr(vel_bc_n(n,3), i)
             call stoch_rhoc_force_3d(fxp(:,:,:,:), fyp(:,:,:,:), fzp(:,:,:,:), ng_x, &
                                   fp(:,:,:,1:), ng_f, &
                                   bxp(:,:,:,1), byp(:,:,:,1), bzp(:,:,:,1), ng_b, &
                                   dx(n,:),lo,hi, &
                                   the_bc_level(n)%adv_bc_level_array(i,:,:,:))
          end select
       end do

    end do

    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(sflux_fc_temp(n,i))
       end do
    end do

  contains
    
    subroutine mult_by_sqrt_chi_2d(sflux,ng_x,chi_face,ng_z,idim,lo,hi)

      integer        , intent(in   ) :: lo(:),hi(:),ng_x,ng_z,idim
      real(kind=dp_t), intent(inout) ::    sflux(lo(1)-ng_x:,lo(2)-ng_x:,:)
      real(kind=dp_t), intent(in   ) :: chi_face(lo(1)-ng_z:,lo(2)-ng_z:)

      ! local
      integer :: i,j,off(2)

      off = 0
      off(idim) = 1

      do i=lo(1),hi(1)+off(1)
         do j=lo(2),hi(2)+off(2)
            sflux(i,j,:) = sflux(i,j,:) * sqrt(chi_face(i,j))
         end do
      end do

    end subroutine mult_by_sqrt_chi_2d

    subroutine mult_by_sqrt_chi_3d(sflux,ng_x,chi_face,ng_z,idim,lo,hi)

      integer        , intent(in   ) :: lo(:),hi(:),ng_x,ng_z,idim
      real(kind=dp_t), intent(inout) ::    sflux(lo(1)-ng_x:,lo(2)-ng_x:,lo(3)-ng_x:,:)
      real(kind=dp_t), intent(in   ) :: chi_face(lo(1)-ng_z:,lo(2)-ng_z:,lo(3)-ng_z:)

      ! local
      integer :: i,j,k,off(3)

      off = 0
      off(idim) = 1

      do i=lo(1),hi(1)+off(1)
         do j=lo(2),hi(2)+off(2)
            do k=lo(3),hi(3)+off(3)
               sflux(i,j,k,:) = sflux(i,j,k,:) * sqrt(chi_face(i,j,k))
            end do
         end do
      end do

    end subroutine mult_by_sqrt_chi_3d

    subroutine scale_rhoc_2d(sflux,ng_x,s_fc,ng_s,idim,lo,hi,phys_bc)

      integer        , intent(in   ) :: lo(:),hi(:),ng_x,ng_s,idim
      integer        , intent(in   ) :: phys_bc(:,:)
      real(kind=dp_t), intent(inout) ::  sflux(lo(1)-ng_x:,lo(2)-ng_x:,:)
      real(kind=dp_t), intent(in   ) :: s_fc(lo(1)-ng_s:,lo(2)-ng_s:,:)

      integer :: i,j
      real(kind=dp_t) :: rho, fac, s_fac

      if (idim .eq. 1) then

         do j = lo(2),hi(2)
            do i = lo(1),hi(1)+1
               rho = s_fc(i,j,1) ! Density at face
               fac = s_fc(i,j,2)/s_fc(i,j,1) ! Concentration at face
               call concentration_amplitude(s_fac, fac, rho)
               if (s_fac .lt. 0.0d0) then
                  if(warn_bad_c) call bl_warn ('negative c*(1-c) set to zero')
                  s_fac = 0.0d0
               end if
               sflux(i,j,:) = sqrt(s_fac)*sflux(i,j,:)
            end do
         end do

         if (phys_bc(1,1) .eq. NO_SLIP_WALL .or. phys_bc(1,1) .eq. SLIP_WALL) then
            sflux(lo(1),lo(2):hi(2),:) = 0.d0
         end if

         if (phys_bc(1,2) .eq. NO_SLIP_WALL .or. phys_bc(1,2) .eq. SLIP_WALL) then
            sflux(hi(1)+1,lo(2):hi(2),:) = 0.d0
         end if

         if (phys_bc(1,1) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1),lo(2):hi(2),:) = sqrt(2.0d0)*sflux(lo(1),lo(2):hi(2),:)
         end if

         if (phys_bc(1,2) .eq. NO_SLIP_RESERVOIR) then
            sflux(hi(1)+1,lo(2):hi(2),:) = sqrt(2.0d0)*sflux(hi(1)+1,lo(2):hi(2),:)
         end if

      else

         do j = lo(2),hi(2)+1
            do i = lo(1),hi(1)
               rho = s_fc(i,j,1) ! Density at face
               fac = s_fc(i,j,2)/s_fc(i,j,1) ! Concentration at face
               call concentration_amplitude(s_fac, fac, rho)
               if (s_fac .lt. 0.0d0) then
                  if(warn_bad_c) call bl_warn ('negative c*(1-c) set to zero')
                  s_fac = 0.0d0
               end if
               sflux(i,j,:) = sqrt(s_fac)*sflux(i,j,:)
            end do
         end do

         if (phys_bc(2,1) .eq. NO_SLIP_WALL .or. phys_bc(2,1) .eq. SLIP_WALL) then
            sflux(lo(1):hi(1),lo(2),:) = 0.d0
         end if

         if (phys_bc(2,2) .eq. NO_SLIP_WALL .or. phys_bc(2,2) .eq. SLIP_WALL) then
            sflux(lo(1):hi(1),hi(2)+1,:) = 0.d0
         end if

         if (phys_bc(2,1) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1):hi(1),lo(2),:) = sqrt(2.0d0)*sflux(lo(1):hi(1),lo(2),:)
         end if

         if (phys_bc(2,2) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1):hi(1),hi(2)+1,:) = sqrt(2.0d0)*sflux(lo(1):hi(1),hi(2)+1,:)
         end if

      end if

    end subroutine scale_rhoc_2d

    subroutine scale_rhoc_3d(sflux,ng_x,s_fc,ng_s,idim,lo,hi,phys_bc)

      integer        , intent(in   ) :: lo(:),hi(:),ng_x,ng_s,idim
      integer        , intent(in   ) :: phys_bc(:,:)
      real(kind=dp_t), intent(inout) ::  sflux(lo(1)-ng_x:,lo(2)-ng_x:,lo(3)-ng_x:,:)
      real(kind=dp_t), intent(in   ) :: s_fc(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)

      integer :: i,j,k
      real(kind=dp_t) :: rho, fac, s_fac

      if (idim .eq. 1) then

         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)+1
                  rho = s_fc(i,j,k,1) ! Density at face
                  fac = s_fc(i,j,k,2)/s_fc(i,j,k,1) ! Concentration at face
                  call concentration_amplitude(s_fac, fac, rho)    
                  if (s_fac .lt. 0.0d0) then
                     if(warn_bad_c) call bl_warn ('negative c*(1-c) set to zero')
                     s_fac = 0.0d0
                  end if
                  sflux(i,j,k,:) = sqrt(s_fac)*sflux(i,j,k,:)
               end do
            end do
         end do

         if (phys_bc(1,1) .eq. NO_SLIP_WALL .or. phys_bc(1,1) .eq. SLIP_WALL) then
            sflux(lo(1),lo(2):hi(2),lo(3):hi(3),:) = 0.d0
         end if

         if (phys_bc(1,2) .eq. NO_SLIP_WALL .or. phys_bc(1,2) .eq. SLIP_WALL) then
            sflux(hi(1)+1,lo(2):hi(2),lo(3):hi(3),:) = 0.d0
         end if

         if (phys_bc(1,1) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1),lo(2):hi(2),lo(3):hi(3),:) = sqrt(2.0d0)*sflux(lo(1),lo(2):hi(2),lo(3):hi(3),:)
         end if

         if (phys_bc(1,2) .eq. NO_SLIP_RESERVOIR) then
            sflux(hi(1)+1,lo(2):hi(2),lo(3):hi(3),:) = sqrt(2.0d0)*sflux(hi(1)+1,lo(2):hi(2),lo(3):hi(3),:)
         end if

      else if (idim .eq. 2) then

         do k = lo(3),hi(3)
            do j = lo(2),hi(2)+1
               do i = lo(1),hi(1)
                  rho = s_fc(i,j,k,1) ! Density at face
                  fac = s_fc(i,j,k,2)/s_fc(i,j,k,1) ! Concentration at face
                  call concentration_amplitude(s_fac, fac, rho)
                  if (s_fac .lt. 0.0d0) then
                     if(warn_bad_c) call bl_warn ('negative c*(1-c) set to zero')
                     s_fac = 0.0d0
                  end if
                  sflux(i,j,k,:) = sqrt(s_fac)*sflux(i,j,k,:)
               end do
            end do
         end do

         if (phys_bc(2,1) .eq. NO_SLIP_WALL .or. phys_bc(2,1) .eq. SLIP_WALL) then
            sflux(lo(1):hi(1),lo(2),lo(3):hi(3),:) = 0.d0
         end if

         if (phys_bc(2,2) .eq. NO_SLIP_WALL .or. phys_bc(2,2) .eq. SLIP_WALL) then
            sflux(lo(1):hi(1),hi(2)+1,lo(3):hi(3),:) = 0.d0
         end if

         if (phys_bc(2,1) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1):hi(1),lo(2),lo(3):hi(3),:) = sqrt(2.0d0)*sflux(lo(1):hi(1),lo(2),lo(3):hi(3),:)
         end if

         if (phys_bc(2,2) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1):hi(1),hi(2)+1,lo(3):hi(3),:) = sqrt(2.0d0)*sflux(lo(1):hi(1),hi(2)+1,lo(3):hi(3),:)
         end if

      else

         do k = lo(3),hi(3)+1
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  rho = s_fc(i,j,k,1) ! Density at face
                  fac = s_fc(i,j,k,2)/s_fc(i,j,k,1) ! Concentration at face
                  call concentration_amplitude(s_fac, fac, rho)
                  if (s_fac .lt. 0.0d0) then
                     if(warn_bad_c) call bl_warn ('negative c*(1-c) set to zero')
                     s_fac = 0.0d0
                  end if
                  sflux(i,j,k,:) = sqrt(s_fac)*sflux(i,j,k,:)
               end do
            end do
         end do

         if (phys_bc(3,1) .eq. NO_SLIP_WALL .or. phys_bc(3,1) .eq. SLIP_WALL) then
            sflux(lo(1):hi(1),lo(2):hi(2),lo(3),:) = 0.d0
         end if

         if (phys_bc(3,2) .eq. NO_SLIP_WALL .or. phys_bc(3,2) .eq. SLIP_WALL) then
            sflux(lo(1):hi(1),lo(2):hi(2),hi(3)+1,:) = 0.d0
         end if

         if (phys_bc(3,1) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1):hi(1),lo(2):hi(2),lo(3),:) = sqrt(2.0d0)*sflux(lo(1):hi(1),lo(2):hi(2),lo(3),:)
         end if

         if (phys_bc(3,2) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1):hi(1),lo(2):hi(2),hi(3)+1,:) = sqrt(2.0d0)*sflux(lo(1):hi(1),lo(2):hi(2),hi(3)+1,:)
         end if

      end if

    end subroutine scale_rhoc_3d

    subroutine concentration_amplitude(s_fac, fac, rho)

      real(kind=dp_t), intent(in) :: rho, fac
      real(kind=dp_t), intent(out) :: s_fac
      
      ! For ideal mixture use: rho*c*(1-c)*(c*m2+(1-c)*m1)
      s_fac = rho*fac*(1.0d0-fac)* (fac*molmass(2)+(1.0d0-fac)*molmass(1))  
    
    end subroutine concentration_amplitude
    
    subroutine stoch_rhoc_force_2d(xflux,yflux,ng_x,div,ng_f,vel_bc_nx,vel_bc_ny,ng_b,dx,lo,hi,adv_bc)

      integer        , intent(in   ) :: lo(:),hi(:),ng_f,ng_x,ng_b
      real(kind=dp_t), intent(in   ) ::     xflux(lo(1)-ng_x:,lo(2)-ng_x:,:)
      real(kind=dp_t), intent(in   ) ::     yflux(lo(1)-ng_x:,lo(2)-ng_x:,:)
      real(kind=dp_t), intent(inout) ::       div(lo(1)-ng_f:,lo(2)-ng_f:,:)
      real(kind=dp_t), intent(inout) :: vel_bc_nx(lo(1)-ng_b:,lo(2)-ng_b:)
      real(kind=dp_t), intent(inout) :: vel_bc_ny(lo(1)-ng_b:,lo(2)-ng_b:)
      real(kind=dp_t), intent(in   ) ::   dx(:)
      integer        , intent(in   ) :: adv_bc(:,:,:)

      integer :: i,j
      real(kind=dp_t) :: S_fac

      S_fac = 1.d0/rhobar(1)-1.d0/rhobar(2)

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            div(i,j,:) = div(i,j,:) + &
                 (xflux(i+1,j,:) - xflux(i,j,:)) / dx(1) + &
                 (yflux(i,j+1,:) - yflux(i,j,:)) / dx(2)
         end do
      end do
   
      ! update umac bc based on flux at boundary
      if (adv_bc(1,1,1) .eq. DIR_VEL .and. adv_bc(1,1,4) .eq. EXT_DIR) then
         vel_bc_nx(lo(1),lo(2):hi(2)) = vel_bc_nx(lo(1),lo(2):hi(2)) &
              + S_fac*xflux(lo(1),lo(2):hi(2),1)
      end if
      if (adv_bc(1,2,1) .eq. DIR_VEL .and. adv_bc(1,2,4) .eq. EXT_DIR) then
         vel_bc_nx(hi(1)+1,lo(2):hi(2)) = vel_bc_nx(hi(1)+1,lo(2):hi(2)) &
              + S_fac*xflux(hi(1)+1,lo(2):hi(2),1)
      end if

      ! update vmac bc based on flux at boundary
      if (adv_bc(2,1,2) .eq. DIR_VEL .and. adv_bc(2,1,4) .eq. EXT_DIR) then
         vel_bc_ny(lo(1):hi(1),lo(2)) = vel_bc_ny(lo(1):hi(1),lo(2)) + &
              S_fac*yflux(lo(1):hi(1),lo(2),1)
      end if
      if (adv_bc(2,2,2) .eq. DIR_VEL .and. adv_bc(2,2,4) .eq. EXT_DIR) then
         vel_bc_ny(lo(1):hi(1),hi(2)+1) = vel_bc_ny(lo(1):hi(1),hi(2)+1) + &
              S_fac*yflux(lo(1):hi(1),hi(2)+1,1)
      end if

    end subroutine stoch_rhoc_force_2d

    subroutine stoch_rhoc_force_3d(xflux,yflux,zflux,ng_x,div,ng_f, &
                                vel_bc_nx,vel_bc_ny,vel_bc_nz,ng_b,dx,lo,hi,adv_bc)

      integer        , intent(in   ) :: lo(:),hi(:),ng_f,ng_x,ng_b
      real(kind=dp_t), intent(in   ) ::     xflux(lo(1)-ng_x:,lo(2)-ng_x:,lo(3)-ng_x:,:)
      real(kind=dp_t), intent(in   ) ::     yflux(lo(1)-ng_x:,lo(2)-ng_x:,lo(3)-ng_x:,:)
      real(kind=dp_t), intent(in   ) ::     zflux(lo(1)-ng_x:,lo(2)-ng_x:,lo(3)-ng_x:,:)
      real(kind=dp_t), intent(inout) ::       div(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:,:)
      real(kind=dp_t), intent(inout) :: vel_bc_nx(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
      real(kind=dp_t), intent(inout) :: vel_bc_ny(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
      real(kind=dp_t), intent(inout) :: vel_bc_nz(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: adv_bc(:,:,:)

      integer :: i,j,k
      real(kind=dp_t) :: S_fac

      S_fac = 1.d0/rhobar(1)-1.d0/rhobar(2)

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               div(i,j,k,:) = div(i,j,k,:) + &
                    (xflux(i+1,j,k,:) - xflux(i,j,k,:)) / dx(1) + &
                    (yflux(i,j+1,k,:) - yflux(i,j,k,:)) / dx(2) + &
                    (zflux(i,j,k+1,:) - zflux(i,j,k,:)) / dx(3)
            end do
         end do
      end do

      ! update umac BC based on flux at boundary
      if (adv_bc(1,1,1) .eq. DIR_VEL .and. adv_bc(1,1,5) .eq. EXT_DIR) then
         vel_bc_nx(lo(1),lo(2):hi(2),lo(3):hi(3)) = vel_bc_nx(lo(1),lo(2):hi(2),lo(3):hi(3)) &
              + S_fac*xflux(lo(1),lo(2):hi(2),lo(3):hi(3),1)
      end if
      if (adv_bc(1,2,1) .eq. DIR_VEL .and. adv_bc(1,2,5) .eq. EXT_DIR) then
         vel_bc_nx(hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = vel_bc_nx(hi(1)+1,lo(2):hi(2),lo(3):hi(3)) &
              + S_fac*xflux(hi(1)+1,lo(2):hi(2),lo(3):hi(3),1)
      end if

      ! update vmac BC based on flux at boundary
      if (adv_bc(2,1,2) .eq. DIR_VEL .and. adv_bc(2,1,5) .eq. EXT_DIR) then
         vel_bc_ny(lo(1):hi(1),lo(2),lo(3):hi(3)) = vel_bc_ny(lo(1):hi(1),lo(2),lo(3):hi(3)) + &
              S_fac*yflux(lo(1):hi(1),lo(2),lo(3):hi(3),1)
      end if
      if (adv_bc(2,2,2) .eq. DIR_VEL .and. adv_bc(2,2,5) .eq. EXT_DIR) then
         vel_bc_ny(lo(1):hi(1),hi(2)+1,lo(3):hi(3)) = vel_bc_ny(lo(1):hi(1),hi(2)+1,lo(3):hi(3)) + &
              S_fac*yflux(lo(1):hi(1),hi(2)+1,lo(3):hi(3),1)
      end if

      ! update wmacBC  based on flux at boundary
      if (adv_bc(3,1,3) .eq. DIR_VEL .and. adv_bc(3,1,5) .eq. EXT_DIR) then
         vel_bc_nz(lo(1):hi(1),lo(2):hi(2),lo(3)) = vel_bc_nz(lo(1):hi(1),lo(2):hi(2),lo(3)) + &
              S_fac*zflux(lo(1):hi(1),lo(2):hi(2),lo(3),1)
      end if
      if (adv_bc(3,2,3) .eq. DIR_VEL .and. adv_bc(3,2,5) .eq. EXT_DIR) then
         vel_bc_nz(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = vel_bc_nz(lo(1):hi(1),lo(2):hi(2),hi(3)+1) + &
              S_fac*zflux(lo(1):hi(1),lo(2):hi(2),hi(3)+1,1)
      end if

    end subroutine stoch_rhoc_force_3d

  end subroutine stochastic_rhoc_fluxdiv

  ! fill the stochastic multifabs with random numbers
  subroutine fill_rhoc_stochastic(mla)

    type(ml_layout), intent(in   ) :: mla

    ! local
    integer :: nlevs,i,dm,comp

    nlevs = mla%nlevel
    dm = mla%dim

    do comp=1,n_rngs
       do i=1,dm
          call multifab_fill_random(sflux_fc(:,i,comp))
       end do
    end do

  end subroutine fill_rhoc_stochastic

  ! call this once at the beginning of simulation to allocate multifabs
  ! that will hold random numbers
  subroutine init_rhoc_stochastic(mla,n_rngs_in)

    type(ml_layout), intent(in   ) :: mla
    integer        , intent(in   ) :: n_rngs_in

    ! local
    integer :: n,nlevs,i,dm,comp
    
    n_rngs = n_rngs_in

    nlevs = mla%nlevel
    dm = mla%dim

    allocate(sflux_fc(mla%nlevel, mla%dim, n_rngs))

    do n=1,nlevs
       do comp=1,n_rngs
          do i=1,dm
             ! we need one face-centered flux for each concentration
             call multifab_build_edge(sflux_fc(n,i,comp),mla%la(n),1,filtering_width,i)
          end do
       end do ! end loop over n_rngs
    end do ! end loop over nlevs

  end subroutine init_rhoc_stochastic

  ! call this once at the end of simulation to deallocate memory
  subroutine destroy_rhoc_stochastic(mla)

    type(ml_layout), intent(in   ) :: mla

    ! local
    integer :: n,nlevs,i,dm,comp
    
    nlevs = mla%nlevel
    dm = mla%dim

    do n=1,nlevs
       do comp=1,n_rngs
          do i=1,dm
             call multifab_destroy(sflux_fc(n,i,comp))
          end do
       end do
    end do
    
    deallocate(sflux_fc)

  end subroutine destroy_rhoc_stochastic

  subroutine add_c_fluctuations(mla, dx, variance, prim, cctemp)
    type(ml_layout), intent(in   ) :: mla
    real(dp_t)     , intent(in   ) :: variance, dx(:,:)
    type(multifab) , intent(inout) :: prim(:)
    type(multifab) , intent(inout) :: cctemp(:)

    ! local
    integer :: n,dm,box,nlevs
    real(kind=dp_t), pointer :: fp(:,:,:,:), fpvar(:,:,:,:)

    dm = mla%dim
    nlevs = mla%nlevel
    
    do n=1,nlevs

       do box = 1, nfabs(prim(n))
          
          ! Only interior cells here:
          fp => dataptr(prim(n), box, get_ibox(prim(n),box))
          fpvar => dataptr(cctemp(n), box, get_ibox(cctemp(n),box))
          
          call NormalRNGs(fpvar(:,:,:,1), size(fpvar(:,:,:,1))) ! Fill the whole grid with random numbers

          ! Concentration fluctuations prefactor M*rho^(-1)*c*(1-c)
          fpvar(:,:,:,1) = sqrt( fp(:,:,:,2) * (1.0d0-fp(:,:,:,2)) / fp(:,:,:,1) * &
            (fp(:,:,:,2)*molmass(2)+(1.0d0-fp(:,:,:,2))*molmass(1)) ) * fpvar(:,:,:,1)
            
          ! Now add the fluctuations to the mean:  
          fp(:,:,:,2) = fp(:,:,:,2) + sqrt(abs(variance)/product(dx(n,1:dm))) * fpvar(:,:,:,1)
          ! And calculate the density from the EOS:
          fp(:,:,:,1) = 1.0d0/(fp(:,:,:,2)/rhobar(1)+(1.0d0-fp(:,:,:,2))/rhobar(2))
          
       end do   

    end do    
    
  end subroutine add_c_fluctuations

end module stochastic_rhoc_fluxdiv_module
