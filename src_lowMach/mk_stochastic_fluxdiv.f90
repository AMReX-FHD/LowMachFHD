module mk_stochastic_fluxdiv_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use BoxLibRNGs
  use analysis_module
  use convert_stag_module
  use bc_module
  use multifab_physbc_stag_module
  use probin_common_module , only: visc_type, diff_type
  use probin_lowmach_module, only: nscal, rhobar, visc_coef, diff_coef, variance_coef, &
                                   conc_scal, stoch_stress_form, filtering_width, mol_mass, &
                                   kT

  implicit none

  private

  public :: mk_stochastic_s_fluxdiv, mk_stochastic_m_fluxdiv, fill_stochastic, &
       init_stochastic, destroy_stochastic, add_momentum_fluctuations, &
       add_concentration_fluctuations

  ! Stochastic fluxes for momentum are generated on:
  ! -cell-centered grid for diagonal components
  ! -node-centered (2D) or edge-centered (3D) grid for off-diagonal components
  type(multifab), allocatable, save :: mflux_cc(:), mflux_nd(:), mflux_ed(:,:)

  ! Stochastic fluxes for scalars are face-centered
  type(multifab), allocatable, save :: sflux_fc(:,:)
  
  logical   , save :: warn_bad_c=.false. ! Should we issue warnings about c<0 or c>1
  
contains

  ! Note that here we *increment* stoch_s_force so it must be initialized externally!
  subroutine mk_stochastic_s_fluxdiv(mla,the_bc_level,stoch_s_force,s_fc, &
                                     chi_fc,dx,dt,vel_bc_n)
    
    type(ml_layout), intent(in   ) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(multifab) , intent(inout) :: stoch_s_force(:)
    type(multifab) , intent(in   ) :: s_fc(:,:)
    type(multifab) , intent(in   ) :: chi_fc(:,:)
    real(dp_t)     , intent(in   ) :: dx(:,:),dt
    type(multifab) , intent(inout) :: vel_bc_n(:,:)

    ! local
    integer :: n,nlevs,i,dm,m
    integer :: ng_x,ng_z,ng_s,ng_f,ng_b

    real(dp_t) :: variance

    type(multifab) :: sflux_fc_temp(mla%nlevel,mla%dim)

    real(kind=dp_t), pointer :: fp(:,:,:,:), sp(:,:,:,:), dp(:,:,:,:)
    real(kind=dp_t), pointer :: fxp(:,:,:,:), fyp(:,:,:,:), fzp(:,:,:,:)
    real(kind=dp_t), pointer :: bxp(:,:,:,:), byp(:,:,:,:), bzp(:,:,:,:)
    integer :: lo(mla%dim), hi(mla%dim)

    ! there are no scalars with stochastic terms
    if (nscal < 2) return

    nlevs = mla%nlevel
    dm = mla%dim

    do n=1,nlevs

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! create a place to store temporary copy of the random numbers
       ! these copies will be scaled and used to create the stochastic flux divergence
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       do i=1,dm
          ! we need one face-centered flux for each concentration
          call multifab_build_edge(sflux_fc_temp(n,i),mla%la(n),nscal-1,sflux_fc(n,i)%ng,i)
          ! make a copy of the random numbers
          call multifab_copy_c(sflux_fc_temp(n,i),1,sflux_fc(n,i),1,nscal-1,sflux_fc_temp(n,i)%ng)
       end do

    end do

    ng_x = sflux_fc_temp(1,1)%ng
    ng_z = chi_fc(1,1)%ng
    ng_s = s_fc(1,1)%ng
    ng_f = stoch_s_force(1)%ng
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
             call multifab_filter(sflux_fc_temp(n,i), dm)
          end if

       end do

       ! calculate divergence and add to stoch_s_force
       do i=1,nfabs(stoch_s_force(n))
          fp  => dataptr(stoch_s_force(n),i)
          fxp => dataptr(sflux_fc_temp(n,1), i)
          fyp => dataptr(sflux_fc_temp(n,2), i)
          bxp => dataptr(vel_bc_n(n,1), i)
          byp => dataptr(vel_bc_n(n,2), i)
          lo =  lwb(get_box(stoch_s_force(n), i))
          hi =  upb(get_box(stoch_s_force(n), i))
          select case (dm)
          case (2)
             call stoch_s_force_2d(fxp(:,:,1,:), fyp(:,:,1,:), ng_x, &
                                   fp(:,:,1,1:), ng_f, &
                                   bxp(:,:,1,1), byp(:,:,1,1), ng_b, &
                                   dx(n,:),lo,hi, &
                                   the_bc_level(n)%adv_bc_level_array(i,:,:,:))
          case (3)
             fzp => dataptr(sflux_fc_temp(n,3), i)
             bzp => dataptr(vel_bc_n(n,3), i)
             call stoch_s_force_3d(fxp(:,:,:,:), fyp(:,:,:,:), fzp(:,:,:,:), ng_x, &
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
      real(kind=dp_t) :: rho, fac(nscal-1), s_fac(nscal-1)

      if (idim .eq. 1) then

         do j = lo(2),hi(2)
            do i = lo(1),hi(1)+1
               rho = s_fc(i,j,1) ! Density at face
               fac(:) = s_fc(i,j,2:)/s_fc(i,j,1) ! Concentration at face
               call concentration_amplitude(s_fac, fac, rho)
               if (any (s_fac(:) .lt. 0.0d0)) then
                  if(warn_bad_c) call bl_warn ('negative c*(1-c) set to zero')
                  s_fac(:) = 0.0d0
               end if
               sflux(i,j,:) = sqrt(s_fac(:))*sflux(i,j,:)
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
               fac(:) = s_fc(i,j,2:)/s_fc(i,j,1) ! Concentration at face
               call concentration_amplitude(s_fac, fac, rho)
               if (any (s_fac(:) .lt. 0.0d0)) then
                  if(warn_bad_c) call bl_warn ('negative c*(1-c) set to zero')
                  s_fac(:) = 0.0d0
               end if
               sflux(i,j,:) = sqrt(s_fac(:))*sflux(i,j,:)
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
      real(kind=dp_t) :: rho, fac(nscal-1), s_fac(nscal-1)

      if (idim .eq. 1) then

         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)+1
                  rho = s_fc(i,j,k,1) ! Density at face
                  fac(:) = s_fc(i,j,k,2:)/s_fc(i,j,k,1) ! Concentration at face
                  call concentration_amplitude(s_fac, fac, rho)    
                  if (any (s_fac(:) .lt. 0.0d0)) then
                     if(warn_bad_c) call bl_warn ('negative c*(1-c) set to zero')
                     s_fac(:) = 0.0d0
                  end if
                  sflux(i,j,k,:) = sqrt(s_fac(:))*sflux(i,j,k,:)
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
                  fac(:) = s_fc(i,j,k,2:)/s_fc(i,j,k,1) ! Concentration at face
                  call concentration_amplitude(s_fac, fac, rho)
                  if (any (s_fac(:) .lt. 0.0d0)) then
                     if(warn_bad_c) call bl_warn ('negative c*(1-c) set to zero')
                     s_fac(:) = 0.0d0
                  end if
                  sflux(i,j,k,:) = sqrt(s_fac(:))*sflux(i,j,k,:)
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
                  fac(:) = s_fc(i,j,k,2:)/s_fc(i,j,k,1) ! Concentration at face
                  call concentration_amplitude(s_fac, fac, rho)
                  if (any (s_fac(:) .lt. 0.0d0)) then
                     if(warn_bad_c) call bl_warn ('negative c*(1-c) set to zero')
                     s_fac(:) = 0.0d0
                  end if
                  sflux(i,j,k,:) = sqrt(s_fac(:))*sflux(i,j,k,:)
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

      real(kind=dp_t), intent(in) :: rho, fac(nscal-1)
      real(kind=dp_t), intent(out) :: s_fac(nscal-1)
      
      ! For ideal mixture use: rho*c*(1-c)*M
      s_fac(:) = rho*fac(:)*(1.0d0-fac(:))* (fac(:)*mol_mass(1)+(1.0d0-fac(:))*mol_mass(2))  
    
    end subroutine concentration_amplitude
    
    subroutine stoch_s_force_2d(xflux,yflux,ng_x,div,ng_f,vel_bc_nx,vel_bc_ny,ng_b,dx,lo,hi,adv_bc)

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

    end subroutine stoch_s_force_2d

    subroutine stoch_s_force_3d(xflux,yflux,zflux,ng_x,div,ng_f, &
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

    end subroutine stoch_s_force_3d

  end subroutine mk_stochastic_s_fluxdiv

  ! Note that here we *increment* stoch_m_force so it must be initialized externally!
  subroutine mk_stochastic_m_fluxdiv(mla,the_bc_level,stoch_m_force,eta,eta_ed,dx,dt)
    
    type(ml_layout), intent(in   ) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(multifab) , intent(inout) :: stoch_m_force(:,:)
    type(multifab) , intent(in   ) :: eta(:)
    type(multifab) , intent(in   ) :: eta_ed(:,:)
    real(dp_t)     , intent(in   ) :: dx(:,:),dt

    ! local
    integer :: n,nlevs,dm,i
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
       
       ! make a copy of the random numbers
       call multifab_copy_c(mflux_cc_temp(n),1,mflux_cc(n),1,dm,mflux_cc_temp(n)%ng)
       if (dm .eq. 2) then
          call multifab_copy_c(mflux_nd_temp(n),1,mflux_nd(n),1,2,mflux_nd_temp(n)%ng)
       else if (dm .eq. 3) then
          call multifab_copy_c(mflux_ed_temp(n,1),1,mflux_ed(n,1),1,2,mflux_ed_temp(n,1)%ng)
          call multifab_copy_c(mflux_ed_temp(n,2),1,mflux_ed(n,2),1,2,mflux_ed_temp(n,2)%ng)
          call multifab_copy_c(mflux_ed_temp(n,3),1,mflux_ed(n,3),1,2,mflux_ed_temp(n,3)%ng)
       end if
    end do

    ng_c = mflux_cc_temp(1)%ng
    ng_y = eta(1)%ng
    ng_f = stoch_m_force(1,1)%ng

    do n=1,nlevs

       if (visc_type < 0) then
          ! eta varies in space, add its contribution below in an i/j/k loop
          variance = sqrt(variance_coef*2.d0*kT          /(product(dx(n,1:dm))*dt))
       else
          ! eta is constant in space, include it here
          variance = sqrt(variance_coef*2.d0*kT*visc_coef/(product(dx(n,1:dm))*dt))
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
             call multifab_filter(mflux_nd_temp(n), dm)
             call multifab_filter(mflux_cc_temp(n), dm)
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
             call multifab_filter(mflux_ed_temp(n,1), dm)
             call multifab_filter(mflux_ed_temp(n,2), dm)
             call multifab_filter(mflux_ed_temp(n,3), dm)
             call multifab_filter(mflux_cc_temp(n), dm)
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

  end subroutine mk_stochastic_m_fluxdiv

  ! fill the stochastic multifabs with random numbers
  subroutine fill_stochastic(mla)

    type(ml_layout), intent(in   ) :: mla

    ! local
    integer :: n,nlevs,i,dm

    nlevs = mla%nlevel
    dm = mla%dim

    ! Diagonal components of stochastic stress tensor for momentum
    select case(stoch_stress_form)
    case (0) ! Non-symmetric
       call multifab_fill_random(mflux_cc)
    case default ! Symmetric
       call multifab_fill_random(mflux_cc, variance=2.d0)
    end select

    ! Off-diagonal components of stochastic stress tensor for momentum
    if (dm .eq. 2) then
       ! in 2D, we need 2 random fluxes at each node
       select case(stoch_stress_form)
       case(0) ! Non-symmetric
          call multifab_fill_random(mflux_nd)
       case default ! Symmetric
          call multifab_fill_random(mflux_nd, comp=1)
          do n=1,nlevs
             call multifab_copy_c(mflux_nd(n),2,mflux_nd(n),1)
          end do
       end select
    else if (dm .eq. 3) then
       ! in 3D, we need 2 random fluxes at each edge
       select case(stoch_stress_form)
       case(0) ! Non-symmetric
          call multifab_fill_random(mflux_ed(:,1))
          call multifab_fill_random(mflux_ed(:,2))
          call multifab_fill_random(mflux_ed(:,3))
       case default ! Symmetric
          call multifab_fill_random(mflux_ed(:,1), comp=1)
          call multifab_fill_random(mflux_ed(:,2), comp=1)
          call multifab_fill_random(mflux_ed(:,3), comp=1)
          do n = 1, nlevs
             call multifab_copy_c(mflux_ed(n,1),2,mflux_ed(n,1),1)
             call multifab_copy_c(mflux_ed(n,2),2,mflux_ed(n,2),1)
             call multifab_copy_c(mflux_ed(n,3),2,mflux_ed(n,3),1)
          end do
       end select
    end if

    if(nscal>1) then ! Stochastic diffusive flux
       do i=1,dm
          call multifab_fill_random(sflux_fc(:,i))
       end do
    end if

  end subroutine fill_stochastic

  ! fill a multifab with random numbers
  subroutine multifab_fill_random(mfab, comp, variance, variance_mfab)
    type(multifab), intent(inout)           :: mfab(:)
    integer       , intent(in   ), optional :: comp ! Only one component
    real(dp_t)    , intent(in   ), optional :: variance
    type(multifab), intent(in   ), optional :: variance_mfab(:)

    integer :: n,box

    real(kind=dp_t), pointer :: fp(:,:,:,:), fpvar(:,:,:,:)

    !--------------------------------------
    do n=1,size(mfab)
       do box = 1, nfabs(mfab(n))
          if(present(comp)) then
             fp => dataptr(mfab(n),box,comp,1)
          else
             fp => dataptr(mfab(n),box)
          end if

          call NormalRNGs(fp, size(fp)) ! Fill the whole grid with random numbers

          if(present(variance_mfab)) then ! Must have same distribution
             fpvar => dataptr(variance_mfab(n),box,1,size(fp,4))
             fp=sqrt(fpvar)*fp
          end if
          if(present(variance)) then
             fp=sqrt(variance)*fp
          end if
       end do
    end do

  end subroutine multifab_fill_random

  subroutine multifab_filter(mfab, dm)

    ! Note: This routine assumes ghost values are valid on entry but does NOT 
    !       sync ghost values at the end.
    !       If ghost values need to remain consistent do a multifab_fill_boundary 
    !       after calling multifab_filter!

    type(multifab) , intent(inout) :: mfab
    integer, intent(in) :: dm

    integer :: i,j,k,box
    real(kind=dp_t), pointer :: fp(:,:,:,:), fpvar(:,:,:,:)

    integer :: lo(3), hi(3), lo_g(3), hi_g(3)
    
    lo=1
    hi=1
    lo_g=1
    hi_g=1
  
    !--------------------------------------
    FilterX: do box = 1, nfabs(mfab)

       fp => dataptr(mfab,box) ! Including ghosts
       
       ! Without ghosts:
       lo(1:dm) = lwb(get_ibox(mfab,box))
       hi(1:dm) = upb(get_ibox(mfab,box))
       !write(*,*) "lo=", lo, " hi=", hi
       
       ! With ghosts:
       lo_g(1:dm) = lwb(get_pbox(mfab,box))
       hi_g(1:dm) = upb(get_pbox(mfab,box))       
       !write(*,*) "lo_g=", lo_g, " hi_g=", hi_g

       if(.not.(all(lo_g(1:dm)<=lo(1:dm)-filtering_width) .and. all(hi_g(1:dm)>=hi(1:dm)+filtering_width) ) ) then
         call bl_error("Filtering requires at least filtering_width cells!")
       end if
       
       allocate(fpvar(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3), &
                      lbound(fp,4):ubound(fp,4))) ! Temporary array
              
       ! Filter along x:
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                select case(filtering_width)
                case(1)
                   fpvar(i,j,k,:) = ( fp(i,j,k,:)/2 &
                     +fp(i-1,j,k,:)/4+fp(i+1,j,k,:)/4 )
                case(2)
                   fpvar(i,j,k,:) = ( 5*fp(i,j,k,:)/8 &
                     -fp(i-2,j,k,:)/16+fp(i-1,j,k,:)/4 &
                     -fp(i+2,j,k,:)/16+fp(i+1,j,k,:)/4 )
                case(4)
                   fpvar(i,j,k,:) = ( 93*fp(i,j,k,:)/128 &
                     + 7*fp(i-1,j,k,:)/32 - 7*fp(i-2,j,k,:)/64 + fp(i-3,j,k,:)/32 - fp(i-4,j,k,:)/256 &
                     + 7*fp(i+1,j,k,:)/32 - 7*fp(i+2,j,k,:)/64 + fp(i+3,j,k,:)/32 - fp(i+4,j,k,:)/256 )
                case default
                  call bl_error("Input filtering_width must be 0,1,2 or 4")
                end select      
             end do
          end do
       end do
       fp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = fpvar(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:)

       deallocate(fpvar)

    end do FilterX
   
    !--------------------------------------
    ! Update ghost values to reflect filtering
    call multifab_fill_boundary(mfab)
   
    FilterY: do box = 1, nfabs(mfab)

       fp => dataptr(mfab,box) ! Including ghosts
       
       ! Without ghosts:
       lo(1:dm) = lwb(get_ibox(mfab,box))
       hi(1:dm) = upb(get_ibox(mfab,box))
       ! With ghosts:
       lo_g(1:dm) = lwb(get_pbox(mfab,box))
       hi_g(1:dm) = upb(get_pbox(mfab,box))       

       allocate(fpvar(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3), &
                      lbound(fp,4):ubound(fp,4))) ! Temporary array

       ! Filter along y:
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                select case(filtering_width)
                case(1)
                   fpvar(i,j,k,:) = ( fp(i,j,k,:)/2 &
                     +fp(i,j-1,k,:)/4+fp(i,j+1,k,:)/4 )
                case(2)
                   fpvar(i,j,k,:) = ( 5*fp(i,j,k,:)/8 &
                    -fp(i,j-2,k,:)/16+fp(i,j-1,k,:)/4 &
                    -fp(i,j+2,k,:)/16+fp(i,j+1,k,:)/4 )
                case(4)
                   fpvar(i,j,k,:) = ( 93*fp(i,j,k,:)/128 &
                     + 7*fp(i,j-1,k,:)/32 - 7*fp(i,j-2,k,:)/64 + fp(i,j-3,k,:)/32 - fp(i,j-4,k,:)/256 &
                     + 7*fp(i,j+1,k,:)/32 - 7*fp(i,j+2,k,:)/64 + fp(i,j+3,k,:)/32 - fp(i,j+4,k,:)/256 )
                case default
                  call bl_error("Input filtering_width must be 0,1,2 or 4")
                end select      
             end do
          end do
       end do
       fp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = fpvar(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:)

       deallocate(fpvar)

    end do FilterY
    !--------------------------------------

    if(dm<=2) return ! We are done
    
    !--------------------------------------
    call multifab_fill_boundary(mfab) ! Update ghost values to reflect filtering
   
    FilterZ: do box = 1, nfabs(mfab)

       fp => dataptr(mfab,box) ! Including ghosts
       
       ! Without ghosts:
       lo(1:dm) = lwb(get_ibox(mfab,box))
       hi(1:dm) = upb(get_ibox(mfab,box))
       ! With ghosts:
       lo_g(1:dm) = lwb(get_pbox(mfab,box))
       hi_g(1:dm) = upb(get_pbox(mfab,box))       

       allocate(fpvar(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3), &
                      lbound(fp,4):ubound(fp,4))) ! Temporary array

       ! Filter along z:
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                select case(filtering_width)
                case(1)
                   fpvar(i,j,k,:) = ( fp(i,j,k,:)/2 &
                     +fp(i,j,k-1,:)/4+fp(i,j,k+1,:)/4 )
                case(2)
                   fpvar(i,j,k,:) = ( 5*fp(i,j,k,:)/8 &
                    -fp(i,j,k-2,:)/16+fp(i,j,k-1,:)/4 &
                    -fp(i,j,k+2,:)/16+fp(i,j,k+1,:)/4 )
                case(4)
                   fpvar(i,j,k,:) = ( 93*fp(i,j,k,:)/128 &
                     + 7*fp(i,j,k-1,:)/32 - 7*fp(i,j,k-2,:)/64 + fp(i,j,k-3,:)/32 - fp(i,j,k-4,:)/256 &
                     + 7*fp(i,j,k+1,:)/32 - 7*fp(i,j,k+2,:)/64 + fp(i,j,k+3,:)/32 - fp(i,j,k+4,:)/256 )
                case default
                  call bl_error("Input filtering_width must be 0,1,2 or 4")
                end select      
             end do
          end do
       end do
       fp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = fpvar(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:)

       deallocate(fpvar)

    end do FilterZ
    !--------------------------------------

  end subroutine multifab_filter

  ! call this once at the beginning of simulation to allocate multifabs
  ! that will hold random numbers
  subroutine init_stochastic(mla)

    type(ml_layout), intent(in   ) :: mla

    ! local
    integer :: n,nlevs,i,dm
    logical :: nodal_temp(mla%dim)
    
    nlevs = mla%nlevel
    dm = mla%dim

    allocate(sflux_fc(mla%nlevel,mla%dim))
    allocate(mflux_cc(mla%nlevel))
    allocate(mflux_nd(mla%nlevel))
    allocate(mflux_ed(mla%nlevel,3))

    do n=1,nlevs
       if(nscal>1) then
          do i=1,dm
             ! we need one face-centered flux for each concentration
            call multifab_build_edge(sflux_fc(n,i),mla%la(n),nscal-1,filtering_width,i)
         end do
       end if
       ! we need dm cell-centered fluxes for momentum
       call multifab_build(mflux_cc(n),mla%la(n),dm,max(1,filtering_width))
       if (dm .eq. 2) then
          ! in 2D, we need 2 random fluxes at each node
          nodal_temp = .true.
          call multifab_build(mflux_nd(n),mla%la(n),2,filtering_width,nodal_temp)
       else if (dm .eq. 3) then
          ! in 3D, we need 2 random fluxes at each edge
          nodal_temp(1) = .true.
          nodal_temp(2) = .true.
          nodal_temp(3) = .false.
          call multifab_build(mflux_ed(n,1),mla%la(n),2,filtering_width,nodal_temp)
          nodal_temp(1) = .true.
          nodal_temp(2) = .false.
          nodal_temp(3) = .true.
          call multifab_build(mflux_ed(n,2),mla%la(n),2,filtering_width,nodal_temp)
          nodal_temp(1) = .false.
          nodal_temp(2) = .true.
          nodal_temp(3) = .true.
          call multifab_build(mflux_ed(n,3),mla%la(n),2,filtering_width,nodal_temp)
       end if
    end do

  end subroutine init_stochastic

  ! call this once at the end of simulation to deallocate memory
  subroutine destroy_stochastic(mla)

    type(ml_layout), intent(in   ) :: mla

    ! local
    integer :: n,nlevs,i,dm
    
    nlevs = mla%nlevel
    dm = mla%dim

    do n=1,nlevs
       if(nscal>1) then
          do i=1,dm
             call multifab_destroy(sflux_fc(n,i))
         end do
       end if
       call multifab_destroy(mflux_cc(n))
       if (dm .eq. 2) then
          call multifab_destroy(mflux_nd(n))
       else if (dm .eq. 3) then
          call multifab_destroy(mflux_ed(n,1))
          call multifab_destroy(mflux_ed(n,2))
          call multifab_destroy(mflux_ed(n,3))
       end if
    end do
    
    deallocate(mflux_cc,mflux_nd,mflux_ed,sflux_fc)

  end subroutine destroy_stochastic

 ! Add equilibrium fluctuations to the momentum (valid and ghost regions)
 subroutine add_momentum_fluctuations(mla,dx,variance,s_cc,s_face,m_face,mactemp)

   type(ml_layout), intent(in   ) :: mla
   real(dp_t)     , intent(in   ) :: variance, dx(:,:)
   type(multifab) , intent(in   ) :: s_cc(:), s_face(:,:)
   type(multifab) , intent(inout) :: m_face(:,:)  
   type(multifab) , intent(inout) :: mactemp(:,:) ! Temporary multifab

   ! local
   integer :: n,i,dm,nlevs
   real(dp_t) :: av_mom(mla%dim)

   dm = mla%dim
   nlevs = mla%nlevel

   ! Generate random numbers first and store them in umac temporarily  
   do n=1,nlevs
      do i=1,dm
         call multifab_fill_random(mactemp(n:n,i), &
              variance=abs(variance)*kT/product(dx(n,1:dm)), variance_mfab=s_face(n:n,i))
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
      
      call sum_mass_momentum(mla, cons=s_cc, m=m_face, av_momentum=av_mom)
      do i=1,dm
         call setval(mactemp(1,i), -av_mom(i))
         call saxpy(m_face(1,i), 1.0_dp_t, mactemp(1,i))
      end do
   end if            

 end subroutine add_momentum_fluctuations

  subroutine add_concentration_fluctuations(mla, dx, variance, prim, cctemp)
    type(ml_layout), intent(in   ) :: mla
    real(dp_t)     , intent(in   ) :: variance, dx(:,:)
    type(multifab) , intent(inout) :: prim(:)
    type(multifab) , intent(inout) :: cctemp(:)

    ! local
    integer :: i,n,dm,box,nlevs
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
            (fp(:,:,:,2)*mol_mass(1)+(1.0d0-fp(:,:,:,2))*mol_mass(2)) ) * fpvar(:,:,:,1)
            
          ! Now add the fluctuations to the mean:  
          fp(:,:,:,2) = fp(:,:,:,2) + sqrt(abs(variance)/product(dx(n,1:dm))) * fpvar(:,:,:,1)
          ! And calculate the density from the EOS:
          fp(:,:,:,1) = 1.0d0/(fp(:,:,:,2)/rhobar(1)+(1.0d0-fp(:,:,:,2))/rhobar(2))
          
       end do   

    end do    
    
  end subroutine add_concentration_fluctuations

end module mk_stochastic_fluxdiv_module
