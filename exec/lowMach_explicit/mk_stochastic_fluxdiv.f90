module mk_stochastic_fluxdiv_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use ml_layout_module
  use ml_restriction_module
  use define_bc_module
  use bc_module
  
  use BoxLibRNGs
  use analysis_module, only : sum_mass_momentum
  use probin_module, only : variance_coeff, conc_scal, visc_coef, diff_coef, &
      stoch_stress_form, nodal, nscal, kT, mol_mass, rhobar, filtering_width

  implicit none

  private

  public :: mk_stochastic_fluxdiv, create_random_increments, destroy_random_increments, &
            add_momentum_fluctuations, add_concentration_fluctuations, multifab_fill_random

  ! A. Donev: Stochastic fluxes for momentum are generated on:
  ! -cell-centered grid for diagonal components
  ! -node-centered (2D) or edge-centered (3D) grid for off-diagonal components
  type(multifab), allocatable, save :: mflux_cc(:,:), mflux_nd(:,:)
  type(multifab), allocatable, save :: mflux_xy(:,:), mflux_xz(:,:), mflux_yz(:,:)
  ! Stochastic fluxes for scalars are face-centered
  type(multifab), allocatable, save :: sflux(:,:,:)
  
  integer, save :: ntracers=0 ! How many randomly-forced scalars are there
  logical, save :: warn_bad_c=.false. ! Should we issue warnings about c<0 or c>1
  
contains

  ! This is an interface to mk_stochastic_fluxdiv_work
  ! Note that here we *increment* stoch_m_force and stoch_s_force so they must be initialized externally!
  subroutine mk_stochastic_fluxdiv(mla,the_bc_level,stoch_m_force,stoch_s_force,s_face, &
                                   eta,eta_nodal,eta_edge,chi_face,umac,dx,dt,weights)

    type(ml_layout), intent(in   ) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(multifab) , intent(inout) :: stoch_m_force(:,:)
    type(multifab) , intent(inout) :: stoch_s_force(:)
    type(multifab) , intent(in   ) :: s_face(:,:)
    type(multifab) , intent(in   ) :: eta(:)
    type(multifab) , intent(in   ) :: eta_nodal(:)
    type(multifab) , intent(in   ) :: eta_edge(:,:)
    type(multifab) , intent(in   ) :: chi_face(:,:)
    type(multifab) , intent(inout) :: umac(:,:)
    real(dp_t)     , intent(in   ) :: dx(:,:)
    real(dp_t)     , intent(in   ) :: dt
    real(dp_t), intent(in), optional :: weights(:) ! If present, reuse previously-generated rngs

    integer :: i,n,dm,nlevs,box,idim,rng
    logical :: reuse

    nlevs = mla%nlevel
    dm    = mla%dim
    
    reuse=.false.
    
    if(present(weights)) then ! Make a weighted sum of previously-generated random numbers
    if(size(weights)>0) then
       !write(*,*) "REUSING Weiner increments:", weights
       reuse=.true.
       
       call multifab_weighted_sum(mflux_cc, weights)
       do idim = 1,dm
         if(ntracers>0) then
            call multifab_weighted_sum(sflux(:,idim,:), weights)
         end if   
       end do

       if (dm .eq. 2) then
          do n=1,nlevs
             call multifab_weighted_sum(mflux_nd, weights)
          end do
       else if (dm .eq. 3) then
          do n=1,nlevs
             call multifab_weighted_sum(mflux_xy, weights)
             call multifab_weighted_sum(mflux_xz, weights)
             call multifab_weighted_sum(mflux_yz, weights)
          end do       
       end if
       
    end if          
    end if

    call mk_stochastic_fluxdiv_work(mla,the_bc_level,stoch_m_force,stoch_s_force,s_face,dx,dt, &
                             mflux_cc(:,0),mflux_nd(:,0), &
                             mflux_xy(:,0),mflux_xz(:,0),mflux_yz(:,0), &
                             sflux(:,:,0),eta,eta_nodal,eta_edge,chi_face,umac,reuse=reuse)
  
  end subroutine mk_stochastic_fluxdiv

  subroutine mk_stochastic_fluxdiv_work(mla,the_bc_level,stoch_m_force,stoch_s_force,s_face, &
                                 dx,dt,mflux_cc,mflux_nd,mflux_xy,mflux_xz,mflux_yz, &
                                 sflux,eta,eta_nodal,eta_edge,chi_face,umac,reuse)
    
    type(ml_layout), intent(in   ) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(multifab) , intent(inout) :: stoch_m_force(:,:)
    type(multifab) , intent(inout) :: stoch_s_force(:)
    type(multifab) , intent(in   ) :: s_face(:,:)
    real(dp_t)     , intent(in   ) :: dx(:,:)
    real(dp_t)     , intent(in   ) :: dt
    ! The random numbers are generated outside of this routine and passed in
    type(multifab) :: mflux_cc(mla%nlevel), mflux_nd(mla%nlevel)
    type(multifab) :: mflux_xy(mla%nlevel), mflux_xz(mla%nlevel), mflux_yz(mla%nlevel)
    ! Stochastic fluxes for scalars are face-centered
    type(multifab) :: sflux(mla%nlevel,mla%dim)
    type(multifab) , intent(in   ) :: eta(:)
    type(multifab) , intent(in   ) :: eta_nodal(:)
    type(multifab) , intent(in   ) :: eta_edge(:,:)
    type(multifab) , intent(in   ) :: chi_face(:,:)
    type(multifab) , intent(inout) :: umac(:,:)
    logical, intent(in) :: reuse ! Reuse old rngs or generate them on the fly

    ! local variables
    ! --------------------------
    integer :: i,n,dm,nlevs,box,idim,ng_c,ng_n,ng_e,ng_f,ng_s,ng_x,ng_y,ng_z,ng_m,ng_w
    real(kind=dp_t), pointer :: fp(:,:,:,:), sp(:,:,:,:), dp(:,:,:,:)
    real(kind=dp_t), pointer :: fxp(:,:,:,:), fyp(:,:,:,:), fzp(:,:,:,:)
    real(kind=dp_t), pointer :: dxp(:,:,:,:), dyp(:,:,:,:), dzp(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:), vmp(:,:,:,:), wmp(:,:,:,:)
    real(kind=dp_t), pointer :: ep1(:,:,:,:), ep2(:,:,:,:), ep3(:,:,:,:)

    real(dp_t) :: variance

    integer :: lo(mla%dim), hi(mla%dim)
  
    !--------------------------------------
    ! Begin execution:
    nlevs = mla%nlevel
    dm    = mla%dim

    ng_f  = stoch_m_force(1,1)%ng
    ng_s  = s_face(1,1)%ng
    
    ng_c = mflux_cc(1)%ng
    ng_n = mflux_nd(1)%ng
    ng_e = mflux_xy(1)%ng

    ng_x = sflux(1,1)%ng
    ng_y = eta(1)%ng
    ng_z = chi_face(1,1)%ng

    ng_m = umac(1,1)%ng
        
    !--------------------------------------    
    ! Fill the stochastic flux multifabs with random numbers:
    
    if (dm .eq. 2) then

       ng_w = eta_nodal(1)%ng

       do n=1,nlevs           
          
          if (visc_coef >= 0) then
             variance = sqrt(variance_coeff*2.d0*kT*visc_coef/(product(dx(n,1:dm))*dt))
          else
             ! exclude the visc_coef contribution since we first need to average to nodes - 
             ! --- to be added below
             variance = sqrt(variance_coeff*2.d0*kT/(product(dx(n,1:dm))*dt))
          end if
          
          do box = 1, nfabs(mflux_cc(n))
             
             fp  => dataptr(mflux_cc(n),box)
             sp  => dataptr(mflux_nd(n),box)
             dp  => dataptr(eta(n),box)
             ep1 => dataptr(eta_nodal(n),box)
             lo = lwb(get_box(mflux_cc(n),box))
             hi = upb(get_box(mflux_cc(n),box))
             
             if(.not.reuse) then ! Fill the whole grid with random numbers
               select case(stoch_stress_form)
               case(0) ! Non-symmetric
                call NormalRNGs(fp, size(fp)) 
                call NormalRNGs(sp, size(sp))
               case default ! Symmetric
                call NormalRNGs(fp, size(fp))
                fp=sqrt(2.0d0)*fp
                call NormalRNGs(sp(:,:,:,1), size(sp(:,:,:,1)))
                sp(:,:,:,2)=sp(:,:,:,1)
               end select              
             end if
             
             fp = variance*fp
             sp = variance*sp
             if (visc_coef < 0) then
                call mult_by_sqrt_eta_2d(fp(:,:,1,:),ng_c,sp(:,:,1,:),ng_n, &
                                         dp(:,:,1,1),ng_y,ep1(:,:,1,1),ng_w,lo,hi)
             end if
             
             call mflux_bc_2d(sp(:,:,1,:),ng_n,lo,hi, &
                              the_bc_level(n)%phys_bc_level_array(box,:,:))
          end do
          
          ! Now sync up the random numbers at the boundaries:
          call multifab_internal_sync(mflux_nd(n))
          call multifab_fill_boundary(mflux_nd(n))
          call multifab_fill_boundary(mflux_cc(n))
          
          if(filtering_width>0) then
             call multifab_filter(mflux_nd(n), dm)
             call multifab_filter(mflux_cc(n), dm)
             call multifab_fill_boundary(mflux_cc(n)) ! First ghost cell is used in divergence
          end if
          
       enddo

    else if (dm .eq. 3) then

       ng_w = eta_edge(1,1)%ng

       do n=1,nlevs              
          
          if (visc_coef >= 0) then
             variance = sqrt(variance_coeff*2.d0*kT*visc_coef/(product(dx(n,1:dm))*dt))
          else
             ! exclude the visc_coef contribution since we first need to average to nodes - 
             ! --- to be added below
             variance = sqrt(variance_coeff*2.d0*kT/(product(dx(n,1:dm))*dt))
          end if
          
          do box = 1, nfabs(mflux_cc(n))
             
             fp  => dataptr(mflux_cc(n),box)
             fxp => dataptr(mflux_xy(n),box)
             fyp => dataptr(mflux_xz(n),box)
             fzp => dataptr(mflux_yz(n),box)
             dp => dataptr(eta(n),box)
             ep1 => dataptr(eta_edge(n,1),box)
             ep2 => dataptr(eta_edge(n,2),box)
             ep3 => dataptr(eta_edge(n,3),box)
             lo = lwb(get_box(mflux_cc(n),box))
             hi = upb(get_box(mflux_cc(n),box))
             
             if(.not.reuse) then
               select case(stoch_stress_form)
               case(0) ! Non-symmetric
                call NormalRNGs(fp, size(fp)) ! Fill the whole grid with random numbers
                call NormalRNGs(fxp, size(fxp))
                call NormalRNGs(fyp, size(fyp))
                call NormalRNGs(fzp, size(fzp))
               case default ! Symmetric
                call NormalRNGs(fp, size(fp)) ! Fill the whole grid with random numbers
                fp=sqrt(2.0d0)*fp
                call NormalRNGs(fxp(:,:,:,1), size(fxp(:,:,:,1)))
                fxp(:,:,:,2)=fxp(:,:,:,1)
                call NormalRNGs(fyp(:,:,:,1), size(fyp(:,:,:,1)))
                fyp(:,:,:,2)=fyp(:,:,:,1)
                call NormalRNGs(fzp(:,:,:,1), size(fzp(:,:,:,1)))
                fzp(:,:,:,2)=fzp(:,:,:,1)                
               end select
             end if

             fp  = variance*fp
             fxp = variance*fxp
             fyp = variance*fyp
             fzp = variance*fzp
             if (visc_coef < 0) then
                call mult_by_sqrt_eta_3d(fp(:,:,:,:),ng_c, &
                                         fxp(:,:,:,:),fyp(:,:,:,:),fzp(:,:,:,:),ng_e, &
                                         dp(:,:,:,1),ng_y, &
                                         ep1(:,:,:,1),ep2(:,:,:,1),ep3(:,:,:,1),ng_w,lo,hi)
             end if
             
             call mflux_bc_3d(fxp(:,:,:,:),fyp(:,:,:,:),fzp(:,:,:,:),ng_e,lo,hi, &
                              the_bc_level(n)%phys_bc_level_array(box,:,:))
          end do
          
          ! Now sync up the random numbers at the boundaries:
          call multifab_internal_sync(mflux_xy(n))
          call multifab_internal_sync(mflux_xz(n))
          call multifab_internal_sync(mflux_yz(n))
          call multifab_fill_boundary(mflux_xy(n))
          call multifab_fill_boundary(mflux_xz(n))
          call multifab_fill_boundary(mflux_yz(n))
          call multifab_fill_boundary(mflux_cc(n))

          if(filtering_width>0) then
             call multifab_filter(mflux_xy(n), dm)
             call multifab_filter(mflux_xz(n), dm)
             call multifab_filter(mflux_yz(n), dm)
             call multifab_filter(mflux_cc(n), dm)
             call multifab_fill_boundary(mflux_cc(n)) ! First ghost cell is used in divergence
          end if
          
       enddo

    end if

    !--------------------------------------
    ! Now calculate the MAC divergence to get the actual forcing:
    do n=1,nlevs

       do i=1,nfabs(stoch_m_force(n,1))
          fp => dataptr(mflux_cc(n), i)
          dxp => dataptr(stoch_m_force(n,1),i)
          dyp => dataptr(stoch_m_force(n,2),i)
          lo =  lwb(get_box(stoch_m_force(n,1), i))
          hi =  upb(get_box(stoch_m_force(n,1), i))
          select case (dm)
          case (2)
             sp => dataptr(mflux_nd(n), i)
             call stoch_m_force_2d(fp(:,:,1,:), sp(:,:,1,:), dxp(:,:,1,1), dyp(:,:,1,1), &
                                   ng_c, ng_n, ng_f, dx(n,:), lo, hi)
          case (3)
             dzp => dataptr(stoch_m_force(n,3), i)
             fxp => dataptr(mflux_xy(n), i)
             fyp => dataptr(mflux_xz(n), i)
             fzp => dataptr(mflux_yz(n), i)
             call stoch_m_force_3d(fp(:,:,:,:), fxp(:,:,:,:), fyp(:,:,:,:), fzp(:,:,:,:), &
                                   dxp(:,:,:,1), dyp(:,:,:,1), dzp(:,:,:,1), &
                                   ng_c, ng_e, ng_f, dx(n,:), lo, hi)

          end select
       end do

       if(.false.) then
          ! This might in principle be needed in case of roundoff differences among processors:
          do idim=1,dm
             call multifab_internal_sync(stoch_m_force(n,idim))
             call multifab_fill_boundary(stoch_m_force(n,idim))
          end do
       end if

    end do
         
    if(ntracers>0) then

       do n=1,nlevs       
       
          if (diff_coef >= 0) then
             variance = sqrt(variance_coeff*conc_scal*2.d0*diff_coef/(product(dx(n,1:dm))*dt))
          else
             ! exclude the diff_coef contribution since we first need to average to faces - 
             ! --- to be added below
             variance = sqrt(variance_coeff*conc_scal*2.d0/(product(dx(n,1:dm))*dt))
          end if
             
          do idim = 1, dm       
             do box = 1, nfabs(sflux(n,idim))
                fp  => dataptr(sflux(n,idim),box)
                sp => dataptr(s_face(n,idim),box)
                dp => dataptr(chi_face(n,idim),box)
                lo =  lwb(get_box(sflux(n,idim),box))
                hi =  upb(get_box(sflux(n,idim),box))
                if(.not.reuse) then
                   call NormalRNGs(fp, size(fp)) ! Fill the whole grid with random numbers
                end if
                fp = variance*fp
                ! Include multiplicative rho*c*(1-c) scaling for random rho*c flux
                ! also sets rho*c flux on walls to zero
                select case (dm)
                case (2)
                   if (diff_coef < 0) then
                      call mult_by_sqrt_chi_2d(fp(:,:,1,:),ng_x,dp(:,:,1,1),ng_z,idim, &
                                               lo,hi,ntracers)
                   end if
                  call scale_rhoc_2d(fp(:,:,1,:),ng_x,sp(:,:,1,1:),ng_s,idim,lo,hi, &
                                     the_bc_level(n)%phys_bc_level_array(box,:,:),ntracers)
                case (3)
                   if (diff_coef < 0) then
                      call mult_by_sqrt_chi_3d(fp(:,:,:,:),ng_x,dp(:,:,:,1),ng_z,idim, &
                                               lo,hi,ntracers)

                   end if
                  call scale_rhoc_3d(fp(:,:,:,:),ng_x,sp(:,:,:,1:),ng_s,idim,lo,hi, &
                                     the_bc_level(n)%phys_bc_level_array(box,:,:),ntracers)
                end select
             end do
             
             ! Now sync up the random numbers at the boundaries:
             call multifab_internal_sync(sflux(n,idim))
             call multifab_fill_boundary(sflux(n,idim))          

             if(filtering_width>0) then
                call multifab_filter(sflux(n,idim), dm)
             end if
             
          end do
          
       enddo

       !--------------------------------------
       ! Now calculate the MAC divergence to get the actual forcing:
       do n = 1,nlevs
       
          do i = 1, nfabs(stoch_s_force(n))
             fp  => dataptr(stoch_s_force(n),i)
             fxp => dataptr(sflux(n,1), i)
             fyp => dataptr(sflux(n,2), i)
             ump => dataptr(umac(n,1), i)
             vmp => dataptr(umac(n,2), i)
             lo =  lwb(get_box(stoch_s_force(n), i))
             hi =  upb(get_box(stoch_s_force(n), i))
             select case (dm)
             case (2)
                fp(:,:,1,1)=0.0_dp_t ! No stochastic forcing for density!
                call stoch_s_force_2d(fxp(:,:,1,:), fyp(:,:,1,:), ng_x, &
                                      fp(:,:,1,2:), ng_f, &
                                      ump(:,:,1,1), vmp(:,:,1,1), ng_m, &
                                      dx(n,:),lo,hi, &
                                      the_bc_level(n)%adv_bc_level_array(i,:,:,:))
             case (3)
                fzp => dataptr(sflux(n,3), i)
                wmp => dataptr(umac(n,3), i)
                fp(:,:,:,1)=0.0_dp_t ! No stochastic forcing for density!
                call stoch_s_force_3d(fxp(:,:,:,:), fyp(:,:,:,:), fzp(:,:,:,:), ng_x, &
                                      fp(:,:,:,2:), ng_f, &
                                      ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_m, &
                                      dx(n,:),lo,hi, &
                                      the_bc_level(n)%adv_bc_level_array(i,:,:,:))
             end select
          end do

          if(.false.) then
             ! This might in principle be needed in case of roundoff differences among processors:
             call multifab_fill_boundary(stoch_s_force(n))
          end if
             
       end do

    end if

  end subroutine mk_stochastic_fluxdiv_work

  ! Create the stochastic flux multifabs, possibly allowing for storing random numbers
  subroutine create_random_increments(mla,n_rngs)
   type(ml_layout), intent(in   ) :: mla
   integer, intent(in) :: n_rngs ! How many random numbers to store per time step
      ! Could be zero if one does not need to store any rngs
      
   ! Local variables
   integer :: i,n,dm,nlevs,box,idim,rng
   logical :: nodal_temp(mla%dim)
   
   ! The zeroth component is used to store the actual stochastic flux
   ! The rest are used to store random numbers that may be reused later
   allocate(mflux_cc(mla%nlevel,0:n_rngs))
   allocate(mflux_nd(mla%nlevel,0:n_rngs))
   allocate(mflux_xy(mla%nlevel,0:n_rngs))
   allocate(mflux_xz(mla%nlevel,0:n_rngs))
   allocate(mflux_yz(mla%nlevel,0:n_rngs))
   allocate(sflux(mla%nlevel,mla%dim,0:n_rngs))

    ! Create temporary multifabs
    !--------------------------------------
    ! Begin execution:
    nlevs = mla%nlevel
    dm    = mla%dim    
    ntracers = nscal - 1 ! AD FIXME: Not sure what we want here exactly...
    
    do rng=0, n_rngs
    do n = 1, nlevs
       call multifab_build(mflux_cc(n,rng),mla%la(n),dm,max(1,filtering_width))
       if(ntracers>0) then
          do idim = 1,dm
            call multifab_build_edge(sflux(n,idim,rng),mla%la(n),ntracers,filtering_width,idim)
         end do
       end if

       if (dm .eq. 2) then
          ! in 2D, we need 2 random fluxes at each node
          call multifab_build(mflux_nd(n,rng),mla%la(n),2,filtering_width,nodal)
       else if (dm .eq. 3) then
          ! in 3D, we need 2 random fluxes at each edge
          nodal_temp(1) = .true.
          nodal_temp(2) = .true.
          nodal_temp(3) = .false.
          call multifab_build(mflux_xy(n,rng),mla%la(n),2,filtering_width,nodal_temp)
          nodal_temp(1) = .true.
          nodal_temp(2) = .false.
          nodal_temp(3) = .true.
          call multifab_build(mflux_xz(n,rng),mla%la(n),2,filtering_width,nodal_temp)
          nodal_temp(1) = .false.
          nodal_temp(2) = .true.
          nodal_temp(3) = .true.
          call multifab_build(mflux_yz(n,rng),mla%la(n),2,filtering_width,nodal_temp)
       end if
   end do
   end do
   
   ! Generate and store the actual random numbers here
   do rng=1, n_rngs
      
      ! Diagonal components of stochastic stress tensor:
      select case(stoch_stress_form)
      case(0) ! Non-symmetric
         call multifab_fill_random(mflux_cc(:,rng)) ! Variance=1
      case default ! Symmetric
         call multifab_fill_random(mflux_cc(:,rng), variance=2.0d0)
      end select   

      ! Off-diagonal components of stochastic stress tensor:
      if (dm .eq. 2) then
         ! in 2D, we need 2 random fluxes at each node
         select case(stoch_stress_form)
         case(0) ! Non-symmetric
            call multifab_fill_random(mflux_nd(:,rng))
         case default ! Symmetric
            call multifab_fill_random(mflux_nd(:,rng), comp=1)
            do n = 1, nlevs
               call multifab_copy_c(mflux_nd(n,rng),2, mflux_nd(n,rng),1)
            end do   
         end select
      else if (dm .eq. 3) then
         ! in 3D, we need 2 random fluxes at each edge
         select case(stoch_stress_form)
         case(0) ! Non-symmetric
            call multifab_fill_random(mflux_xy(:,rng))
            call multifab_fill_random(mflux_xz(:,rng))
            call multifab_fill_random(mflux_yz(:,rng))
         case default ! Symmetric
            call multifab_fill_random(mflux_xy(:,rng), comp=1)
            call multifab_fill_random(mflux_xz(:,rng), comp=1)
            call multifab_fill_random(mflux_yz(:,rng), comp=1)
            do n = 1, nlevs
               call multifab_copy_c(mflux_xy(n,rng),2, mflux_xy(n,rng),1)
               call multifab_copy_c(mflux_xz(n,rng),2, mflux_xz(n,rng),1)
               call multifab_copy_c(mflux_yz(n,rng),2, mflux_yz(n,rng),1)
            end do   
         end select             
      end if
            
      if(ntracers>0) then ! Stochastic diffusive flux
         do idim = 1,dm
            call multifab_fill_random(sflux(:,idim,rng))
         end do   
      end if
      
   end do   
  
  end subroutine
  
  subroutine destroy_random_increments(mla, n_rngs)
    type(ml_layout), intent(in   ) :: mla
   integer, intent(in) :: n_rngs

    integer :: i,n,dm,nlevs,box,idim,rng

    nlevs = mla%nlevel
    dm    = mla%dim    
  
    !--------------------------------------
    ! Destroy the multifabs:
    do rng=0, n_rngs 
    do n = 1, nlevs
       call multifab_destroy(mflux_cc(n,rng))
       do idim = 1,dm
         if(ntracers>0) then
            call multifab_destroy(sflux(n,idim,rng))
         end if   
       end do

       if (dm .eq. 2) then
          call multifab_destroy(mflux_nd(n,rng))
       else if (dm .eq. 3) then
          call multifab_destroy(mflux_xy(n,rng))
          call multifab_destroy(mflux_xz(n,rng))
          call multifab_destroy(mflux_yz(n,rng))
       end if
    end do
    end do

    deallocate(mflux_cc,mflux_nd,mflux_xy,mflux_xz,mflux_yz,sflux)
    
  end subroutine
  
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

  subroutine mult_by_sqrt_chi_2d(sflux,ng_x,chi_face,ng_z,idim,lo,hi,ntracers)

    integer        , intent(in   ) :: lo(:),hi(:),ng_x,ng_z,idim,ntracers
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

  subroutine mult_by_sqrt_chi_3d(sflux,ng_x,chi_face,ng_z,idim,lo,hi,ntracers)

    integer        , intent(in   ) :: lo(:),hi(:),ng_x,ng_z,idim,ntracers
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
  
  subroutine scale_rhoc_2d(sflux,ng_x,s_face,ng_s,idim,lo,hi,phys_bc,ntracers)

    integer        , intent(in   ) :: lo(:),hi(:),ng_x,ng_s,idim,ntracers
    integer        , intent(in   ) :: phys_bc(:,:)
    real(kind=dp_t), intent(inout) ::  sflux(lo(1)-ng_x:,lo(2)-ng_x:,:)
    real(kind=dp_t), intent(in   ) :: s_face(lo(1)-ng_s:,lo(2)-ng_s:,:)

    integer :: i,j
    real(kind=dp_t) :: rho, fac(ntracers), s_fac(ntracers)

    if (idim .eq. 1) then

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1
             rho = s_face(i,j,1) ! Density at face
             fac(:) = s_face(i,j,2:)/s_face(i,j,1) ! Concentration at face
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

       if (phys_bc(1,1) .eq. INLET) then
          sflux(lo(1),lo(2):hi(2),:) = sqrt(2.0d0)*sflux(lo(1),lo(2):hi(2),:)
       end if

       if (phys_bc(1,2) .eq. INLET) then
          sflux(hi(1)+1,lo(2):hi(2),:) = sqrt(2.0d0)*sflux(hi(1)+1,lo(2):hi(2),:)
       end if

    else

       do j = lo(2),hi(2)+1
          do i = lo(1),hi(1)
             rho = s_face(i,j,1) ! Density at face
             fac(:) = s_face(i,j,2:)/s_face(i,j,1) ! Concentration at face
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

       if (phys_bc(2,1) .eq. INLET) then
          sflux(lo(1):hi(1),lo(2),:) = sqrt(2.0d0)*sflux(lo(1):hi(1),lo(2),:)
       end if
       
       if (phys_bc(2,2) .eq. INLET) then
          sflux(lo(1):hi(1),hi(2)+1,:) = sqrt(2.0d0)*sflux(lo(1):hi(1),hi(2)+1,:)
       end if

    end if

  end subroutine scale_rhoc_2d

  subroutine scale_rhoc_3d(sflux,ng_x,s_face,ng_s,idim,lo,hi,phys_bc,ntracers)

    integer        , intent(in   ) :: lo(:),hi(:),ng_x,ng_s,idim,ntracers
    integer        , intent(in   ) :: phys_bc(:,:)
    real(kind=dp_t), intent(inout) ::  sflux(lo(1)-ng_x:,lo(2)-ng_x:,lo(3)-ng_x:,:)
    real(kind=dp_t), intent(in   ) :: s_face(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)

    integer :: i,j,k
    real(kind=dp_t) :: rho, fac(ntracers), s_fac(ntracers)
    
    if (idim .eq. 1) then

       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)+1
                rho = s_face(i,j,k,1) ! Density at face
                fac(:) = s_face(i,j,k,2:)/s_face(i,j,k,1) ! Concentration at face
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

       if (phys_bc(1,1) .eq. INLET) then
          sflux(lo(1),lo(2):hi(2),lo(3):hi(3),:) = sqrt(2.0d0)*sflux(lo(1),lo(2):hi(2),lo(3):hi(3),:)
       end if

       if (phys_bc(1,2) .eq. INLET) then
          sflux(hi(1)+1,lo(2):hi(2),lo(3):hi(3),:) = sqrt(2.0d0)*sflux(hi(1)+1,lo(2):hi(2),lo(3):hi(3),:)
       end if

    else if (idim .eq. 2) then

       do k = lo(3),hi(3)
          do j = lo(2),hi(2)+1
             do i = lo(1),hi(1)
                rho = s_face(i,j,k,1) ! Density at face
                fac(:) = s_face(i,j,k,2:)/s_face(i,j,k,1) ! Concentration at face
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

       if (phys_bc(2,1) .eq. INLET) then
          sflux(lo(1):hi(1),lo(2),lo(3):hi(3),:) = sqrt(2.0d0)*sflux(lo(1):hi(1),lo(2),lo(3):hi(3),:)
       end if
       
       if (phys_bc(2,2) .eq. INLET) then
          sflux(lo(1):hi(1),hi(2)+1,lo(3):hi(3),:) = sqrt(2.0d0)*sflux(lo(1):hi(1),hi(2)+1,lo(3):hi(3),:)
       end if

    else

       do k = lo(3),hi(3)+1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                rho = s_face(i,j,k,1) ! Density at face
                fac(:) = s_face(i,j,k,2:)/s_face(i,j,k,1) ! Concentration at face
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

       if (phys_bc(3,1) .eq. INLET) then
          sflux(lo(1):hi(1),lo(2):hi(2),lo(3),:) = sqrt(2.0d0)*sflux(lo(1):hi(1),lo(2):hi(2),lo(3),:)
       end if
       
       if (phys_bc(3,2) .eq. INLET) then
          sflux(lo(1):hi(1),lo(2):hi(2),hi(3)+1,:) = sqrt(2.0d0)*sflux(lo(1):hi(1),lo(2):hi(2),hi(3)+1,:)
       end if

    end if

  end subroutine scale_rhoc_3d

  subroutine mflux_bc_2d(mflux_nd,ng_n,lo,hi,phys_bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng_n
    real(kind=dp_t), intent(inout) :: mflux_nd(lo(1)-ng_n:,lo(2)-ng_n:,:)
    integer        , intent(in   ) :: phys_bc(:,:)

    ! y-mom fluxes that live on x-domain boundaries
    if (phys_bc(1,1) .eq. NO_SLIP_WALL .or. phys_bc(1,1) .eq. INLET) then
       mflux_nd(lo(1),lo(2):hi(2)+1,:) = sqrt(2.d0)*mflux_nd(lo(1),lo(2):hi(2)+1,:)
    else if (phys_bc(1,1) .eq. SLIP_WALL) then
       mflux_nd(lo(1),lo(2):hi(2)+1,:) = 0.d0
    end if
    if (phys_bc(1,2) .eq. NO_SLIP_WALL .or. phys_bc(1,2) .eq. INLET) then
       mflux_nd(hi(1)+1,lo(2):hi(2)+1,:) = sqrt(2.d0)*mflux_nd(hi(1)+1,lo(2):hi(2)+1,:)
    else if (phys_bc(1,2) .eq. SLIP_WALL) then
       mflux_nd(hi(1)+1,lo(2):hi(2)+1,:) = 0.d0
    end if

    ! x-mom fluxes that live on y-domain boundaries
    if (phys_bc(2,1) .eq. NO_SLIP_WALL .or. phys_bc(2,1) .eq. INLET) then
       mflux_nd(lo(1):hi(1)+1,lo(2),:) = sqrt(2.d0)*mflux_nd(lo(1):hi(1)+1,lo(2),:)
    else if (phys_bc(2,1) .eq. SLIP_WALL) then
       mflux_nd(lo(1):hi(1)+1,lo(2),:) = 0.d0
    end if
    if (phys_bc(2,2) .eq. NO_SLIP_WALL .or. phys_bc(2,2) .eq. INLET) then
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
    if (phys_bc(1,1) .eq. NO_SLIP_WALL .or. phys_bc(1,1) .eq. INLET) then
       mflux_xy(lo(1),lo(2):hi(2)+1,lo(3):hi(3),:) = sqrt(2.d0)*mflux_xy(lo(1),lo(2):hi(2)+1,lo(3):hi(3),:)
       mflux_xz(lo(1),lo(2):hi(2),lo(3):hi(3)+1,:) = sqrt(2.d0)*mflux_xz(lo(1),lo(2):hi(2),lo(3):hi(3)+1,:)
    else if (phys_bc(1,1) .eq. SLIP_WALL) then
       mflux_xy(lo(1),lo(2):hi(2)+1,lo(3):hi(3),:) = 0.d0
       mflux_xz(lo(1),lo(2):hi(2),lo(3):hi(3)+1,:) = 0.d0
    end if
    if (phys_bc(1,2) .eq. NO_SLIP_WALL .or. phys_bc(1,2) .eq. INLET) then
       mflux_xy(hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3),:) = sqrt(2.d0)*mflux_xy(hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3),:)
       mflux_xz(hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1,:) = sqrt(2.d0)*mflux_xz(hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1,:)
    else if (phys_bc(1,2) .eq. SLIP_WALL) then
       mflux_xy(hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3),:) = 0.d0
       mflux_xz(hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1,:) = 0.d0
    end if

    ! x-mom and z-mom fluxes that live on y-domain boundaries
    if (phys_bc(2,1) .eq. NO_SLIP_WALL .or. phys_bc(2,1) .eq. INLET) then
       mflux_xy(lo(1):hi(1)+1,lo(2),lo(3):hi(3),:) = sqrt(2.d0)*mflux_xy(lo(1):hi(1)+1,lo(2),lo(3):hi(3),:)
       mflux_yz(lo(1):hi(1),lo(2),lo(3):hi(3)+1,:) = sqrt(2.d0)*mflux_yz(lo(1):hi(1),lo(2),lo(3):hi(3)+1,:)
    else if (phys_bc(2,1) .eq. SLIP_WALL) then
       mflux_xy(lo(1):hi(1)+1,lo(2),lo(3):hi(3),:) = 0.d0
       mflux_yz(lo(1):hi(1),lo(2),lo(3):hi(3)+1,:) = 0.d0
    end if
    if (phys_bc(2,2) .eq. NO_SLIP_WALL .or. phys_bc(2,2) .eq. INLET) then
       mflux_xy(lo(1):hi(1)+1,hi(2)+1,lo(3):hi(3),:) = sqrt(2.d0)*mflux_xy(lo(1):hi(1)+1,hi(2)+1,lo(3):hi(3),:)
       mflux_yz(lo(1):hi(1),hi(2)+1,lo(3):hi(3)+1,:) = sqrt(2.d0)*mflux_yz(lo(1):hi(1),hi(2)+1,lo(3):hi(3)+1,:)
    else if (phys_bc(2,2) .eq. SLIP_WALL) then
       mflux_xy(lo(1):hi(1)+1,hi(2)+1,lo(3):hi(3),:) = 0.d0
       mflux_yz(lo(1):hi(1),hi(2)+1,lo(3):hi(3)+1,:) = 0.d0
    end if

    ! x-mom and y-mom fluxes that live on z-domain boundaries
    if (phys_bc(3,1) .eq. NO_SLIP_WALL .or. phys_bc(3,1) .eq. INLET) then
       mflux_xz(lo(1):hi(1)+1,lo(2):hi(2),lo(3),:) = sqrt(2.d0)*mflux_xz(lo(1):hi(1)+1,lo(2):hi(2),lo(3),:)
       mflux_yz(lo(1):hi(1),lo(2):hi(2)+1,lo(3),:) = sqrt(2.d0)*mflux_yz(lo(1):hi(1),lo(2):hi(2)+1,lo(3),:)
    else if (phys_bc(3,1) .eq. SLIP_WALL) then
       mflux_xz(lo(1):hi(1)+1,lo(2):hi(2),lo(3),:) = 0.d0
       mflux_yz(lo(1):hi(1),lo(2):hi(2)+1,lo(3),:) = 0.d0
    end if
    if (phys_bc(3,2) .eq. NO_SLIP_WALL .or. phys_bc(3,2) .eq. INLET) then
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

 subroutine stoch_s_force_2d(xflux,yflux,ng_x,div,ng_f,umac,vmac,ng_m,dx,lo,hi,adv_bc)

   use probin_module, only: rhobar

   integer        , intent(in   ) :: lo(:),hi(:),ng_f,ng_m,ng_x
   real(kind=dp_t), intent(in   ) :: xflux(lo(1)-ng_x:,lo(2)-ng_x:,:)
   real(kind=dp_t), intent(in   ) :: yflux(lo(1)-ng_x:,lo(2)-ng_x:,:)
   real(kind=dp_t), intent(inout) ::   div(lo(1)-ng_f:,lo(2)-ng_f:,:)
   real(kind=dp_t), intent(inout) ::  umac(lo(1)-ng_m:,lo(2)-ng_m:)
   real(kind=dp_t), intent(inout) ::  vmac(lo(1)-ng_m:,lo(2)-ng_m:)
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
   
   ! update umac based on diffusive flux at boundary
   if (adv_bc(1,1,1) .eq. EXT_DIR .and. adv_bc(1,1,4) .eq. EXT_DIR) then
      umac(lo(1),lo(2):hi(2)) = umac(lo(1),lo(2):hi(2)) &
           + S_fac*xflux(lo(1),lo(2):hi(2),1)
   end if
   if (adv_bc(1,2,1) .eq. EXT_DIR .and. adv_bc(1,2,4) .eq. EXT_DIR) then
      umac(hi(1)+1,lo(2):hi(2)) = umac(hi(1)+1,lo(2):hi(2)) &
           + S_fac*xflux(hi(1)+1,lo(2):hi(2),1)
   end if

   ! update vmac based on diffusive flux at boundary
   if (adv_bc(2,1,2) .eq. EXT_DIR .and. adv_bc(2,1,4) .eq. EXT_DIR) then
      vmac(lo(1):hi(1),lo(2)) = vmac(lo(1):hi(1),lo(2)) + &
           S_fac*yflux(lo(1):hi(1),lo(2),1)
   end if
   if (adv_bc(2,2,2) .eq. EXT_DIR .and. adv_bc(2,2,4) .eq. EXT_DIR) then
      vmac(lo(1):hi(1),hi(2)+1) = vmac(lo(1):hi(1),hi(2)+1) + &
           S_fac*yflux(lo(1):hi(1),hi(2)+1,1)
   end if

 end subroutine stoch_s_force_2d

 subroutine stoch_s_force_3d(xflux,yflux,zflux,ng_x,div,ng_f,umac,vmac,wmac,ng_m, &
                             dx,lo,hi,adv_bc)

   use probin_module, only: rhobar

   integer        , intent(in   ) :: lo(:),hi(:),ng_f,ng_m,ng_x
   real(kind=dp_t), intent(in   ) :: xflux(lo(1)-ng_x:,lo(2)-ng_x:,lo(3)-ng_x:,:)
   real(kind=dp_t), intent(in   ) :: yflux(lo(1)-ng_x:,lo(2)-ng_x:,lo(3)-ng_x:,:)
   real(kind=dp_t), intent(in   ) :: zflux(lo(1)-ng_x:,lo(2)-ng_x:,lo(3)-ng_x:,:)
   real(kind=dp_t), intent(inout) ::   div(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:,:)
   real(kind=dp_t), intent(inout) ::  umac(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
   real(kind=dp_t), intent(inout) ::  vmac(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
   real(kind=dp_t), intent(inout) ::  wmac(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
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
   
   ! update umac based on diffusive flux at boundary
   if (adv_bc(1,1,1) .eq. EXT_DIR .and. adv_bc(1,1,5) .eq. EXT_DIR) then
      umac(lo(1),lo(2):hi(2),lo(3):hi(3)) = umac(lo(1),lo(2):hi(2),lo(3):hi(3)) &
           + S_fac*xflux(lo(1),lo(2):hi(2),lo(3):hi(3),1)
   end if
   if (adv_bc(1,2,1) .eq. EXT_DIR .and. adv_bc(1,2,5) .eq. EXT_DIR) then
      umac(hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = umac(hi(1)+1,lo(2):hi(2),lo(3):hi(3)) &
           + S_fac*xflux(hi(1)+1,lo(2):hi(2),lo(3):hi(3),1)
   end if

   ! update vmac based on diffusive flux at boundary
   if (adv_bc(2,1,2) .eq. EXT_DIR .and. adv_bc(2,1,5) .eq. EXT_DIR) then
      vmac(lo(1):hi(1),lo(2),lo(3):hi(3)) = vmac(lo(1):hi(1),lo(2),lo(3):hi(3)) + &
           S_fac*yflux(lo(1):hi(1),lo(2),lo(3):hi(3),1)
   end if
   if (adv_bc(2,2,2) .eq. EXT_DIR .and. adv_bc(2,2,5) .eq. EXT_DIR) then
      vmac(lo(1):hi(1),hi(2)+1,lo(3):hi(3)) = vmac(lo(1):hi(1),hi(2)+1,lo(3):hi(3)) + &
           S_fac*yflux(lo(1):hi(1),hi(2)+1,lo(3):hi(3),1)
   end if

   ! update wmac based on diffusive flux at boundary
   if (adv_bc(3,1,3) .eq. EXT_DIR .and. adv_bc(3,1,5) .eq. EXT_DIR) then
      wmac(lo(1):hi(1),lo(2):hi(2),lo(3)) = wmac(lo(1):hi(1),lo(2):hi(2),lo(3)) + &
           S_fac*zflux(lo(1):hi(1),lo(2):hi(2),lo(3),1)
   end if
   if (adv_bc(3,2,3) .eq. EXT_DIR .and. adv_bc(3,2,5) .eq. EXT_DIR) then
      wmac(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = wmac(lo(1):hi(1),lo(2):hi(2),hi(3)+1) + &
           S_fac*zflux(lo(1):hi(1),lo(2):hi(2),hi(3)+1,1)
   end if

 end subroutine stoch_s_force_3d

 !--------------------------------------------------
 ! Equilibrium fluctuations 
 !--------------------------------------------------

  ! The amplitude of the stochastic forcing in the concentration equation is
  subroutine concentration_amplitude(s_fac, fac, rho)
    real(kind=dp_t), intent(in) :: rho, fac(ntracers)
    real(kind=dp_t), intent(out) :: s_fac(ntracers)
    
    ! For ideal mixture use: rho*c*(1-c)*M
    s_fac(:) = rho*fac(:)*(1.0d0-fac(:))* (fac(:)*mol_mass(1)+(1.0d0-fac(:))*mol_mass(2))  
    
  end subroutine

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

 !--------------------------------------------------
 ! Utility functions
 !--------------------------------------------------

  subroutine multifab_weighted_sum(mfab, weights)
    type(multifab) , intent(inout) :: mfab(:,0:)
    real(dp_t), intent(in) :: weights(:)

    integer :: i,n,dm,rng
  
    !--------------------------------------
    ! Ghost cells need not be set here since syncs will set them later
    do n=1,size(mfab,1)
       call setval(mfab(n,0), 0.d0)
       ! We want to do mfab(n,0)=sum(weights(i)*mfab(n,i),i=1..n_rngs)
       do rng=1,size(weights)
         call saxpy(mfab(n,0), weights(rng), mfab(n,rng))
       end do
    end do   

  end subroutine multifab_weighted_sum
 
  subroutine multifab_fill_random(mfab, comp, variance, variance_mfab)
    type(multifab) , intent(inout) :: mfab(:)
    integer, intent(in), optional :: comp ! Only one component
    real(dp_t), intent(in), optional :: variance
    type(multifab) , intent(in), optional :: variance_mfab(:)

    integer :: i,n,dm,box
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

  ! Note: This routine assumes ghost values are valid on entry but does NOT sync ghost values at the end
  ! If ghost values need to remain consistent do a multifab_fill_boundary after calling multifab_filter!
  subroutine multifab_filter(mfab, dm)
    type(multifab) , intent(inout) :: mfab
    integer, intent(in) :: dm

    integer :: i,j,k,n,box
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


end module mk_stochastic_fluxdiv_module
